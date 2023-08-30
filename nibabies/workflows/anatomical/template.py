"""Prepare anatomical images for processing."""
from __future__ import annotations

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms


def init_anat_template_wf(
    *,
    contrast: str,
    num_files: int,
    omp_nthreads: int,
    longitudinal: bool = False,
    bspline_fitting_distance: int = 200,
    sloppy: bool = False,
    has_mask: bool = False,
    has_aseg: bool = False,
    name: str = "anat_template_wf",
) -> LiterateWorkflow:
    """
    Generate a canonically-oriented, structural average from all input images.
    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes
            from nibabies.workflows.anatomical.template import init_anat_template_wf
            wf = init_anat_template_wf(
                longitudinal=False, omp_nthreads=1, num_files=1, contrast="T1w"
            )
    Parameters
    ----------
    contrast : :obj:`str`
        Name of contrast, for reporting purposes, e.g., T1w, T2w, PDw
    num_files : :obj:`int`
        Number of images
    longitudinal : :obj:`bool`
        Create unbiased structural average, regardless of number of inputs
        (may increase runtime)
    omp_nthreads : :obj:`int`
        Maximum number of threads an individual process may use
    bspline_fitting_distance : :obj:`float`
        Distance in mm between B-Spline control points for N4 INU estimation.
    sloppy : :obj:`bool`
        Run in *sloppy* mode.
    name : :obj:`str`, optional
        Workflow name (default: anat_template_wf)
    Inputs
    ------
    anat_files
        List of structural images
    anat_mask
    mask_reference
    anat_aseg
    aseg_reference

    Outputs
    -------
    anat_ref
        Structural reference averaging input images
    anat_valid_list
        List of structural images accepted for combination
    anat_realign_xfm
        List of affine transforms to realign input images to final reference
    out_report
        Conformation report
    anat_mask
        Mask (if provided), resampled to the anatomical reference
    anat_aseg
        Aseg (if provided), resampled to the anatomical reference
    """
    from nipype.interfaces.ants import N4BiasFieldCorrection
    from nipype.interfaces.image import Reorient
    from niworkflows.interfaces.freesurfer import PatchedLTAConvert as LTAConvert
    from niworkflows.interfaces.freesurfer import StructuralReference
    from niworkflows.interfaces.images import Conform, TemplateDimensions
    from niworkflows.interfaces.nibabel import IntensityClip
    from niworkflows.interfaces.nitransforms import ConcatenateXFMs
    from niworkflows.utils.misc import add_suffix

    from nibabies.utils.misc import get_file

    wf = LiterateWorkflow(name=name)
    if num_files > 1:
        import nipype.interfaces.freesurfer as fs

        fs_ver = fs.Info().looseversion() or "<ver>"
        wf.__desc__ = f"""\
An anatomical {contrast}-reference map was computed after registration of
{num_files} {contrast} images (after INU-correction) using
`mri_robust_template` [FreeSurfer {fs_ver}, @fs_template].
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "anat_files",
                # Each derivative requires a reference file, which will be used to find which
                # transform to apply in the case when multiple runs are present
                "anat_mask",
                "mask_reference",
                "anat_aseg",
                "aseg_reference",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "anat_ref",
                "anat_valid_list",
                "anat_realign_xfm",
                "out_report",
                "anat_mask",
                "anat_aseg",
            ],
        ),
        name="outputnode",
    )

    # 0. Reorient anatomical image(s) to RAS and resample to common voxel space
    anat_ref_dimensions = pe.Node(TemplateDimensions(), name="anat_ref_dimensions")
    anat_conform = pe.MapNode(Conform(), iterfield="in_file", name="anat_conform")

    # fmt:off
    wf.connect([
        (inputnode, anat_ref_dimensions, [('anat_files', 't1w_list')]),
        (anat_ref_dimensions, anat_conform, [
            ('t1w_valid_list', 'in_file'),
            ('target_zooms', 'target_zooms'),
            ('target_shape', 'target_shape')]),
        (anat_ref_dimensions, outputnode, [('out_report', 'out_report'),
                                           ('t1w_valid_list', 'anat_valid_list')]),
    ])
    # fmt:on

    if has_mask:
        mask_conform = pe.Node(Conform(), name='mask_conform')
        # fmt:off
        wf.connect([
            (inputnode, mask_conform, [('anat_mask', 'in_file')]),
            (anat_ref_dimensions, mask_conform, [
                ('target_zooms', 'target_zooms'),
                ('target_shape', 'target_shape')]),
        ])
        # fmt:on

    if has_aseg:
        aseg_conform = pe.Node(Conform(), name='aseg_conform')
        # fmt:off
        wf.connect([
            (inputnode, aseg_conform, [('anat_aseg', 'in_file')]),
            (anat_ref_dimensions, aseg_conform, [
                ('target_zooms', 'target_zooms'),
                ('target_shape', 'target_shape')]),
        ])
        # fmt:on

    if num_files == 1:
        get1st = pe.Node(niu.Select(index=[0]), name="get1st")
        outputnode.inputs.anat_realign_xfm = [
            get_file("smriprep", "data/itkIdentityTransform.txt")
        ]

        # fmt:off
        wf.connect([
            (anat_conform, get1st, [('out_file', 'inlist')]),
            (get1st, outputnode, [('out', 'anat_ref')]),
        ])
        # fmt:on
        if has_mask:
            wf.connect(mask_conform, 'out_file', outputnode, 'anat_mask')
        if has_aseg:
            wf.connect(aseg_conform, 'out_file', outputnode, 'anat_aseg')
        return wf

    anat_conform_xfm = pe.MapNode(
        LTAConvert(in_lta="identity.nofile", out_lta=True),
        iterfield=["source_file", "target_file"],
        name="anat_conform_xfm",
    )
    # 1. Template (only if several anatomical images)
    # 1a. Correct for bias field: the bias field is an additive factor
    #     in log-transformed intensity units. Therefore, it is not a linear
    #     combination of fields and N4 fails with merged images.
    # 1b. Align and merge if several anatomical images are provided
    clip_preinu = pe.MapNode(IntensityClip(p_min=50), iterfield="in_file", name="clip_preinu")
    n4_correct = pe.MapNode(
        N4BiasFieldCorrection(
            dimension=3,
            save_bias=False,
            copy_header=True,
            n_iterations=[50] * (5 - 2 * sloppy),
            convergence_threshold=1e-7,
            shrink_factor=4,
            bspline_fitting_distance=bspline_fitting_distance,
        ),
        iterfield="input_image",
        name="n4_correct",
        n_procs=1,
    )
    # StructuralReference is fs.RobustTemplate if > 1 volume, copying otherwise
    anat_merge = pe.Node(
        StructuralReference(
            auto_detect_sensitivity=True,
            initial_timepoint=1,  # For deterministic behavior
            intensity_scaling=True,  # 7-DOF (rigid + intensity)
            subsample_threshold=200,
            fixed_timepoint=not longitudinal,
            no_iteration=not longitudinal,
            transform_outputs=True,
        ),
        mem_gb=2 * num_files - 1,
        name="anat_merge",
    )

    # 2. Reorient template to RAS, if needed (mri_robust_template may set to LIA)
    anat_reorient = pe.Node(Reorient(), name="anat_reorient")

    merge_xfm = pe.MapNode(
        niu.Merge(2),
        name="merge_xfm",
        iterfield=["in1", "in2"],
        run_without_submitting=True,
    )
    concat_xfms = pe.MapNode(
        ConcatenateXFMs(inverse=True),
        name="concat_xfms",
        iterfield=["in_xfms"],
        run_without_submitting=True,
    )

    if has_mask:
        mask_ref_idx = pe.Node(
            niu.Function(function=get_reference), name='mask_ref_idx', run_without_submitting=True
        )
        mask_xfm = pe.Node(niu.Select(), name='mask_xfm', run_without_submitting=True)
        applyxfm_mask = pe.Node(
            ApplyTransforms(interpolation='MultiLabel'), name='applyxfm_mask', mem_gb=1
        )
        mask_reorient = pe.Node(Reorient(), name="mask_reorient")
        # fmt:off
        wf.connect([
            (inputnode, mask_ref_idx, [('mask_reference', 'anat_reference')]),
            (anat_ref_dimensions, mask_ref_idx, [('t1w_valid_list', 'anatomicals')]),
            (concat_xfms, mask_xfm, [('out_xfm', 'inlist')]),
            (mask_ref_idx, mask_xfm, [('out', 'index')]),
            (mask_conform, applyxfm_mask, [('out_file', 'input_image')]),
            (anat_reorient, applyxfm_mask, [('out_file', 'reference_image')]),
            (mask_xfm, applyxfm_mask, [('out', 'transforms')]),
            (applyxfm_mask, mask_reorient, [('output_image', 'in_file')]),
            (mask_reorient, outputnode, [('out_file', 'anat_mask')]),
        ])
        # fmt:on

    if has_aseg:
        aseg_ref_idx = pe.Node(
            niu.Function(function=get_reference), name='aseg_ref_idx', run_without_submitting=True
        )
        aseg_xfm = pe.Node(niu.Select(), name='aseg_xfm', run_without_submitting=True)
        applyxfm_aseg = pe.Node(
            ApplyTransforms(interpolation='MultiLabel'), name='applyxfm_aseg', mem_gb=1
        )
        aseg_reorient = pe.Node(Reorient(), name="aseg_reorient")
        # fmt:off
        wf.connect([
            (inputnode, aseg_ref_idx, [('aseg_reference', 'anat_reference')]),
            (anat_ref_dimensions, aseg_ref_idx, [('t1w_valid_list', 'anatomicals')]),
            (concat_xfms, aseg_xfm, [('out_xfm', 'inlist')]),
            (aseg_ref_idx, aseg_xfm, [('out', 'index')]),
            (aseg_conform, applyxfm_aseg, [('out_file', 'input_image')]),
            (anat_reorient, applyxfm_aseg, [('out_file', 'reference_image')]),
            (aseg_xfm, applyxfm_aseg, [('out', 'transforms')]),
            (applyxfm_aseg, aseg_reorient, [('output_image', 'in_file')]),
            (applyxfm_aseg, outputnode, [('out_file', 'anat_aseg')]),
        ])
        # fmt:on

    def _set_threads(in_list, maximum):
        return min(len(in_list), maximum)

    # fmt:off
    wf.connect([
        (anat_ref_dimensions, anat_conform_xfm, [('t1w_valid_list', 'source_file')]),
        (anat_conform, anat_conform_xfm, [('out_file', 'target_file')]),
        (anat_conform, clip_preinu, [('out_file', 'in_file')]),
        (anat_conform, anat_merge, [
            (('out_file', _set_threads, omp_nthreads), 'num_threads'),
            (('out_file', add_suffix, '_template'), 'out_file')]),
        (clip_preinu, n4_correct, [('out_file', 'input_image')]),
        (n4_correct, anat_merge, [('output_image', 'in_files')]),
        (anat_merge, anat_reorient, [('out_file', 'in_file')]),
        # Combine orientation and template transforms
        (anat_conform_xfm, merge_xfm, [('out_lta', 'in1')]),
        (anat_merge, merge_xfm, [('transform_outputs', 'in2')]),
        (merge_xfm, concat_xfms, [('out', 'in_xfms')]),
        # Output
        (anat_reorient, outputnode, [('out_file', 'anat_ref')]),
        (concat_xfms, outputnode, [('out_xfm', 'anat_realign_xfm')]),
    ])
    # fmt:on
    return wf


def get_reference(anatomicals: list, anat_reference: str) -> int:
    return anatomicals.index(anat_reference)
