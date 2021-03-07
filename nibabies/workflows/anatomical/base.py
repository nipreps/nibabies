from pathlib import Path

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl


def init_infant_anat_wf(
    *,
    age_months,
    ants_affine_init,
    t1w,
    t2w,
    anat_modality,
    bids_root,
    existing_derivatives,
    freesurfer,
    longitudinal,
    omp_nthreads,
    output_dir,
    segmentation_atlases,
    skull_strip_mode,
    skull_strip_template,
    sloppy,
    spaces,
    name="infant_anat_wf",
):
    """

      - T1w reference: realigning and then averaging anatomical images.
      - Brain extraction and INU (bias field) correction.
      - Brain tissue segmentation.
      - Spatial normalization to standard spaces.
      - Surface reconstruction with FreeSurfer_.

    Outputs
    -------

    anat_preproc
        The anatomical reference map, which is calculated as the average of bias-corrected
        and preprocessed anatomical images, defining the anatomical space.
    anat_brain
        Skull-stripped ``anat_preproc``
    anat_mask
        Brain (binary) mask estimated by brain extraction.
    anat_dseg
        Brain tissue segmentation of the preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF).
    anat_tpms
        List of tissue probability maps corresponding to ``t1w_dseg``.
    std_preproc
        T1w reference resampled in one or more standard spaces.
    std_mask
        Mask of skull-stripped template, in MNI space
    std_dseg
        Segmentation, resampled into MNI space
    std_tpms
        List of tissue probability maps in MNI space
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    anat2std_xfm
        Nonlinear spatial transform to resample imaging data given in anatomical space
        into standard space.
    std2anat_xfm
        Inverse transform of the above.
    subject_id
        FreeSurfer subject ID
    anat2fsnative_xfm
        LTA-style affine matrix translating from T1w to
        FreeSurfer-conformed subject space
    fsnative2anat_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed
        subject space to T1w
    surfaces
        GIFTI surfaces (gray/white boundary, midthickness, pial, inflated)
    """
    from nipype.interfaces.base import Undefined
    from nipype.interfaces.ants.base import Info as ANTsInfo
    from niworkflows.interfaces.header import ValidateImage
    from smriprep.workflows.anatomical import init_anat_template_wf, _pop
    from smriprep.workflows.norm import init_anat_norm_wf
    from smriprep.workflows.outputs import (
        init_anat_reports_wf,
        init_anat_derivatives_wf,
    )

    from ...utils.misc import fix_multi_source_name
    from .brain_extraction import init_infant_brain_extraction_wf
    from .registration import init_coregistration_wf
    from .segmentation import init_anat_seg_wf
    from .surfaces import init_infant_surface_recon_wf
    from .outputs import init_coreg_report_wf

    # for now, T1w only
    num_t1w = len(t1w) if t1w else 0
    num_t2w = len(t2w) if t2w else 0

    wf = pe.Workflow(name=name)
    desc = """Anatomical data preprocessing

: """
    desc += f"""\
A total of {num_t1w} T1w and {num_t2w} T2w images were found within the input
BIDS dataset."""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=["t1w", "t2w", "subject_id", "subjects_dir"]
        ),  # FLAIR / ROI?
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "anat_preproc",
                "anat_brain",
                "anat_mask",
                "anat_dseg",
                "anat_tpms",
                "anat_ref_xfms",
                "std_preproc",
                "std_brain",
                "std_dseg",
                "std_tpms",
                "subjects_dir",
                "subject_id",
                "anat2std_xfm",
                "std2anat_xfm",
                "anat2fsnative_xfm",
                "fsnative2anat_xfm",
                "surfaces",
                "anat_aseg",
                "anat_aparc",
            ]
        ),
        name="outputnode",
    )

    # Connect reportlets workflows
    anat_reports_wf = init_anat_reports_wf(
        freesurfer=freesurfer,
        output_dir=output_dir,
    )

    if existing_derivatives:
        raise NotImplementedError("Reusing derivatives is not yet supported.")

    desc += (
        """
All of the T1-weighted images were corrected for intensity non-uniformity (INU)
"""
        if num_t1w > 1
        else """\
The T1-weighted (T1w) image was corrected for intensity non-uniformity (INU)
"""
    )
    desc += """\
with `N4BiasFieldCorrection` [@n4], distributed with ANTs {ants_ver} \
[@ants, RRID:SCR_004757]"""
    desc += (
        ".\n"
        if num_t1w > 1
        else ", and used as T1w-reference throughout the workflow.\n"
    )

    desc += """\
The T1w-reference was then skull-stripped with a modified implementation of
the `antsBrainExtraction.sh` workflow (from ANTs), using {skullstrip_tpl}
as target template.
Brain tissue segmentation of cerebrospinal fluid (CSF),
white-matter (WM) and gray-matter (GM) was performed on
the brain-extracted T1w using ANTs JointFusion, distributed with ANTs {ants_ver}.
"""

    wf.__desc__ = desc.format(
        ants_ver=ANTsInfo.version() or "(version unknown)",
        skullstrip_tpl=skull_strip_template.fullname,
    )
    # Define output workflows
    anat_reports_wf = init_anat_reports_wf(freesurfer=freesurfer, output_dir=output_dir)
    # HACK: remove resolution from TFSelect
    anat_reports_wf.get_node("tf_select").inputs.resolution = Undefined

    anat_derivatives_wf = init_anat_derivatives_wf(
        bids_root=bids_root,
        freesurfer=freesurfer,
        num_t1w=num_t1w,
        output_dir=output_dir,
        spaces=spaces,
    )
    # HACK: remove resolution from TFSelect
    anat_derivatives_wf.get_node("select_tpl").inputs.resolution = Undefined

    # Multiple T1w files -> generate average reference
    t1w_template_wf = init_anat_template_wf(
        longitudinal=False,
        omp_nthreads=omp_nthreads,
        num_t1w=num_t1w,
    )

    t2w_template_wf = init_t2w_template_wf(
        longitudinal=longitudinal,
        omp_nthreads=omp_nthreads,
        num_t2w=num_t2w,
    )

    anat_validate = pe.Node(
        ValidateImage(),
        name="anat_validate",
        run_without_submitting=True,
    )

    # INU + Brain Extraction
    if skull_strip_mode != "force":
        raise NotImplementedError("Skull stripping is currently required.")

    brain_extraction_wf = init_infant_brain_extraction_wf(
        age_months=age_months,
        ants_affine_init=ants_affine_init,
        skull_strip_template=skull_strip_template.space,
        template_specs=skull_strip_template.spec,
        omp_nthreads=omp_nthreads,
        sloppy=sloppy,
    )

    coregistration_wf = init_coregistration_wf(
        omp_nthreads=omp_nthreads,
        sloppy=sloppy,
    )
    coreg_report_wf = init_coreg_report_wf(
        output_dir=output_dir,
    )

    # Ensure single outputs
    be_buffer = pe.Node(
        niu.IdentityInterface(fields=["anat_preproc", "anat_brain"]), name="be_buffer"
    )

    # Segmentation - initial implementation should be simple: JLF
    anat_seg_wf = init_anat_seg_wf(
        age_months=age_months,
        anat_modality=anat_modality.capitalize(),
        template_dir=segmentation_atlases,
        sloppy=sloppy,
        omp_nthreads=omp_nthreads,
    )

    # Spatial normalization (requires segmentation)
    anat_norm_wf = init_anat_norm_wf(
        debug=sloppy,
        omp_nthreads=omp_nthreads,
        templates=spaces.get_spaces(nonstandard=False, dim=(3,)),
    )
    # HACK: remove resolution from TFSelect
    anat_norm_wf.get_node("tf_select").inputs.resolution = Undefined
    # HACK: requires patched niworkflows to allow setting resolution to none
    anat_norm_wf.get_node("registration").inputs.template_resolution = None

    # fmt: off
    wf.connect([
        (inputnode, t1w_template_wf, [
            ("t1w", "inputnode.t1w"),
        ]),
        (inputnode, t2w_template_wf, [
            ("t2w", "inputnode.t2w"),
        ]),
        (t1w_template_wf, outputnode, [
            ("outputnode.t1w_realign_xfm", "anat_ref_xfms"),
        ]),
        (t1w_template_wf, anat_validate, [
            ("outputnode.t1w_ref", "in_file"),
        ]),
        (anat_validate, coregistration_wf, [
            ("out_file", "inputnode.in_t1w"),
        ]),
        (t2w_template_wf, brain_extraction_wf, [
            ("outputnode.t2w_ref", "inputnode.in_t2w"),
        ]),
        (brain_extraction_wf, coregistration_wf, [
            ("outputnode.t2w_preproc", "inputnode.in_t2w_preproc"),
            ("outputnode.out_mask", "inputnode.in_mask"),
            ("outputnode.out_probmap", "inputnode.in_probmap"),
        ]),
        (coregistration_wf, be_buffer, [
            (("outputnode.t1w_preproc", _pop), "anat_preproc"),
            (("outputnode.t1w_brain", _pop), "anat_brain"),
            (("outputnode.t1w_mask", _pop), "anat_mask"),
        ]),
        (be_buffer, outputnode, [
            ("anat_preproc", "anat_preproc"),
            ("anat_brain", "anat_brain"),
            ("anat_mask", "anat_mask"),
        ]),
        (be_buffer, anat_seg_wf, [
            ("anat_brain", "inputnode.anat_brain"),
        ]),
        (anat_seg_wf, outputnode, [
            ("outputnode.anat_dseg", "anat_dseg"),
        ]),
        (anat_seg_wf, anat_norm_wf, [
            ("outputnode.anat_dseg", "inputnode.moving_segmentation"),
            ("outputnode.anat_tpms", "inputnode.moving_tpms"),
        ]),
        (be_buffer, anat_norm_wf, [
            ("anat_preproc", "inputnode.moving_image"),
            ("anat_mask", "inputnode.moving_mask"),
        ]),
        (anat_norm_wf, outputnode, [
            ("poutputnode.standardized", "std_preproc"),
            ("poutputnode.std_mask", "std_mask"),
            ("poutputnode.std_dseg", "std_dseg"),
            ("poutputnode.std_tpms", "std_tpms"),
            ("outputnode.template", "template"),
            ("outputnode.anat2std_xfm", "anat2std_xfm"),
            ("outputnode.std2anat_xfm", "std2anat_xfm"),
        ]),
        (inputnode, anat_norm_wf, [
            (("t1w", fix_multi_source_name), "inputnode.orig_t1w"),  # anat_validate? not used...
        ]),
    ])

    wf.connect([
        # reports
        (inputnode, anat_reports_wf, [
            ("t1w", "inputnode.source_file"),
        ]),
        (inputnode, coreg_report_wf, [
            ("t1w", "inputnode.source_file"),
        ]),
        (coregistration_wf, coreg_report_wf, [
            ("outputnode.t1w_preproc", "inputnode.t1w_preproc"),
            ("outputnode.t2w_preproc", "inputnode.t2w_preproc"),
            ("outputnode.t1w_mask", "inputnode.in_mask"),
        ]),
        (outputnode, anat_reports_wf, [
            ("anat_preproc", "inputnode.t1w_preproc"),
            ("anat_mask", "inputnode.t1w_mask"),
            ("anat_dseg", "inputnode.t1w_dseg"),
            ("std_preproc", "inputnode.std_t1w"),
            ("std_mask", "inputnode.std_mask"),
        ]),
        (t1w_template_wf, anat_reports_wf, [
            ("outputnode.out_report", "inputnode.t1w_conform_report"),
        ]),
        (anat_norm_wf, anat_reports_wf, [
            ("poutputnode.template", "inputnode.template"),
        ]),
        # derivatives
        (t1w_template_wf, anat_derivatives_wf, [
            ("outputnode.t1w_valid_list", "inputnode.source_files"),
            ("outputnode.t1w_realign_xfm", "inputnode.t1w_ref_xfms"),
        ]),
        (be_buffer, anat_derivatives_wf, [
            ("anat_mask", "inputnode.t1w_mask"),
            ("anat_preproc", "inputnode.t1w_preproc"),
        ]),
        (anat_norm_wf, anat_derivatives_wf, [
            ("outputnode.template", "inputnode.template"),
            ("outputnode.anat2std_xfm", "inputnode.anat2std_xfm"),
            ("outputnode.std2anat_xfm", "inputnode.std2anat_xfm"),
        ]),
        (anat_seg_wf, anat_derivatives_wf, [
            ("outputnode.anat_dseg", "inputnode.t1w_dseg"),
            ("outputnode.anat_tpms", "inputnode.t1w_tpms"),
        ]),
    ])
    # fmt:on

    if not freesurfer:
        return wf

    # FreeSurfer surfaces
    surface_recon_wf = init_infant_surface_recon_wf(
        age_months=age_months,
        use_aseg=bool(segmentation_atlases),
    )

    # fmt:off
    wf.connect([
        (anat_seg_wf, surface_recon_wf, [
            ("outputnode.anat_aseg", "inputnode.anat_aseg"),
        ]),
        (inputnode, surface_recon_wf, [
            ("subject_id", "inputnode.subject_id"),
            ("subjects_dir", "inputnode.subjects_dir"),
            ("t2w", "inputnode.t2w"),
        ]),
        (anat_validate, surface_recon_wf, [
            ("out_file", "inputnode.anat_orig"),
        ]),
        (be_buffer, surface_recon_wf, [
            ("anat_brain", "inputnode.anat_skullstripped"),
            ("anat_preproc", "inputnode.anat_preproc"),
        ]),
        (surface_recon_wf, outputnode, [
            ("outputnode.subjects_dir", "subjects_dir"),
            ("outputnode.subject_id", "subject_id"),
            ("outputnode.anat2fsnative_xfm", "anat2fsnative_xfm"),
            ("outputnode.fsnative2anat_xfm", "fsnative2anat_xfm"),
            ("outputnode.surfaces", "surfaces"),
            ("outputnode.anat_aparc", "anat_aparc"),
            ("outputnode.anat_aseg", "anat_aseg"),
        ]),
        (surface_recon_wf, anat_reports_wf, [
            ("outputnode.subject_id", "inputnode.subject_id"),
            ("outputnode.subjects_dir", "inputnode.subjects_dir"),
        ]),
        (surface_recon_wf, anat_derivatives_wf, [
            ("outputnode.anat_aseg", "inputnode.t1w_fs_aseg"),
            ("outputnode.anat_aparc", "inputnode.t1w_fs_aparc"),
            ("outputnode.anat2fsnative_xfm", "inputnode.t1w2fsnative_xfm"),
            ("outputnode.fsnative2anat_xfm", "inputnode.fsnative2t1w_xfm"),
            ("outputnode.surfaces", "inputnode.surfaces"),
        ]),
    ])
    # fmt: on
    return wf


def init_t2w_template_wf(
    longitudinal, omp_nthreads, num_t2w, name="anat_t2w_template_wf"
):
    """
    Adapts :py:func:`~smriprep.workflows.anatomical.init_anat_template_wf` for T2w image reference
    """
    from pkg_resources import resource_filename as pkgr
    from nipype.interfaces import freesurfer as fs, image, ants
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.freesurfer import (
        StructuralReference,
        PatchedLTAConvert as LTAConvert,
    )
    from niworkflows.interfaces.images import TemplateDimensions, Conform
    from niworkflows.interfaces.nitransforms import ConcatenateXFMs
    from niworkflows.utils.misc import add_suffix

    wf = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=["t2w"]), name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=["t2w_ref", "t2w_valid_list", "t2_realign_xfm", "out_report"]
        ),
        name="outputnode",
    )

    # 0. Reorient T2w image(s) to RAS and resample to common voxel space
    t2w_ref_dimensions = pe.Node(TemplateDimensions(), name="t2w_ref_dimensions")
    t2w_conform = pe.MapNode(Conform(), iterfield="in_file", name="t2w_conform")
    # fmt:off
    wf.connect([
        (inputnode, t2w_ref_dimensions, [('t2w', 't1w_list')]),
        (t2w_ref_dimensions, t2w_conform, [
            ('t1w_valid_list', 'in_file'),
            ('target_zooms', 'target_zooms'),
            ('target_shape', 'target_shape')]),
        (t2w_ref_dimensions, outputnode, [('out_report', 'out_report'),
                                          ('t1w_valid_list', 't2w_valid_list')]),
    ])
    # fmt:on

    if num_t2w == 1:
        get1st = pe.Node(niu.Select(index=[0]), name="get1st")
        outputnode.inputs.t2w_realign_xfm = [
            pkgr("smriprep", "data/itkIdentityTransform.txt")
        ]
        # fmt:off
        wf.connect([
            (t2w_conform, get1st, [('out_file', 'inlist')]),
            (get1st, outputnode, [('out', 't2w_ref')]),
        ])
        # fmt:on
        return wf

    wf.__desc__ = f"""\
A T2w-reference map was computed after registration of
{num_t2w} T2w images (after INU-correction) using
`mri_robust_template` [FreeSurfer {fs.Info().looseversion() or "<ver>"}, @fs_template].
"""

    t2w_conform_xfm = pe.MapNode(
        LTAConvert(in_lta="identity.nofile", out_lta=True),
        iterfield=["source_file", "target_file"],
        name="t2w_conform_xfm",
    )

    # 1a. Correct for bias field: the bias field is an additive factor
    #     in log-transformed intensity units. Therefore, it is not a linear
    #     combination of fields and N4 fails with merged images.
    # 1b. Align and merge if several T1w images are provided
    n4_correct = pe.MapNode(
        ants.N4BiasFieldCorrection(dimension=3, copy_header=True),
        iterfield="input_image",
        name="n4_correct",
        n_procs=1,
    )  # n_procs=1 for reproducibility

    # StructuralReference is fs.RobustTemplate if > 1 volume, copying otherwise
    t2w_merge = pe.Node(
        StructuralReference(
            auto_detect_sensitivity=True,
            initial_timepoint=1,  # For deterministic behavior
            intensity_scaling=True,  # 7-DOF (rigid + intensity)
            subsample_threshold=200,
            fixed_timepoint=not longitudinal,
            no_iteration=not longitudinal,
            transform_outputs=True,
        ),
        mem_gb=2 * num_t2w - 1,
        name="t2w_merge",
    )

    # 2. Reorient template to RAS, if needed (mri_robust_template may set to LIA)
    t2w_reorient = pe.Node(image.Reorient(), name="t2w_reorient")

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

    def _set_threads(in_list, maximum):
        return min(len(in_list), maximum)

    # fmt:off
    wf.connect([
        (t2w_ref_dimensions, t2w_conform_xfm, [('t1w_valid_list', 'source_file')]),
        (t2w_conform, t2w_conform_xfm, [('out_file', 'target_file')]),
        (t2w_conform, n4_correct, [('out_file', 'input_image')]),
        (t2w_conform, t2w_merge, [
            (('out_file', _set_threads, omp_nthreads), 'num_threads'),
            (('out_file', add_suffix, '_template'), 'out_file')]),
        (n4_correct, t2w_merge, [('output_image', 'in_files')]),
        (t2w_merge, t2w_reorient, [('out_file', 'in_file')]),
        # Combine orientation and template transforms
        (t2w_conform_xfm, merge_xfm, [('out_lta', 'in1')]),
        (t2w_merge, merge_xfm, [('transform_outputs', 'in2')]),
        (merge_xfm, concat_xfms, [('out', 'in_xfms')]),
        # Output
        (t2w_reorient, outputnode, [('out_file', 't2w_ref')]),
        (concat_xfms, outputnode, [('out_xfm', 't2w_realign_xfm')]),
    ])
    # fmt:on

    return wf
