"""Prepare anatomical images for processing."""
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu


def init_anat_average_wf(
    *,
    bspline_fitting_distance=200,
    longitudinal=False,
    omp_nthreads=None,
    num_maps=1,
    name="anat_average_wf",
    sloppy=False,
):
    """
    Adapts :py:func:`~smriprep.workflows.anatomical.init_anat_template_wf` for T2w image reference
    """
    from pkg_resources import resource_filename as pkgr
    from nipype.interfaces.ants import N4BiasFieldCorrection
    from nipype.interfaces.image import Reorient

    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.header import ValidateImage
    from niworkflows.interfaces.nibabel import IntensityClip, SplitSeries
    from niworkflows.interfaces.freesurfer import (
        StructuralReference,
        PatchedLTAConvert as LTAConvert,
    )
    from niworkflows.interfaces.images import TemplateDimensions, Conform
    from niworkflows.interfaces.nitransforms import ConcatenateXFMs
    from niworkflows.utils.misc import add_suffix

    wf = Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=["in_files"]), name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=["out_file", "valid_list", "realign_xfms", "out_report"]
        ),
        name="outputnode",
    )

    # 1. Validate each of the input images
    validate = pe.MapNode(
        ValidateImage(),
        iterfield="in_file",
        name="validate",
        run_without_submitting=True,
    )

    # 2. Ensure we don't have two timepoints and implicitly squeeze image
    split = pe.MapNode(SplitSeries(), iterfield="in_file", name="split")

    # 3. INU correction of all independent volumes
    clip_preinu = pe.MapNode(
        IntensityClip(p_min=50), iterfield="in_file", name="clip_preinu"
    )
    correct_inu = pe.MapNode(
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
        n_procs=omp_nthreads,
        name="correct_inu",
    )
    clip_postinu = pe.MapNode(
        IntensityClip(p_min=10.0, p_max=99.5), iterfield="in_file", name="clip_postinu"
    )

    # 4. Reorient T2w image(s) to RAS and resample to common voxel space
    ref_dimensions = pe.Node(TemplateDimensions(), name="ref_dimensions")
    conform = pe.MapNode(Conform(), iterfield="in_file", name="conform")
    # fmt:off
    wf.connect([
        (inputnode, validate, [("in_files", "in_file")]),
        (validate, split, [("out_file", "in_file")]),
        (split, clip_preinu, [(("out_files", _flatten), "in_file")]),
        (clip_preinu, correct_inu, [("out_file", "input_image")]),
        (correct_inu, clip_postinu, [("output_image", "in_file")]),
        (clip_postinu, ref_dimensions, [("out_file", "t1w_list")]),
        (ref_dimensions, conform, [
            ("t1w_valid_list", "in_file"),
            ("target_zooms", "target_zooms"),
            ("target_shape", "target_shape")]),
        (ref_dimensions, outputnode, [("out_report", "out_report"),
                                      ("t1w_valid_list", "valid_list")]),
    ])
    # fmt:on

    # 5. Reorient template to RAS, if needed (mri_robust_template may set to LIA)
    ensure_ras = pe.Node(Reorient(), name="ensure_ras")

    if num_maps == 1:
        get1st = pe.Node(niu.Select(index=[0]), name="get1st")
        outputnode.inputs.realign_xfms = [
            pkgr("smriprep", "data/itkIdentityTransform.txt")
        ]
        # fmt:off
        wf.connect([
            (conform, get1st, [("out_file", "inlist")]),
            (get1st, ensure_ras, [("out", "in_file")]),
            (ensure_ras, outputnode, [("out_file", "out_file")]),
        ])
        # fmt:on
        return wf

    from nipype.interfaces import freesurfer as fs

    wf.__desc__ = f"""\
An anatomical reference-map was computed after registration of
{num_maps} images (after INU-correction) using
`mri_robust_template` [FreeSurfer {fs.Info().looseversion() or "<ver>"}, @fs_template].
"""

    conform_xfm = pe.MapNode(
        LTAConvert(in_lta="identity.nofile", out_lta=True),
        iterfield=["source_file", "target_file"],
        name="conform_xfm",
    )

    # 6. StructuralReference is fs.RobustTemplate if > 1 volume, copying otherwise
    merge = pe.Node(
        StructuralReference(
            auto_detect_sensitivity=True,
            initial_timepoint=1,  # For deterministic behavior
            intensity_scaling=True,  # 7-DOF (rigid + intensity)
            subsample_threshold=200,
            fixed_timepoint=not longitudinal,
            no_iteration=not longitudinal,
            transform_outputs=True,
        ),
        mem_gb=2 * num_maps - 1,
        name="merge",
    )

    # 7. Final intensity equalization/conformation
    clip_final = pe.Node(IntensityClip(p_min=2.0, p_max=99.9), name="clip_final")

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
        (ref_dimensions, conform_xfm, [("t1w_valid_list", "source_file")]),
        (conform, conform_xfm, [("out_file", "target_file")]),
        (conform, merge, [
            ("out_file", "in_file"),
            (("out_file", _set_threads, omp_nthreads), "num_threads"),
            (("out_file", add_suffix, "_template"), "out_file")]),
        (merge, ensure_ras, [("out_file", "in_file")]),
        # Combine orientation and template transforms
        (conform_xfm, merge_xfm, [("out_lta", "in1")]),
        (merge, merge_xfm, [("transform_outputs", "in2")]),
        (merge_xfm, concat_xfms, [("out", "in_xfms")]),
        # Output
        (ensure_ras, clip_final, [("out_file", "in_file")]),
        (clip_final, outputnode, [("out_file", "out_file")]),
        (concat_xfms, outputnode, [("out_xfm", "realign_xfms")]),
    ])
    # fmt:on

    return wf


def _flatten(inlist):
    from bids.utils import listify

    return [el for items in listify(inlist) for el in listify(items)]
