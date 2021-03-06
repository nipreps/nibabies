# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Nipype translation of ANTs' workflows."""
# import numpy as np
# general purpose
from pkg_resources import resource_filename as pkgr_fn

# nipype
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu


def init_coregistration_wf(
    bspline_fitting_distance=200,
    mem_gb=3.0,
    name="coregistration_wf",
    omp_nthreads=None,
    sloppy=False,
):
    """
    Set-up a T2w-to-T1w within-baby co-registration framework.

    Parameters
    ----------
    mem_gb : :obj:`float`
        Base memory fingerprint unit.
    name : :obj:`str`
        This particular workflow's unique name (Nipype requirement).
    omp_nthreads : :obj:`int`
        The number of threads for individual processes in this workflow.
    sloppy : :obj:`bool`
        Run in *sloppy* mode.

    Inputs
    ------
    in_t1w : :obj:`str`
        The unprocessed input T1w image.
    in_t2w_preproc : :obj:`str`
        The preprocessed input T2w image, from the brain extraction workflow.
    in_mask : :obj:`str`
        The brainmask, as obtained in T2w space.
    in_probmap : :obj:`str`
        The probabilistic brainmask, as obtained in T2w space.

    Outputs
    -------
    t1w_preproc : :obj:`str`
        The preprocessed T1w image (INU and clipping).
    t1w_brain : :obj:`str`
        The preprocessed, brain-extracted T1w image.
    t1w_mask : :obj:`str`
        The binary brainmask projected from the T2w.
    t1w2t2w_xfm : :obj:`str`
        The T1w-to-T2w mapping.

    """
    from nipype.interfaces.ants import N4BiasFieldCorrection
    from niworkflows.interfaces.fixes import (
        FixHeaderRegistration as Registration,
        FixHeaderApplyTransforms as ApplyTransforms,
    )
    from niworkflows.interfaces.nibabel import ApplyMask, Binarize

    from ...interfaces.nibabel import IntensityClip

    workflow = pe.Workflow(name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=["in_t1w", "in_t2w_preproc", "in_mask", "in_probmap"]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=["t1w_preproc", "t1w_brain", "t1w_mask", "t1w2t2w_xfm"]
        ),
        name="outputnode",
    )

    pre_n4_clip = pe.Node(IntensityClip(), name="pre_n4_clip")
    init_n4 = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            save_bias=False,
            copy_header=True,
            n_iterations=[50] * (4 - sloppy),
            convergence_threshold=1e-7,
            shrink_factor=4,
            bspline_fitting_distance=bspline_fitting_distance,
        ),
        n_procs=omp_nthreads,
        name="init_n4",
    )
    post_n4_clip = pe.Node(IntensityClip(p_min=2.0, p_max=100), name="post_n4_clip")

    fixed_masks_arg = pe.Node(
        niu.Merge(4), name="fixed_masks_arg", run_without_submitting=True
    )
    fixed_masks_arg.inputs.in1 = "NULL"
    fixed_masks_arg.inputs.in2 = "NULL"
    fixed_masks_arg.inputs.in3 = "NULL"

    # dilate t2w mask for easier t1->t2 registration
    # Unclear whether this is necessary
    # dil_brainmask = pe.Node(
    #     ImageMath(operation="MD", op2="8", copy_header=True), name="dil_brainmask"
    # )

    # Set up T2w -> T1w within-subject registration
    coreg = pe.Node(
        Registration(from_file=pkgr_fn("nibabies.data", "within_subject_t1t2.json")),
        name="coreg",
        n_procs=omp_nthreads,
        mem_gb=mem_gb,
    )
    coreg.inputs.float = sloppy

    map_mask = pe.Node(
        ApplyTransforms(interpolation="Gaussian"), name="map_mask", mem_gb=1
    )
    thr_mask = pe.Node(Binarize(thresh_low=0.80), name="thr_mask")

    final_n4 = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            bspline_fitting_distance=bspline_fitting_distance,
            save_bias=True,
            copy_header=True,
            n_iterations=[50] * 5,
            convergence_threshold=1e-7,
            rescale_intensities=True,
            shrink_factor=4,
        ),
        n_procs=omp_nthreads,
        name="final_n4",
    )
    apply_mask = pe.Node(ApplyMask(), name="apply_mask")

    # fmt:off
    workflow.connect([
        (inputnode, map_mask, [("in_t1w", "reference_image")]),
        (inputnode, pre_n4_clip, [("in_t1w", "in_file")]),
        (inputnode, coreg, [("in_t2w_preproc", "fixed_image")]),
        (inputnode, map_mask, [("in_probmap", "input_image")]),
        (inputnode, fixed_masks_arg, [("in_mask", "in4")]),
        (pre_n4_clip, init_n4, [("out_file", "input_image")]),
        (init_n4, post_n4_clip, [("output_image", "in_file")]),
        (post_n4_clip, coreg, [("out_file", "moving_image")]),
        (fixed_masks_arg, coreg, [("out", "fixed_image_masks")]),
        (coreg, map_mask, [
            ("reverse_transforms", "transforms"),
            ("reverse_invert_flags", "invert_transform_flags"),
        ]),
        (map_mask, thr_mask, [("output_image", "in_file")]),
        (pre_n4_clip, final_n4, [("out_file", "input_image")]),
        (map_mask, final_n4, [("output_image", "weight_image")]),
        (final_n4, apply_mask, [("output_image", "in_file")]),
        (thr_mask, apply_mask, [("out_mask", "in_mask")]),
        (final_n4, outputnode, [("output_image", "t1w_preproc")]),
        (thr_mask, outputnode, [("out_mask", "t1w_mask")]),
        (apply_mask, outputnode, [("out_file", "t1w_corrected_brain")]),
        (coreg, map_mask, [("forward_transforms", "t1w2t2w_xfm")]),
    ])
    # fmt:on
    return workflow
