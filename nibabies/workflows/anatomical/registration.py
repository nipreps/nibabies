# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Within-baby registration of a T1w into a T2w image."""
from typing import Optional

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from pkg_resources import resource_filename as pkgr_fn


def init_coregistration_wf(
    *,
    bspline_fitting_distance: int = 200,
    mem_gb: float = 3.0,
    omp_nthreads: Optional[int] = None,
    sloppy: bool = False,
    debug: bool = False,
    precomputed_mask: bool = False,
    name: str = "coregistration_wf",
):
    """
    Set-up a T2w-to-T1w within-baby co-registration framework.

    See the ANTs' registration config file (under ``nibabies/data``) for further
    details.
    The main surprise in it is that, for some participants, accurate registration
    requires extra degrees of freedom (one affine level and one SyN level) to ensure
    that the T1w and T2w images align well.
    I attribute this requirement to the following potential reasons:

      * The T1w image and the T2w image were acquired in different sessions, apart in
        time enough for growth to happen.
        Although this is, in theory possible, it doesn't seem the images we have tested
        on are acquired on different sessions.
      * The skull is still so malleable that a change of position of the baby inside the
        coil made an actual change on the overall shape of their head.
      * Nonlinear distortions of the T1w and T2w images are, for some reason, more notorious
        for babies than they are for adults.
        We would need to look into each sequence's details to confirm this.

    Parameters
    ----------
    bspline_fitting_distance : :obj:`float`
        Distance in mm between B-Spline control points for N4 INU estimation.
    mem_gb : :obj:`float`
        Base memory fingerprint unit.
    name : :obj:`str`
        This particular workflow's unique name (Nipype requirement).
    omp_nthreads : :obj:`int`
        The number of threads for individual processes in this workflow.
    sloppy : :obj:`bool`
        Run in *sloppy* mode.
    debug : :obj:`bool`
        Produce intermediate registration files
    precomputed_mask : :obj:`bool`
        A precomputed mask for the T1w is available. In this case, generate a
        quick mask to assist in coregistration, but use the precomputed mask
        as the final output.


    Inputs
    ------
    in_t1w : :obj:`str`
        The preprocessed input T1w image (Denoising/INU/Clipping)
    in_t2w : :obj:`str`
        The preprocessed input T2w image (Denoising/INU/Clipping)
    in_mask : :obj:`str`
        The brainmask.
        If `precomputed_mask` is False, will be in T2w space.
        If `precomputed_mask` is True, will be in T1w space.
    in_probmap : :obj:`str`
        The probabilistic brainmask, as obtained in T2w space.

    Outputs
    -------
    t1w_coreg : :obj:`str`
        The preprocessed T1w image (INU and clipping), in its native space.
    t2w_coreg : :obj:`str`
        The preprocessed T2w image (INU and clipping), aligned into the T1w's space.
    t1w_brain : :obj:`str`
        The preprocessed, brain-extracted T1w image.
    t1w_mask : :obj:`str`
        The binary brainmask in T1w space.
    t1w2t2w_xfm : :obj:`str`
        The T1w-to-T2w mapping.

    """
    from nipype.interfaces.ants import N4BiasFieldCorrection
    from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
    from niworkflows.interfaces.fixes import FixHeaderRegistration as Registration
    from niworkflows.interfaces.nibabel import ApplyMask, Binarize, BinaryDilation

    from nibabies.utils.misc import get_file

    workflow = pe.Workflow(name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["in_t1w", "in_t2w", "in_mask", "in_probmap"]),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "t1w_preproc",
                "t1w_brain",
                "t1w_mask",
                "t1w2t2w_xfm",
                "t2w_preproc",
            ]
        ),
        name="outputnode",
    )

    # Dilate t2w mask for easier t1->t2 registration
    fixed_masks_arg = pe.Node(niu.Merge(3), name="fixed_masks_arg", run_without_submitting=True)
    reg_mask = pe.Node(BinaryDilation(radius=8, iterations=3), name="reg_mask")
    refine_mask = pe.Node(BinaryDilation(radius=8, iterations=1), name="refine_mask")

    # Set up T1w -> T2w within-subject registration
    coreg = pe.Node(
        Registration(from_file=get_file("nibabies", "data/within_subject_t1t2.json")),
        name="coreg",
        n_procs=omp_nthreads,
        mem_gb=mem_gb,
    )
    coreg.inputs.float = sloppy
    if debug:
        coreg.inputs.args = "--write-interval-volumes 5"
        coreg.inputs.output_inverse_warped_image = sloppy
        coreg.inputs.output_warped_image = sloppy

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
    # Move the T2w into T1w space, and apply the mask to the T1w
    map_t2w = pe.Node(ApplyTransforms(interpolation="BSpline"), name="map_t2w", mem_gb=1)
    apply_mask = pe.Node(ApplyMask(), name="apply_mask")

    # fmt: off
    workflow.connect([
        (inputnode, final_n4, [("in_t1w", "input_image")]),
        (inputnode, coreg, [("in_t1w", "moving_image"),
                            ("in_t2w", "fixed_image")]),
        (reg_mask, fixed_masks_arg, [
            ("out_file", "in1"),
            ("out_file", "in2")]),
        (refine_mask, fixed_masks_arg, [("out_file", "in3")]),
        (inputnode, map_t2w, [
            ("in_t1w", "reference_image"),
            ("in_t2w", "input_image")]),
        (fixed_masks_arg, coreg, [("out", "fixed_image_masks")]),
        (coreg, map_t2w, [
            ("reverse_transforms", "transforms"),
            ("reverse_invert_flags", "invert_transform_flags"),
        ]),
        (final_n4, apply_mask, [("output_image", "in_file")]),
        (final_n4, outputnode, [("output_image", "t1w_preproc")]),
        (map_t2w, outputnode, [("output_image", "t2w_preproc")]),
        (apply_mask, outputnode, [("out_file", "t1w_brain")]),
        (coreg, outputnode, [("forward_transforms", "t1w2t2w_xfm")]),
    ])
    # fmt: on

    if precomputed_mask:
        # The input mask is already in T1w space.
        # Generate a quick, rough mask of the T2w to be used to facilitate co-registration.
        from sdcflows.interfaces.brainmask import BrainExtraction

        masker = pe.Node(BrainExtraction(), name="t2w_masker")
        # fmt:off
        workflow.connect([
            (inputnode, masker, [("in_t2w", "in_file")]),
            (masker, reg_mask, [("out_mask", "in_file")]),
            (masker, refine_mask, [("out_mask", "in_file")]),
            (inputnode, apply_mask, [("in_mask", "in_mask")]),
            (inputnode, outputnode, [("in_mask", "t1w_mask")]),
        ])
        # fmt:on
    else:
        # The T2w mask from the brain extraction workflow will be mapped to T1w space
        map_mask = pe.Node(ApplyTransforms(interpolation="Gaussian"), name="map_mask", mem_gb=1)
        thr_mask = pe.Node(Binarize(thresh_low=0.80), name="thr_mask")
        # fmt:off
        workflow.connect([
            (inputnode, map_mask, [
                ("in_t1w", "reference_image"),
                ("in_probmap", "input_image")]),
            (inputnode, reg_mask, [("in_mask", "in_file")]),
            (inputnode, refine_mask, [("in_mask", "in_file")]),
            (coreg, map_mask, [
                ("reverse_transforms", "transforms"),
                ("reverse_invert_flags", "invert_transform_flags")]),
            (map_mask, thr_mask, [("output_image", "in_file")]),
            (map_mask, final_n4, [("output_image", "weight_image")]),
            (thr_mask, outputnode, [("out_mask", "t1w_mask")]),
            (thr_mask, apply_mask, [("out_mask", "in_mask")]),
        ])
        # fmt:on
    return workflow
