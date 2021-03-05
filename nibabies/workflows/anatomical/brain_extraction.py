# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Nipype translation of ANTs' workflows."""
# import numpy as np
# general purpose
from pkg_resources import resource_filename as pkgr_fn

# nipype
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu


def init_infant_brain_extraction_wf(
    age_months=None,
    ants_affine_init=False,
    bspline_fitting_distance=200,
    sloppy=False,
    skull_strip_template="UNCInfant",
    template_specs=None,
    mem_gb=3.0,
    name="infant_brain_extraction_wf",
    omp_nthreads=None,
):
    """
    Build an atlas-based brain extraction pipeline for infant T1w/T2w MRI data.

    Pros/Cons of available templates
    --------------------------------
    * MNIInfant
     + More cohorts available for finer-grain control
     + T1w/T2w images available
     - Template masks are poor

    * UNCInfant
     + Accurate masks
     - No T2w image available


    Parameters
    ----------
    ants_affine_init : :obj:`bool`, optional
        Set-up a pre-initialization step with ``antsAI`` to account for mis-oriented images.

    """
    from nipype.interfaces.ants import N4BiasFieldCorrection, ImageMath

    # niworkflows
    from niworkflows.interfaces.header import ValidateImage
    from niworkflows.interfaces.nibabel import ApplyMask, Binarize
    from niworkflows.interfaces.fixes import (
        FixHeaderRegistration as Registration,
        FixHeaderApplyTransforms as ApplyTransforms,
    )
    from templateflow.api import get as get_template

    from ...interfaces.nibabel import IntensityClip
    from ...utils.misc import cohort_by_months

    # handle template specifics
    template_specs = template_specs or {}
    if skull_strip_template == "MNIInfant":
        template_specs["resolution"] = 2 if sloppy else 1

    if not template_specs.get("cohort"):
        if age_months is None:
            raise KeyError(
                f"Age or cohort for {skull_strip_template} must be provided!"
            )
        template_specs["cohort"] = cohort_by_months(skull_strip_template, age_months)

    tpl_target_path = get_template(
        skull_strip_template,
        suffix="T1w",  # no T2w template
        desc=None,
        **template_specs,
    )
    if not tpl_target_path:
        raise RuntimeError(
            f"An instance of template <tpl-{skull_strip_template}> with T1w suffix "
            "could not be found."
        )

    tpl_brainmask_path = get_template(
        skull_strip_template, label="brain", suffix="probseg", **template_specs
    ) or get_template(
        skull_strip_template, desc="brain", suffix="mask", **template_specs
    )

    tpl_regmask_path = get_template(
        skull_strip_template,
        label="BrainCerebellumExtraction",
        suffix="mask",
        **template_specs,
    )

    # main workflow
    workflow = pe.Workflow(name)

    inputnode = pe.Node(niu.IdentityInterface(fields=["in_t2w"]), name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=["t2w_preproc", "t2w_brain", "out_mask", "out_probmap"]
        ),
        name="outputnode",
    )

    # validate images
    val_t2w = pe.Node(ValidateImage(), name="val_t2w")

    # truncate target intensity for N4 correction
    clip_tmpl = pe.Node(IntensityClip(in_file=_pop(tpl_target_path)), name="clip_tmpl")
    clip_t2w = pe.Node(IntensityClip(), name="clip_t2w")

    # INU correction of the t1w
    init_t2w_n4 = pe.Node(
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
        name="init_t2w_n4",
    )

    clip_t2w_inu = pe.Node(IntensityClip(p_min=2.0, p_max=100.0), name="clip_t2w_inu")

    # Spatial normalization step
    lap_tmpl = pe.Node(ImageMath(operation="Laplacian", op2="0.4 1"), name="lap_tmpl")
    lap_t2w = pe.Node(ImageMath(operation="Laplacian", op2="0.4 1"), name="lap_t2w")

    # Merge image nodes
    mrg_tmpl = pe.Node(niu.Merge(2), name="mrg_tmpl", run_without_submitting=True)
    mrg_t2w = pe.Node(niu.Merge(2), name="mrg_t2w", run_without_submitting=True)

    norm_lap_tmpl = pe.Node(IntensityClip(), name="norm_lap_tmpl")
    norm_lap_t2w = pe.Node(IntensityClip(), name="norm_lap_t2w")

    # Set up initial spatial normalization
    ants_params = "testing" if sloppy else "precise"
    norm = pe.Node(
        Registration(
            from_file=pkgr_fn(
                "niworkflows.data", f"antsBrainExtraction_{ants_params}.json"
            )
        ),
        name="norm",
        n_procs=omp_nthreads,
        mem_gb=mem_gb,
    )
    norm.inputs.float = sloppy
    if tpl_regmask_path:
        norm.inputs.fixed_image_masks = tpl_regmask_path

    map_mask_t2w = pe.Node(
        ApplyTransforms(interpolation="Gaussian", float=True),
        name="map_mask_t2w",
        mem_gb=1,
    )

    # map template brainmask to t2w space
    map_mask_t2w.inputs.input_image = str(tpl_brainmask_path)

    thr_t2w_mask = pe.Node(Binarize(thresh_low=0.80), name="thr_t2w_mask")

    bspline_grid = pe.Node(
        niu.Function(function=_bspline_distance), name="bspline_grid"
    )

    # Refine INU correction
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
        # 1. Massage T2w
        (inputnode, val_t2w, [("in_t2w", "in_file")]),
        (val_t2w, clip_t2w, [("out_file", "in_file")]),
        (clip_t2w, init_t2w_n4, [("out_file", "input_image")]),
        (init_t2w_n4, clip_t2w_inu, [("output_image", "in_file")]),
        (clip_t2w_inu, lap_t2w, [("out_file", "op1")]),
        (lap_t2w, norm_lap_t2w, [("output_image", "in_file")]),
        (clip_t2w_inu, mrg_t2w, [("out_file", "in1")]),
        (norm_lap_t2w, mrg_t2w, [("out_file", "in2")]),
        # 2. Prepare template
        (clip_tmpl, lap_tmpl, [("out_file", "op1")]),
        (lap_tmpl, norm_lap_tmpl, [("output_image", "in_file")]),
        (clip_tmpl, mrg_tmpl, [("out_file", "in1")]),
        (norm_lap_tmpl, mrg_tmpl, [("out_file", "in2")]),
        # 3. Set normalization node inputs
        (mrg_tmpl, norm, [("out", "fixed_image")]),
        (mrg_t2w, norm, [("out", "moving_image")]),
        # 4. Map template brainmask into T2w space
        (val_t2w, map_mask_t2w, [("out_file", "reference_image")]),
        (norm, map_mask_t2w, [
            ("reverse_transforms", "transforms"),
            ("reverse_invert_flags", "invert_transform_flags")
        ]),
        (map_mask_t2w, thr_t2w_mask, [("output_image", "in_file")]),
        (thr_t2w_mask, apply_mask, [("out_mask", "in_mask")]),
        (final_n4, apply_mask, [("output_image", "in_file")]),
        # 5. Refine T2w INU correction with brain mask
        (clip_t2w, bspline_grid, [("out_file", "in_file")]),
        (clip_t2w, final_n4, [("out_file", "input_image")]),
        (bspline_grid, final_n4, [("out", "args")]),
        (map_mask_t2w, final_n4, [("output_image", "weight_image")]),
        # 9. Outputs
        (final_n4, outputnode, [("output_image", "t2w_preproc")]),
        (map_mask_t2w, outputnode, [("output_image", "out_probmap")]),
        (thr_t2w_mask, outputnode, [("out_mask", "out_mask")]),
        (apply_mask, outputnode, [("out_file", "t2w_brain")]),
    ])
    # fmt:on

    if ants_affine_init:
        from nipype.interfaces.ants.utils import AI

        ants_kwargs = dict(
            metric=("Mattes", 32, "Regular", 0.2),
            transform=("Affine", 0.1),
            search_factor=(20, 0.12),
            principal_axes=False,
            convergence=(10, 1e-6, 10),
            search_grid=(40, (0, 40, 40)),
            verbose=True,
        )

        if ants_affine_init == "random":
            ants_kwargs["metric"] = ("Mattes", 32, "Random", 0.2)
        if ants_affine_init == "search":
            ants_kwargs["search_grid"] = (20, (20, 40, 40))

        init_aff = pe.Node(
            AI(**ants_kwargs),
            name="init_aff",
            n_procs=omp_nthreads,
        )
        if tpl_regmask_path:
            init_aff.inputs.fixed_image_mask = _pop(tpl_regmask_path)

        # fmt:off
        workflow.connect([
            (clip_tmpl, init_aff, [("out_file", "fixed_image")]),
            (clip_t2w_inu, init_aff, [("out_file", "moving_image")]),
            (init_aff, norm, [("output_transform", "initial_moving_transform")]),
        ])
        # fmt:on

    return workflow


def _pop(in_files):
    if isinstance(in_files, (list, tuple)):
        return in_files[0]
    return in_files


def _bspline_distance(in_file, spacings=(20, 20, 20)):
    import numpy as np
    import nibabel as nb

    img = nb.load(in_file)
    extent = (np.array(img.shape[:3]) - 1) * img.header.get_zooms()[:3]
    retval = [f"{v}" for v in np.ceil(extent / np.array(spacings)).astype(int)]
    return f"-b {'x'.join(retval)}"
