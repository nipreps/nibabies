# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Baby brain extraction from T2w images."""
from pkg_resources import resource_filename as pkgr_fn

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
    debug=False,
    name="infant_brain_extraction_wf",
    omp_nthreads=None,
):
    """
    Build an atlas-based brain extraction pipeline for infant T2w MRI data.

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
    age_months : :obj:`int`
        Age of this participant, in months.
    ants_affine_init : :obj:`bool`, optional
        Set-up a pre-initialization step with ``antsAI`` to account for mis-oriented images.
    bspline_fitting_distance : :obj:`float`
        Distance in mm between B-Spline control points for N4 INU estimation.
    sloppy : :obj:`bool`
        Run in *sloppy* mode.
    skull_strip_template : :obj:`str`
        A TemplateFlow ID indicating which template will be used as target for atlas-based
        segmentation.
    template_specs : :obj:`dict`
        Additional template specifications (e.g., resolution or cohort) to correctly select
        the adequate template instance.
    mem_gb : :obj:`float`
        Base memory fingerprint unit.
    name : :obj:`str`
        This particular workflow's unique name (Nipype requirement).
    omp_nthreads : :obj:`int`
        The number of threads for individual processes in this workflow.
    debug : :obj:`bool`
        Produce intermediate registration files

    Inputs
    ------
    in_t2w : :obj:`str`
        The unprocessed input T2w image.

    Outputs
    -------
    t2w_preproc : :obj:`str`
        The preprocessed T2w image (INU and clipping).
    t2w_brain : :obj:`str`
        The preprocessed, brain-extracted T2w image.
    out_mask : :obj:`str`
        The brainmask projected from the template into the T2w, after
        binarization.
    out_probmap : :obj:`str`
        The same as above, before binarization.

    """
    from nipype.interfaces.ants import N4BiasFieldCorrection, ImageMath

    # niworkflows
    from niworkflows.interfaces.nibabel import ApplyMask, Binarize, IntensityClip
    from niworkflows.interfaces.fixes import (
        FixHeaderRegistration as Registration,
        FixHeaderApplyTransforms as ApplyTransforms,
    )
    from templateflow.api import get as get_template

    from ...interfaces.nibabel import BinaryDilation
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

    # Ensure template comes with a range of intensities ANTs will like
    clip_tmpl = pe.Node(IntensityClip(p_max=99), name="clip_tmpl")
    clip_tmpl.inputs.in_file = _pop(tpl_target_path)

    # Generate laplacian registration targets
    lap_tmpl = pe.Node(ImageMath(operation="Laplacian", op2="0.4 1"), name="lap_tmpl")
    lap_t2w = pe.Node(ImageMath(operation="Laplacian", op2="0.4 1"), name="lap_t2w")
    norm_lap_tmpl = pe.Node(niu.Function(function=_norm_lap), name="norm_lap_tmpl")
    norm_lap_t2w = pe.Node(niu.Function(function=_norm_lap), name="norm_lap_t2w")

    # Merge image nodes
    mrg_tmpl = pe.Node(niu.Merge(2), name="mrg_tmpl", run_without_submitting=True)
    mrg_t2w = pe.Node(niu.Merge(2), name="mrg_t2w", run_without_submitting=True)
    bin_regmask = pe.Node(Binarize(thresh_low=0.20), name="bin_regmask")
    bin_regmask.inputs.in_file = str(tpl_brainmask_path)
    refine_mask = pe.Node(BinaryDilation(radius=3, iterations=2), name="refine_mask")

    fixed_masks = pe.Node(niu.Merge(4), name="fixed_masks", run_without_submitting=True)
    fixed_masks.inputs.in1 = "NULL"
    fixed_masks.inputs.in2 = "NULL"
    fixed_masks.inputs.in3 = "NULL" if not tpl_regmask_path else _pop(tpl_regmask_path)

    # Set up initial spatial normalization
    ants_params = "testing" if sloppy else "precise"
    norm = pe.Node(
        Registration(
            from_file=pkgr_fn(
                "nibabies.data", f"antsBrainExtraction_{ants_params}.json"
            )
        ),
        name="norm",
        n_procs=omp_nthreads,
        mem_gb=mem_gb,
    )
    norm.inputs.float = sloppy
    if debug:
        norm.inputs.args = "--write-interval-volumes 5"

    map_mask_t2w = pe.Node(
        ApplyTransforms(interpolation="Gaussian", float=True),
        name="map_mask_t2w",
        mem_gb=1,
    )

    # map template brainmask to t2w space
    map_mask_t2w.inputs.input_image = str(tpl_brainmask_path)

    thr_t2w_mask = pe.Node(Binarize(thresh_low=0.80), name="thr_t2w_mask")

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
    final_clip = pe.Node(IntensityClip(p_min=5.0, p_max=99.5), name="final_clip")
    apply_mask = pe.Node(ApplyMask(), name="apply_mask")

    # fmt:off
    workflow.connect([
        (inputnode, final_n4, [("in_t2w", "input_image")]),
        # 1. Massage T2w
        (inputnode, mrg_t2w, [("in_t2w", "in1")]),
        (inputnode, lap_t2w, [("in_t2w", "op1")]),
        (inputnode, map_mask_t2w, [("in_t2w", "reference_image")]),
        (bin_regmask, refine_mask, [("out_file", "in_file")]),
        (refine_mask, fixed_masks, [("out_file", "in4")]),
        (lap_t2w, norm_lap_t2w, [("output_image", "in_file")]),
        (norm_lap_t2w, mrg_t2w, [("out", "in2")]),
        # 2. Prepare template
        (clip_tmpl, lap_tmpl, [("out_file", "op1")]),
        (lap_tmpl, norm_lap_tmpl, [("output_image", "in_file")]),
        (clip_tmpl, mrg_tmpl, [("out_file", "in1")]),
        (norm_lap_tmpl, mrg_tmpl, [("out", "in2")]),
        # 3. Set normalization node inputs
        (mrg_tmpl, norm, [("out", "fixed_image")]),
        (mrg_t2w, norm, [("out", "moving_image")]),
        (fixed_masks, norm, [("out", "fixed_image_masks")]),
        # 4. Map template brainmask into T2w space
        (norm, map_mask_t2w, [
            ("reverse_transforms", "transforms"),
            ("reverse_invert_flags", "invert_transform_flags")
        ]),
        (map_mask_t2w, thr_t2w_mask, [("output_image", "in_file")]),
        (thr_t2w_mask, apply_mask, [("out_mask", "in_mask")]),
        (final_n4, apply_mask, [("output_image", "in_file")]),
        # 5. Refine T2w INU correction with brain mask
        (map_mask_t2w, final_n4, [("output_image", "weight_image")]),
        (final_n4, final_clip, [("output_image", "in_file")]),
        # 9. Outputs
        (final_clip, outputnode, [("out_file", "t2w_preproc")]),
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
            (inputnode, init_aff, [("in_t2w", "moving_image")]),
            (init_aff, norm, [("output_transform", "initial_moving_transform")]),
        ])
        # fmt:on

    return workflow


def _pop(in_files):
    if isinstance(in_files, (list, tuple)):
        return in_files[0]
    return in_files


def _norm_lap(in_file):
    from pathlib import Path
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    img = nb.load(in_file)
    data = img.get_fdata()
    data -= np.median(data)
    l_max = np.percentile(data[data > 0], 99.8)
    l_min = np.percentile(data[data < 0], 0.2)
    data[data < 0] *= -1.0 / l_min
    data[data > 0] *= 1.0 / l_max
    data = np.clip(data, a_min=-1.0, a_max=1.0)

    out_file = fname_presuffix(
        Path(in_file).name, suffix="_norm", newpath=str(Path.cwd().absolute())
    )
    hdr = img.header.copy()
    hdr.set_data_dtype("float32")
    img.__class__(data.astype("float32"), img.affine, hdr).to_filename(out_file)
    return out_file
