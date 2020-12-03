# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Nipype translation of ANTs' workflows."""
# import numpy as np
# general purpose
from pkg_resources import resource_filename as pkgr_fn

# nipype
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.ants import N4BiasFieldCorrection, ImageMath
# from nipype.interfaces.ants.utils import AI

# niworkflows
from niworkflows.anat.ants import init_atropos_wf, ATROPOS_MODELS
from niworkflows.interfaces.images import RegridToZooms, ValidateImage
from niworkflows.interfaces.nibabel import ApplyMask, Binarize
from niworkflows.interfaces.fixes import (
    FixHeaderRegistration as Registration,
    FixHeaderApplyTransforms as ApplyTransforms,
)
from niworkflows.interfaces.registration import (
    SimpleBeforeAfterRPT as SimpleBeforeAfter
)

from templateflow.api import get as get_template
from ...utils.filtering import (
    gaussian_filter as _gauss_filter,
    truncation as _trunc
)
from ...utils.misc import cohort_by_months

HIRES_ZOOMS = (1, 1, 1)
LOWRES_ZOOMS = (2, 2, 2)


def init_infant_brain_extraction_wf(
    age_months=None,
    ants_affine_init=False,
    bspline_fitting_distance=200,
    sloppy=False,
    skull_strip_template="UNCInfant",
    template_specs=None,
    interim_checkpoints=True,
    mem_gb=3.0,
    mri_scheme="T1w",
    name="infant_brain_extraction_wf",
    atropos_model=None,
    omp_nthreads=None,
    output_dir=None,
    use_float=True,
    use_t2w=False,
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
    # handle template specifics
    template_specs = template_specs or {}
    if skull_strip_template == 'MNIInfant':
        template_specs['resolution'] = 2 if sloppy else 1

    if not template_specs.get('cohort'):
        if age_months is None:
            raise KeyError(f"Age or cohort for {skull_strip_template} must be provided!")
        template_specs['cohort'] = cohort_by_months(skull_strip_template, age_months)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["t1w", "t2w", "in_mask"]), name="inputnode"
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["t1w_corrected", "t1w_corrected_brain", "t1w_mask"]),
        name="outputnode"
    )

    if not use_t2w:
        raise RuntimeError("A T2w image is currently required.")

    tpl_target_path = get_template(
        skull_strip_template,
        suffix='T1w',  # no T2w template
        desc=None,
        **template_specs,
    )
    if not tpl_target_path:
        raise RuntimeError(
            f"An instance of template <tpl-{skull_strip_template}> with MR scheme "
            f"'{'T1w' or mri_scheme}' could not be found."
        )

    tpl_brainmask_path = get_template(
        skull_strip_template, desc="brain", suffix="probseg", **template_specs
    )
    if not tpl_brainmask_path:
        tpl_brainmask_path = get_template(
            skull_strip_template, desc="brain", suffix="mask", **template_specs
        )

    tpl_regmask_path = get_template(
        skull_strip_template, desc="BrainCerebellumExtraction", suffix="mask", **template_specs
    )

    # validate images
    val_tmpl = pe.Node(ValidateImage(), name='val_tmpl')
    val_tmpl.inputs.in_file = _pop(tpl_target_path)
    # val_t1w = val_tmpl.clone("val_t1w")
    # val_t2w = val_tmpl.clone("val_t2w")

    # Resample both target and template to a controlled, isotropic resolution
    res_tmpl = pe.Node(RegridToZooms(zooms=HIRES_ZOOMS), name="res_tmpl")
    res_t1w = res_tmpl.clone("res_t1w")
    res_t2w = res_tmpl.clone("res_t2w")

    gauss_tmpl = pe.Node(niu.Function(function=_gauss_filter), name="gauss_tmpl")
    gauss_t1w = gauss_tmpl.clone("gauss_t1w")
    gauss_t2w = gauss_tmpl.clone("gauss_t2w")

    # Spatial normalization step
    lap_tmpl = pe.Node(ImageMath(operation="Laplacian", op2="0.4 1"), name="lap_tmpl")
    lap_t1w = lap_tmpl.clone("lap_t1w")
    lap_t2w = lap_tmpl.clone("lap_t2w")

    # Merge image nodes
    mrg_tmpl = pe.Node(niu.Merge(2), name="mrg_tmpl")
    mrg_t2w = mrg_tmpl.clone("mrg_t2w")
    mrg_t1w = mrg_tmpl.clone("mrg_t1w")

    norm_lap_tmpl = pe.Node(niu.Function(function=_trunc), name="norm_lap_tmpl")
    norm_lap_tmpl.inputs.dtype = "float32"
    norm_lap_tmpl.inputs.out_max = 1.0
    norm_lap_tmpl.inputs.percentile = (0.01, 99.99)
    norm_lap_tmpl.inputs.clip_max = None

    norm_lap_t1w = norm_lap_tmpl.clone('norm_lap_t1w')
    norm_lap_t2w = norm_lap_t1w.clone('norm_lap_t2w')

    # Set up initial spatial normalization
    ants_params = "testing" if sloppy else "precise"
    norm = pe.Node(
        Registration(from_file=pkgr_fn(
            "niworkflows.data",
            f"antsBrainExtraction_{ants_params}.json")
        ),
        name="norm",
        n_procs=omp_nthreads,
        mem_gb=mem_gb,
    )
    norm.inputs.float = use_float

    # Set up T2w -> T1w within-subject registration
    norm_subj = pe.Node(
        Registration(from_file=pkgr_fn("nibabies.data", "within_subject.json")),
        name="norm_subj",
        n_procs=omp_nthreads,
        mem_gb=mem_gb,
    )
    norm_subj.inputs.float = use_float

    # main workflow
    wf = pe.Workflow(name)
    # Create a buffer interface as a cache for the actual inputs to registration
    buffernode = pe.Node(niu.IdentityInterface(
        fields=["hires_target", "smooth_target"]), name="buffernode")

    # truncate target intensity for N4 correction
    clip_tmpl = pe.Node(
        niu.Function(function=_trunc),
        name="clip_tmpl",
    )
    clip_t2w = clip_tmpl.clone('clip_t2w')
    clip_t1w = clip_tmpl.clone('clip_t1w')

    # INU correction of the t1w
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
    clip_inu = pe.Node(
        niu.Function(function=_trunc),
        name="clip_inu",
    )

    # Graft a template registration-mask if present
    if tpl_regmask_path:
        hires_mask = pe.Node(
            ApplyTransforms(
                input_image=_pop(tpl_regmask_path),
                transforms="identity",
                interpolation="NearestNeighbor",
                float=True),
            name="hires_mask",
            mem_gb=1
        )
        wf.connect([
            (res_tmpl, hires_mask, [("out_file", "reference_image")]),
        ])

    map_mask_t2w = pe.Node(
        ApplyTransforms(interpolation="Gaussian", float=True),
        name="map_mask_t2w",
        mem_gb=1
    )
    map_mask_t1w = map_mask_t2w.clone("map_mask_t1w")

    # map template brainmask to t2w space
    map_mask_t2w.inputs.input_image = str(tpl_brainmask_path)

    thr_t2w_mask = pe.Node(Binarize(thresh_low=0.80), name="thr_t2w_mask")
    thr_t1w_mask = thr_t2w_mask.clone('thr_t1w_mask')

    bspline_grid = pe.Node(niu.Function(function=_bspline_distance),
                           name="bspline_grid")

    # Refine INU correction
    final_n4 = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
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
    final_mask = pe.Node(ApplyMask(), name="final_mask")

    if atropos_model is None:
        atropos_model = tuple(ATROPOS_MODELS[mri_scheme].values())

    atropos_wf = init_atropos_wf(
        use_random_seed=False,
        omp_nthreads=omp_nthreads,
        mem_gb=mem_gb,
        in_segmentation_model=atropos_model,
    )
    # if tpl_regmask_path:
    #     atropos_wf.get_node('inputnode').inputs.in_mask_dilated = tpl_regmask_path

    sel_wm = pe.Node(niu.Select(index=atropos_model[-1] - 1), name='sel_wm',
                     run_without_submitting=True)

    wf.connect([
        # 1. massage template
        (val_tmpl, clip_tmpl, [("out_file", "in_file")]),
        (clip_tmpl, res_tmpl, [("out", "in_file")]),
        (res_tmpl, gauss_tmpl, [("out_file", "in_file")]),
        (res_tmpl, lap_tmpl, [("out_file", "op1")]),
        (res_tmpl, mrg_tmpl, [("out_file", "in1")]),
        (lap_tmpl, norm_lap_tmpl, [("output_image", "in_file")]),
        (norm_lap_tmpl, mrg_tmpl, [("out", "in2")]),
        # 2. massage T2w
        (inputnode, val_t2w, [('t2w', 'in_file')]),
        (val_t2w, res_t2w, [('out_file', 'in_file')]),
        (res_t2w, clip_t2w, [('out_file', 'in_file')]),
        # no need for init_n4 + clip_inu of T2w
        (clip_t2w, gauss_t2w, [('out', 'in_file')]),
        (clip_t2w, buffernode, [('out', 'hires_target')]),
        (gauss_t2w, buffernode, [('out', 'smooth_target')]),
        (buffernode, lap_t2w, [("smooth_target", "op1")]),
        (lap_t2w, norm_lap_t2w, [("output_image", "in_file")]),
        (norm_lap_t2w, mrg_t2w, [("out", "in2")]),
        # 3. normalize T2w to target template (UNC)
        (mrg_t2w, norm, [("out", "moving_image")]),
        (mrg_tmpl, norm, [("out", "fixed_image")]),
        # 4. map template brainmask to T2w space
        (inputnode, map_mask_t2w, [('t2w', 'reference_image')]),
        (norm, map_mask_t2w, [
            ("reverse_transforms", "transforms"),
            ("reverse_invert_flags", "invert_transform_flags")
        ]),
        (map_mask_t2w, thr_t2w_mask, [("output_image", "in_file")]),
        # 5. massage T1w
        (inputnode, res_t1w, [('t1w', 'in_file')]),
        (res_t1w, clip_t1w, [('out_file', 'in_file')]),
        (clip_t1w, init_n4, [("out", "input_image")]),
        (init_n4, clip_inu, [("output_image", "in_file")]),
        (clip_inu, gauss_t1w, [("out", "in_file")]),
        (clip_inu, mrg_t1w, [("out", "in1")]),
        (gauss_t1w, lap_t1w, [("out", "op1")]),
        (lap_t1w, norm_lap_t1w, [("output_image", "in_file")]),
        (norm_lap_t1w, mrg_t1w, [("out", "in2")]),
        # 6. normalize within subject T2w to T1w
        (mrg_t2w, norm_subj, [("out", "moving_image")]),
        (mrg_t1w, norm_subj, [("out", "fixed_image")]),
        # 7. map mask to T1w space
        (thr_t2w_mask, map_mask_t1w, [("out_mask", "input_image")]),
        (inputnode, map_mask_t1w, [("t1w", "reference_image")]),
        (norm_subj, map_mask_t1w, [
            ("reverse_transforms", "transforms"),
            ("reverse_invert_flags", "invert_transform_flags"),
        ]),
        (map_mask_t1w, thr_t1w_mask, [("output_image", "in_file")]),
        # 8. T1w INU
        (inputnode, final_n4, [("t1w", "input_image")]),
        (inputnode, bspline_grid, [("t1w", "in_file")]),
        (bspline_grid, final_n4, [("out", "args")]),
        (map_mask_t1w, final_n4, [("output_image", "weight_image")]),
        (final_n4, final_mask, [("output_image", "in_file")]),
        (thr_t1w_mask, final_mask, [("out_mask", "in_mask")]),
        # 9. Outputs
        (final_n4, outputnode, [("output_image", "t1w_corrected")]),
        (thr_t1w_mask, outputnode, [("out_mask", "t1w_mask")]),
        (final_mask, outputnode, [("out_file", "t1w_corrected_brain")]),
    ])

    if tpl_regmask_path:
        wf.connect([
            (hires_mask, norm, [
                ("output_image", "fixed_image_masks")]),
        ])

    # if interim_checkpoints:
    #     final_apply = pe.Node(
    #         ApplyTransforms(
    #             interpolation="BSpline",
    #             float=True),
    #         name="final_apply",
    #         mem_gb=1
    #     )
    #     final_report = pe.Node(SimpleBeforeAfter(
    #         before_label=f"tpl-{skull_strip_template}",
    #         after_label="target",
    #         out_report="final_report.svg"),
    #         name="final_report"
    #     )
    #     wf.connect([
    #         (inputnode, final_apply, [("in_file", "reference_image")]),
    #         (res_tmpl, final_apply, [("out_file", "input_image")]),
    #         (norm, final_apply, [
    #             ("reverse_transforms", "transforms"),
    #             ("reverse_invert_flags", "invert_transform_flags")]),
    #         (final_apply, final_report, [("output_image", "before")]),
    #         (outputnode, final_report, [("out_corrected", "after"),
    #                                     ("out_mask", "wm_seg")]),
    #     ])

    # if output_dir:
    #     from nipype.interfaces.io import DataSink
    #     ds_final_inu = pe.Node(DataSink(base_directory=str(output_dir.parent)),
    #                            name="ds_final_inu")
    #     ds_final_msk = pe.Node(DataSink(base_directory=str(output_dir.parent)),
    #                            name="ds_final_msk")
    #     ds_report = pe.Node(DataSink(base_directory=str(output_dir.parent)),
    #                         name="ds_report")

    #     wf.connect([
    #         (outputnode, ds_final_inu, [
    #             ("out_corrected", f"{output_dir.name}.@inu_corrected")]),
    #         (outputnode, ds_final_msk, [
    #             ("out_mask", f"{output_dir.name}.@brainmask")]),
    #         (final_report, ds_report, [
    #             ("out_report", f"{output_dir.name}.@report")]),
    #     ])

    # if not ants_affine_init:
    #     return wf

    # # Initialize transforms with antsAI
    # lowres_tmpl = pe.Node(RegridToZooms(zooms=LOWRES_ZOOMS), name="lowres_tmpl")
    # lowres_target = pe.Node(RegridToZooms(zooms=LOWRES_ZOOMS), name="lowres_target")

    # init_aff = pe.Node(
    #     AI(
    #         metric=("Mattes", 32, "Regular", 0.25),
    #         transform=("Affine", 0.1),
    #         search_factor=(15, 0.1),
    #         principal_axes=False,
    #         convergence=(10, 1e-6, 10),
    #         search_grid=(40, (0, 40, 40)),
    #         verbose=True,
    #     ),
    #     name="init_aff",
    #     n_procs=omp_nthreads,
    # )
    # wf.connect([
    #     (gauss_tmpl, lowres_tmpl, [("out", "in_file")]),
    #     (lowres_tmpl, init_aff, [("out_file", "fixed_image")]),
    #     (gauss_target, lowres_target, [("out", "in_file")]),
    #     (lowres_target, init_aff, [("out_file", "moving_image")]),
    #     (init_aff, norm, [("output_transform", "initial_moving_transform")]),
    # ])

    # if tpl_regmask_path:
    #     lowres_mask = pe.Node(
    #         ApplyTransforms(
    #             input_image=_pop(tpl_regmask_path),
    #             transforms="identity",
    #             interpolation="MultiLabel",
    #             float=True),
    #         name="lowres_mask",
    #         mem_gb=1
    #     )
    #     wf.connect([
    #         (lowres_tmpl, lowres_mask, [("out_file", "reference_image")]),
    #         (lowres_mask, init_aff, [("output_image", "fixed_image_mask")]),
    #     ])

    # if interim_checkpoints:
    #     init_apply = pe.Node(
    #         ApplyTransforms(
    #             interpolation="BSpline",
    #             float=True),
    #         name="init_apply",
    #         mem_gb=1
    #     )
    #     init_report = pe.Node(SimpleBeforeAfter(
    #         before_label=f"tpl-{skull_strip_template}",
    #         after_label="target",
    #         out_report="init_report.svg"),
    #         name="init_report"
    #     )
    #     wf.connect([
    #         (lowres_target, init_apply, [("out_file", "input_image")]),
    #         (res_tmpl, init_apply, [("out_file", "reference_image")]),
    #         (init_aff, init_apply, [("output_transform", "transforms")]),
    #         (init_apply, init_report, [("output_image", "after")]),
    #         (res_tmpl, init_report, [("out_file", "before")]),
    #     ])

    #     if output_dir:
    #         ds_init_report = pe.Node(DataSink(base_directory=str(output_dir.parent)),
    #                 name="ds_init_report")
    #         wf.connect(
    #             init_report, "out_report",
    #             ds_init_report, f"{output_dir.name}.@init_report"
    #         )
    return wf


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
