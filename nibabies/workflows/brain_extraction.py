# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Nipype translation of ANTs' workflows."""
import numpy as np
# general purpose
from pkg_resources import resource_filename as pkgr_fn

# nipype
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.ants import N4BiasFieldCorrection
from nipype.interfaces.ants.utils import AI

# niworkflows
from niworkflows.interfaces.ants import ImageMath
from niworkflows.interfaces.images import RegridToZooms
from niworkflows.interfaces.nibabel import ApplyMask, Binarize
from niworkflows.interfaces.fixes import (
    FixHeaderRegistration as Registration,
    FixHeaderApplyTransforms as ApplyTransforms,
)
from niworkflows.interfaces.registration import (
    SimpleBeforeAfterRPT as SimpleBeforeAfter
)

from templateflow.api import get as get_template
from ..utils.filtering import (
    gaussian_filter as _gauss_filter,
    truncation as _trunc
)

LOWRES_ZOOMS = HIRES_ZOOMS = (4, 4, 4)


def init_infant_brain_extraction_wf(
    ants_affine_init=False,
    bspline_fitting_distance=8,
    debug=False,
    in_template="MNIInfant",
    template_specs=None,
    interim_checkpoints=True,
    mem_gb=3.0,
    mri_scheme="T2w",
    name="infant_brain_extraction_wf",
    omp_nthreads=None,
    output_dir=None,
    use_float=True,
):
    """
    Build an atlas-based brain extraction pipeline for infant T2w MRI data.

    Parameters
    ----------
    ants_affine_init : :obj:`bool`, optional
        Set-up a pre-initialization step with ``antsAI`` to account for mis-oriented images.

    """
    inputnode = pe.Node(
        niu.IdentityInterface(fields=["in_files", "in_mask"]), name="inputnode"
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["out_corrected", "out_brain", "out_mask"]),
        name="outputnode"
    )

    template_specs = template_specs or {}
    # Find a suitable target template in TemplateFlow
    tpl_target_path = get_template(
        in_template,
        suffix=mri_scheme,
        **template_specs
    )
    if not tpl_target_path:
        raise RuntimeError(
            f"An instance of template <tpl-{in_template}> with MR scheme '{mri_scheme}'"
            " could not be found.")

    tpl_brainmask_path = get_template(
        in_template, desc="brain", suffix="probseg", **template_specs
    )
    if not tpl_brainmask_path:
        tpl_brainmask_path = get_template(
            in_template, desc="brain", suffix="mask", **template_specs
        )

    tpl_regmask_path = get_template(
        in_template, desc="BrainCerebellumExtraction", suffix="mask", **template_specs
    )

    # Resample both target and template to a controlled, isotropic resolution
    res_tmpl = pe.Node(RegridToZooms(zooms=HIRES_ZOOMS), name="res_tmpl")
    res_target = pe.Node(RegridToZooms(zooms=HIRES_ZOOMS), name="res_target")
    gauss_tmpl = pe.Node(niu.Function(function=_gauss_filter), name="gauss_tmpl")
    gauss_tmpl.inputs.sigma = tuple(np.array(LOWRES_ZOOMS) * 10.0)

    # Spatial normalization step
    lap_tmpl = pe.Node(ImageMath(operation="Laplacian", op2="0.4 1"), name="lap_tmpl")
    lap_target = pe.Node(ImageMath(operation="Laplacian", op2="0.4 1"), name="lap_target")

    # Merge image nodes
    mrg_target = pe.Node(niu.Merge(2), name="mrg_target")
    mrg_tmpl = pe.Node(niu.Merge(2), name="mrg_tmpl")

    norm_lap_tmpl = pe.Node(niu.Function(function=_trunc), name="norm_lap_tmpl")
    norm_lap_tmpl.inputs.dtype = "float32"
    norm_lap_tmpl.inputs.out_max = 1.0
    norm_lap_tmpl.inputs.percentile = (0.01, 99.99)
    norm_lap_tmpl.inputs.clip_max = None

    norm_lap_target = pe.Node(niu.Function(function=_trunc), name="norm_lap_target")
    norm_lap_target.inputs.dtype = "float32"
    norm_lap_target.inputs.out_max = 1.0
    norm_lap_target.inputs.percentile = (0.01, 99.99)
    norm_lap_target.inputs.clip_max = None

    # Set up initial spatial normalization
    ants_params = "testing" if debug else "precise"
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

    # main workflow
    wf = pe.Workflow(name)
    # Create a buffer interface as a cache for the actual inputs to registration
    buffernode = pe.Node(niu.IdentityInterface(
        fields=["hires_target", "smooth_target"]), name="buffernode")

    # truncate target intensity for N4 correction
    clip_target = pe.Node(
        niu.Function(function=_trunc),
        name="clip_target",
    )
    clip_tmpl = pe.Node(
        niu.Function(function=_trunc),
        name="clip_tmpl",
    )
    clip_tmpl.inputs.in_file = _pop(tpl_target_path)

    # INU correction of the target image
    init_n4 = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            save_bias=False,
            copy_header=True,
            n_iterations=[50] * (4 - debug),
            convergence_threshold=1e-7,
            shrink_factor=4,
            bspline_fitting_distance=8,
        ),
        n_procs=omp_nthreads,
        name="init_n4",
    )
    clip_inu = pe.Node(
        niu.Function(function=_trunc),
        name="clip_inu",
    )
    gauss_target = pe.Node(niu.Function(function=_gauss_filter), name="gauss_target")
    gauss_target.inputs.sigma = tuple(np.array(LOWRES_ZOOMS) * 8.0)
    wf.connect([
        # truncation, resampling, and initial N4
        (inputnode, res_target, [(("in_files", _pop), "in_file")]),
        (res_target, clip_target, [("out_file", "in_file")]),
        (clip_tmpl, res_tmpl, [("out", "in_file")]),
        (clip_target, init_n4, [("out", "input_image")]),
        (init_n4, clip_inu, [("output_image", "in_file")]),
        (clip_inu, gauss_target, [("out", "in_file")]),
        (clip_inu, buffernode, [("out", "hires_target")]),
        (gauss_target, buffernode, [("out", "smooth_target")]),
        (res_tmpl, gauss_tmpl, [("out_file", "in_file")]),
    ])

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

    map_brainmask = pe.Node(
        ApplyTransforms(interpolation="Gaussian", float=True),
        name="map_brainmask",
        mem_gb=1
    )
    map_brainmask.inputs.input_image = str(tpl_brainmask_path)

    thr_brainmask = pe.Node(Binarize(thresh_low=0.80),
                            name="thr_brainmask")
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

    wf.connect([
        (inputnode, map_brainmask, [(("in_files", _pop), "reference_image")]),
        (inputnode, final_n4, [(("in_files", _pop), "input_image")]),
        (inputnode, bspline_grid, [(("in_files", _pop), "in_file")]),
        # (bspline_grid, final_n4, [("out", "bspline_fitting_distance")]),
        (bspline_grid, final_n4, [("out", "args")]),
        # merge laplacian and original images
        (buffernode, lap_target, [("smooth_target", "op1")]),
        (buffernode, mrg_target, [("hires_target", "in1")]),
        (lap_target, norm_lap_target, [("output_image", "in_file")]),
        (norm_lap_target, mrg_target, [("out", "in2")]),
        # Template massaging
        (res_tmpl, lap_tmpl, [("out_file", "op1")]),
        (res_tmpl, mrg_tmpl, [("out_file", "in1")]),
        (lap_tmpl, norm_lap_tmpl, [("output_image", "in_file")]),
        (norm_lap_tmpl, mrg_tmpl, [("out", "in2")]),
        # spatial normalization
        (mrg_target, norm, [("out", "moving_image")]),
        (mrg_tmpl, norm, [("out", "fixed_image")]),
        (norm, map_brainmask, [
            ("reverse_transforms", "transforms"),
            ("reverse_invert_flags", "invert_transform_flags")]),
        (map_brainmask, thr_brainmask, [("output_image", "in_file")]),
        # take a second pass of N4
        (map_brainmask, final_n4, [("output_image", "weight_image")]),
        (final_n4, final_mask, [("output_image", "in_file")]),
        (thr_brainmask, final_mask, [("out_mask", "in_mask")]),
        (final_n4, outputnode, [("output_image", "out_corrected")]),
        (thr_brainmask, outputnode, [("out_mask", "out_mask")]),
        (final_mask, outputnode, [("out_file", "out_brain")]),
    ])

    if tpl_regmask_path:
        wf.connect([
            (hires_mask, norm, [
                ("output_image", "fixed_image_masks")]),
        ])

    if interim_checkpoints:
        final_apply = pe.Node(
            ApplyTransforms(
                interpolation="BSpline",
                float=True),
            name="final_apply",
            mem_gb=1
        )
        final_report = pe.Node(SimpleBeforeAfter(
            before_label="tpl-WHS",
            after_label="target"),
            name="final_report"
        )
        wf.connect([
            (inputnode, final_apply, [(("in_files", _pop), "reference_image")]),
            (res_tmpl, final_apply, [("out_file", "input_image")]),
            (norm, final_apply, [
                ("reverse_transforms", "transforms"),
                ("reverse_invert_flags", "invert_transform_flags")]),
            (final_apply, final_report, [("output_image", "before")]),
            (outputnode, final_report, [("out_corrected", "after"),
                                        ("out_mask", "wm_seg")]),
        ])

    if output_dir:
        from nipype.interfaces.io import DataSink
        ds_final_inu = pe.Node(DataSink(base_directory=str(output_dir.parent)),
                               name="ds_final_inu")
        ds_final_msk = pe.Node(DataSink(base_directory=str(output_dir.parent)),
                               name="ds_final_msk")
        ds_report = pe.Node(DataSink(base_directory=str(output_dir.parent)),
                            name="ds_report")

        wf.connect([
            (outputnode, ds_final_inu, [
                ("out_corrected", f"{output_dir.name}.@inu_corrected")]),
            (outputnode, ds_final_msk, [
                ("out_mask", f"{output_dir.name}.@brainmask")]),
            (final_report, ds_report, [
                ("out_report", f"{output_dir.name}.@report")]),
        ])

    if not ants_affine_init:
        return wf

    # Initialize transforms with antsAI
    lowres_tmpl = pe.Node(RegridToZooms(zooms=LOWRES_ZOOMS), name="lowres_tmpl")
    lowres_target = pe.Node(RegridToZooms(zooms=LOWRES_ZOOMS), name="lowres_target")

    init_aff = pe.Node(
        AI(
            metric=("Mattes", 32, "Regular", 1.0),
            transform=("Affine", 0.1),
            search_factor=(10, 0.08),
            principal_axes=False,
            convergence=(40, 1e-6, 10),
            search_grid=(25, (0, 0, 0)) if debug else (2, (1, 1, 1)),
            verbose=True,
        ),
        name="init_aff",
        n_procs=omp_nthreads,
    )
    wf.connect([
        (gauss_tmpl, lowres_tmpl, [("out", "in_file")]),
        (lowres_tmpl, init_aff, [("out_file", "fixed_image")]),
        (gauss_target, lowres_target, [("out", "in_file")]),
        (lowres_target, init_aff, [("out_file", "moving_image")]),
        (init_aff, norm, [("output_transform", "initial_moving_transform")]),
    ])

    if tpl_regmask_path:
        lowres_mask = pe.Node(
            ApplyTransforms(
                input_image=_pop(tpl_regmask_path),
                transforms="identity",
                interpolation="MultiLabel",
                float=True),
            name="lowres_mask",
            mem_gb=1
        )
        wf.connect([
            (lowres_tmpl, lowres_mask, [("out_file", "reference_image")]),
            (lowres_mask, init_aff, [("output_image", "fixed_image_mask")]),
        ])

    if interim_checkpoints:
        init_apply = pe.Node(
            ApplyTransforms(
                interpolation="BSpline",
                float=True),
            name="init_apply",
            mem_gb=1
        )
        init_report = pe.Node(SimpleBeforeAfter(
            before_label="tpl-WHS",
            after_label="target"),
            name="init_report"
        )
        wf.connect([
            (lowres_target, init_apply, [("out_file", "input_image")]),
            (res_tmpl, init_apply, [("out_file", "reference_image")]),
            (init_aff, init_apply, [("output_transform", "transforms")]),
            (init_apply, init_report, [("output_image", "after")]),
            (res_tmpl, init_report, [("out_file", "before")]),
        ])

    return wf


def _pop(in_files):
    if isinstance(in_files, (list, tuple)):
        return in_files[0]
    return in_files


def _bspline_distance(in_file, spacings=(8, 2, 8)):
    import numpy as np
    import nibabel as nb

    img = nb.load(in_file)
    extent = (np.array(img.shape[:3]) - 1) * img.header.get_zooms()[:3]
    retval = [f"{v}" for v in np.ceil(extent / np.array(spacings)).astype(int)]
    return f"-b {'x'.join(retval)}"
