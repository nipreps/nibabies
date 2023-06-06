"""Base anatomical preprocessing."""
import warnings
from pathlib import Path
from typing import Literal, Optional, Union

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow
from niworkflows.utils.spaces import Reference, SpatialReferences

from ... import config


def init_infant_anat_wf(
    *,
    age_months: Optional[int],
    ants_affine_init: bool,
    t1w: list,
    t2w: list,
    anat_modality: str,
    bids_root: Optional[Union[str, Path]],
    existing_derivatives: dict,
    freesurfer: bool,
    hires: Optional[bool],
    longitudinal: bool,
    omp_nthreads: int,
    output_dir: Union[str, Path],
    segmentation_atlases: Optional[Union[str, Path]],
    skull_strip_mode: str,
    skull_strip_template: Reference,
    sloppy: bool,
    spaces: Optional[SpatialReferences],
    cifti_output: Optional[Literal['91k', '170k']],
    name: str = "infant_anat_wf",
) -> LiterateWorkflow:
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
    t1w2fsnative_xfm
        LTA-style affine matrix translating from T1w to
        FreeSurfer-conformed subject space
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed
        subject space to T1w
    surfaces
        GIFTI surfaces (gray/white boundary, midthickness, pial, inflated)
    """
    from nipype.interfaces.ants.base import Info as ANTsInfo
    from smriprep.workflows.norm import init_anat_norm_wf

    from ...utils.misc import fix_multi_source_name
    from .brain_extraction import init_infant_brain_extraction_wf
    from .outputs import (
        init_anat_derivatives_wf,
        init_anat_reports_wf,
        init_coreg_report_wf,
    )
    from .preproc import init_anat_preproc_wf
    from .registration import init_coregistration_wf
    from .segmentation import init_anat_segmentations_wf
    from .surfaces import init_anat_ribbon_wf
    from .template import init_anat_template_wf

    # for now, T1w only
    num_t1w = len(t1w) if t1w else 0
    num_t2w = len(t2w) if t2w else 0

    precomp_mask = existing_derivatives.get("anat_mask")
    precomp_aseg = existing_derivatives.get("anat_aseg")

    # verify derivatives are relatively similar to T1w
    if precomp_mask or precomp_aseg:
        if num_t1w > 1:
            precomp_mask = None
            precomp_aseg = None
            warnings.warn(
                "Multiple T1w files were found; precomputed derivatives will not be used."
            )

        else:
            from ...utils.validation import validate_t1w_derivatives

            validated_derivatives = (
                validate_t1w_derivatives(  # compare derivatives to the first T1w
                    t1w[0], anat_mask=precomp_mask, anat_aseg=precomp_aseg
                )
            )
            precomp_mask = validated_derivatives.get("anat_mask")
            precomp_aseg = validated_derivatives.get("anat_aseg")

    wf = LiterateWorkflow(name=name)
    desc = f"""\n
Anatomical data preprocessing

: A total of {num_t1w} T1w and {num_t2w} T2w images were found within the input
BIDS dataset."""

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["t1w", "t2w", "subject_id", "subjects_dir"]),  # FLAIR / ROI?
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
                "t1w2fsnative_xfm",
                "fsnative2t1w_xfm",
                "surfaces",
                "morphometrics",
                "anat_aseg",
                "anat_mcrib",
                "anat_aparc",
                "anat_ribbon",
                "template",
            ]
        ),
        name="outputnode",
    )

    desc += (
        """\
All of the T1-weighted images were denoised and corrected for intensity non-uniformity (INU)"""
        if num_t1w > 1
        else """\
The T1-weighted (T1w) image was denoised and corrected for intensity non-uniformity (INU)"""
    )

    desc += """\
with `N4BiasFieldCorrection` [@n4], distributed with ANTs {ants_ver} \
[@ants, RRID:SCR_004757]"""
    desc += ".\n" if num_t1w > 1 else ", and used as T1w-reference throughout the workflow.\n"

    desc += (
        "A previously computed mask was used to skull-strip the anatomical image."
        if precomp_mask
        else """\
The T1w-reference was then skull-stripped with a modified implementation of
the `antsBrainExtraction.sh` workflow (from ANTs), using {skullstrip_tpl}
as target template.
"""
    )

    wf.__desc__ = desc.format(
        ants_ver=ANTsInfo.version() or "(version unknown)",
        skullstrip_tpl=skull_strip_template.fullname,
    )
    wf.__postdesc__ = ""

    # Define output workflows
    anat_reports_wf = init_anat_reports_wf(
        freesurfer=freesurfer, output_dir=output_dir, sloppy=sloppy
    )

    anat_derivatives_wf = init_anat_derivatives_wf(
        bids_root=bids_root,
        freesurfer=freesurfer,
        num_t1w=num_t1w,
        output_dir=output_dir,
        spaces=spaces,
        cifti_output=cifti_output,
    )

    # Multiple anatomical files -> generate average reference
    t1w_template_wf = init_anat_template_wf(
        contrast="T1w",
        num_files=num_t1w,
        longitudinal=longitudinal,
        omp_nthreads=omp_nthreads,
        sloppy=sloppy,
        name="t1w_template_wf",
    )

    t2w_template_wf = init_anat_template_wf(
        contrast="T2w",
        num_files=num_t2w,
        longitudinal=longitudinal,
        omp_nthreads=omp_nthreads,
        sloppy=sloppy,
        name="t2w_template_wf",
    )

    # Clean up each anatomical template
    # Denoise, INU, + Clipping
    t1w_preproc_wf = init_anat_preproc_wf(
        precomputed_mask=bool(precomp_mask),
        precomputed_aseg=bool(precomp_aseg),
        name="t1w_preproc_wf",
    )
    t2w_preproc_wf = init_anat_preproc_wf(name="t2w_preproc_wf")

    if skull_strip_mode != "force":
        raise NotImplementedError("Skull stripping is currently required.")

    coregistration_wf = init_coregistration_wf(
        omp_nthreads=omp_nthreads,
        sloppy=sloppy,
        debug="registration" in config.execution.debug,
        precomputed_mask=bool(precomp_mask),
    )
    coreg_report_wf = init_coreg_report_wf(
        output_dir=output_dir,
    )

    # Segmentation - initial implementation should be simple: JLF
    anat_seg_wf = init_anat_segmentations_wf(
        anat_modality=anat_modality.capitalize(),  # TODO: Revisit this option
        template_dir=segmentation_atlases,
        sloppy=sloppy,
        omp_nthreads=omp_nthreads,
        precomp_aseg=precomp_aseg,
    )

    # Spatial normalization (requires segmentation)
    anat_norm_wf = init_anat_norm_wf(
        sloppy=sloppy,
        omp_nthreads=omp_nthreads,
        templates=spaces.get_spaces(nonstandard=False, dim=(3,)),
    )

    # fmt:off
    wf.connect([
        (inputnode, t1w_template_wf, [("t1w", "inputnode.anat_files")]),
        (inputnode, t2w_template_wf, [("t2w", "inputnode.anat_files")]),
        (inputnode, anat_reports_wf, [("t1w", "inputnode.source_file")]),
        (inputnode, coreg_report_wf, [("t1w", "inputnode.source_file")]),
        (inputnode, anat_norm_wf, [(("t1w", fix_multi_source_name), "inputnode.orig_t1w")]),

        (t1w_template_wf, outputnode, [
            ("outputnode.anat_realign_xfm", "anat_ref_xfms")]),
        (t1w_template_wf, t1w_preproc_wf, [("outputnode.anat_ref", "inputnode.in_anat")]),
        (t1w_template_wf, anat_derivatives_wf, [
            ("outputnode.anat_valid_list", "inputnode.source_files"),
            ("outputnode.anat_realign_xfm", "inputnode.t1w_ref_xfms")]),
        (t1w_template_wf, anat_reports_wf, [
            ("outputnode.out_report", "inputnode.t1w_conform_report")]),

        (t2w_template_wf, t2w_preproc_wf, [("outputnode.anat_ref", "inputnode.in_anat")]),
        (t2w_template_wf, anat_derivatives_wf, [
            ("outputnode.anat_valid_list", "inputnode.t2w_source_files")]),

        (t1w_preproc_wf, coregistration_wf, [("outputnode.anat_preproc", "inputnode.in_t1w")]),
        (t1w_preproc_wf, coreg_report_wf, [("outputnode.anat_preproc", "inputnode.t1w_preproc")]),

        (coregistration_wf, coreg_report_wf, [
            ("outputnode.t1w_mask", "inputnode.in_mask"),
            ("outputnode.t2w_preproc", "inputnode.t2w_preproc")]),
        (coregistration_wf, anat_norm_wf, [
            ("outputnode.t1w_preproc", "inputnode.moving_image"),
            ("outputnode.t1w_mask", "inputnode.moving_mask")]),
        (coregistration_wf, anat_seg_wf, [("outputnode.t1w_brain", "inputnode.anat_brain")]),
        (coregistration_wf, anat_derivatives_wf, [
            ("outputnode.t1w_mask", "inputnode.t1w_mask"),
            ("outputnode.t1w_preproc", "inputnode.t1w_preproc"),
            ("outputnode.t2w_preproc", "inputnode.t2w_preproc"),
         ]),
        (coregistration_wf, outputnode, [
            ("outputnode.t1w_preproc", "anat_preproc"),
            ("outputnode.t1w_brain", "anat_brain"),
            ("outputnode.t1w_mask", "anat_mask"),
        ]),

        (anat_seg_wf, outputnode, [
            ("outputnode.anat_dseg", "anat_dseg"),
            ("outputnode.anat_tpms", "anat_tpms")]),
        (anat_seg_wf, anat_derivatives_wf, [
            ("outputnode.anat_dseg", "inputnode.t1w_dseg"),
            ("outputnode.anat_tpms", "inputnode.t1w_tpms"),
        ]),
        (anat_seg_wf, anat_norm_wf, [
            ("outputnode.anat_dseg", "inputnode.moving_segmentation"),
            ("outputnode.anat_tpms", "inputnode.moving_tpms")]),

        (anat_norm_wf, anat_reports_wf, [("poutputnode.template", "inputnode.template")]),
        (anat_norm_wf, outputnode, [
            ("poutputnode.standardized", "std_preproc"),
            ("poutputnode.std_mask", "std_mask"),
            ("poutputnode.std_dseg", "std_dseg"),
            ("poutputnode.std_tpms", "std_tpms"),
            ("outputnode.template", "template"),
            ("outputnode.anat2std_xfm", "anat2std_xfm"),
            ("outputnode.std2anat_xfm", "std2anat_xfm")]),
        (anat_norm_wf, anat_derivatives_wf, [
            ("outputnode.template", "inputnode.template"),
            ("outputnode.anat2std_xfm", "inputnode.anat2std_xfm"),
            ("outputnode.std2anat_xfm", "inputnode.std2anat_xfm")]),

        (outputnode, anat_reports_wf, [
            ("anat_preproc", "inputnode.t1w_preproc"),
            ("anat_mask", "inputnode.t1w_mask"),
            ("anat_dseg", "inputnode.t1w_dseg"),
            ("std_preproc", "inputnode.std_t1w"),
            ("std_mask", "inputnode.std_mask"),
        ]),
    ])

    if precomp_mask:
        # Ensure the mask is conformed along with the T1w
        t1w_preproc_wf.inputs.inputnode.in_mask = precomp_mask
        # fmt:off
        wf.connect([
            (t1w_preproc_wf, coregistration_wf, [("outputnode.anat_mask", "inputnode.in_mask")]),
            (t2w_preproc_wf, coregistration_wf, [("outputnode.anat_preproc", "inputnode.in_t2w")])
        ])
        # fmt:on
    else:
        # Run brain extraction on the T2w
        brain_extraction_wf = init_infant_brain_extraction_wf(
            age_months=age_months,
            ants_affine_init=ants_affine_init,
            skull_strip_template=skull_strip_template.space,
            template_specs=skull_strip_template.spec,
            omp_nthreads=omp_nthreads,
            sloppy=sloppy,
            debug="registration" in config.execution.debug,
        )
        # fmt:off
        wf.connect([
            (t2w_preproc_wf, brain_extraction_wf, [
                ("outputnode.anat_preproc", "inputnode.t2w_preproc")]),
            (brain_extraction_wf, coregistration_wf, [
                ("outputnode.t2w_preproc", "inputnode.in_t2w"),
                ("outputnode.out_mask", "inputnode.in_mask"),
                ("outputnode.out_probmap", "inputnode.in_probmap")]),
        ])
        # fmt:on

    if precomp_aseg:
        # Ensure the segmentation is conformed along with the T1w
        t1w_preproc_wf.inputs.inputnode.in_aseg = precomp_aseg
        wf.connect(t1w_preproc_wf, "outputnode.anat_aseg", anat_seg_wf, "inputnode.anat_aseg")

    if not freesurfer:
        return wf

    if config.workflow.surface_recon_method == 'freesurfer':
        from smriprep.workflows.surfaces import init_surface_recon_wf

        surface_recon_wf = init_surface_recon_wf(omp_nthreads=omp_nthreads, hires=hires)
    elif config.workflow.surface_recon_method == 'infantfs':
        from .surfaces import init_infantfs_surface_recon_wf

        # if running with precomputed aseg, or JLF, pass the aseg along to FreeSurfer
        use_aseg = bool(precomp_aseg or segmentation_atlases)
        surface_recon_wf = init_infantfs_surface_recon_wf(
            age_months=age_months,
            use_aseg=use_aseg,
        )

    elif config.workflow.surface_recon_method == 'mcribs':
        from nipype.interfaces.ants import DenoiseImage

        from .surfaces import init_mcribs_surface_recon_wf

        # Denoise raw T2w, since using the template / preproc resulted in intersection errors
        denoise_raw_t2w = pe.Node(
            DenoiseImage(dimension=3, noise_model="Rician"), name='denoise_raw_t2w'
        )

        surface_recon_wf = init_mcribs_surface_recon_wf(
            omp_nthreads=omp_nthreads,
            use_aseg=bool(precomp_aseg),
            use_mask=bool(precomp_mask),
            mcribs_dir=str(config.execution.mcribs_dir),  # Needed to preserve runs
        )

        # Transformed gives
        if precomp_aseg:
            surface_recon_wf.inputs.inputnode.ants_segs = precomp_aseg
        if precomp_mask:
            surface_recon_wf.inputs.inputnode.anat_mask = precomp_mask
        # fmt:off
        wf.connect([
            (inputnode, denoise_raw_t2w, [('t2w', 'input_image')]),
            (denoise_raw_t2w, surface_recon_wf, [('output_image', 'inputnode.t2w')])
        ])
        # fmt:on

    if config.workflow.surface_recon_method in ('freesurfer', 'infantfs'):
        # fmt:off
        wf.connect([
            (t2w_preproc_wf, surface_recon_wf, [
                ("outputnode.anat_preproc", "inputnode.t2w")]),
            (anat_seg_wf, surface_recon_wf, [
                ("outputnode.anat_aseg", "inputnode.ants_segs")]),
        ])
        # fmt:on

    # Anatomical ribbon file using HCP signed-distance volume method
    anat_ribbon_wf = init_anat_ribbon_wf()

    # fmt:off
    wf.connect([
        (inputnode, surface_recon_wf, [
            ("subject_id", "inputnode.subject_id"),
            ("subjects_dir", "inputnode.subjects_dir")]),
        (t1w_template_wf, surface_recon_wf, [
            ("outputnode.anat_ref", "inputnode.t1w"),
        ]),
        (coregistration_wf, surface_recon_wf, [
            ("outputnode.t1w_brain", "inputnode.skullstripped_t1"),
            ("outputnode.t1w_preproc", "inputnode.corrected_t1"),
        ]),
        (surface_recon_wf, outputnode, [
            ("outputnode.subjects_dir", "subjects_dir"),
            ("outputnode.subject_id", "subject_id"),
            ("outputnode.t1w2fsnative_xfm", "t1w2fsnative_xfm"),
            ("outputnode.fsnative2t1w_xfm", "fsnative2t1w_xfm"),
            ("outputnode.surfaces", "surfaces"),
            ("outputnode.morphometrics", "morphometrics"),
            ("outputnode.out_aparc", "anat_aparc"),
            ("outputnode.out_aseg", "anat_aseg"),
        ]),
        (coregistration_wf, anat_ribbon_wf, [
            ("outputnode.t1w_mask", "inputnode.t1w_mask"),
        ]),
        (surface_recon_wf, anat_ribbon_wf, [
            ("outputnode.surfaces", "inputnode.surfaces"),
        ]),
        (anat_ribbon_wf, outputnode, [
            ("outputnode.anat_ribbon", "anat_ribbon")
        ]),
        (anat_ribbon_wf, anat_derivatives_wf, [
            ("outputnode.anat_ribbon", "inputnode.anat_ribbon"),
        ]),
        (surface_recon_wf, anat_reports_wf, [
            ("outputnode.subject_id", "inputnode.subject_id"),
            ("outputnode.subjects_dir", "inputnode.subjects_dir"),
        ]),
        (surface_recon_wf, anat_derivatives_wf, [
            ("outputnode.out_aseg", "inputnode.t1w_fs_aseg"),
            ("outputnode.out_aparc", "inputnode.t1w_fs_aparc"),
            ("outputnode.t1w2fsnative_xfm", "inputnode.t1w2fsnative_xfm"),
            ("outputnode.fsnative2t1w_xfm", "inputnode.fsnative2t1w_xfm"),
            ("outputnode.surfaces", "inputnode.surfaces"),
            ("outputnode.morphometrics", "inputnode.morphometrics"),
        ]),
    ])
    # fmt: on

    if cifti_output:
        from smriprep.workflows.surfaces import init_morph_grayords_wf

        morph_grayords_wf = init_morph_grayords_wf(grayord_density=cifti_output)
        anat_derivatives_wf.get_node('inputnode').inputs.cifti_density = cifti_output
        # fmt:off
        wf.connect([
            (surface_recon_wf, morph_grayords_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
            ]),
            (morph_grayords_wf, anat_derivatives_wf, [
                ("outputnode.cifti_morph", "inputnode.cifti_morph"),
                ("outputnode.cifti_metadata", "inputnode.cifti_metadata"),
            ]),
        ])
        # fmt:on

    return wf
