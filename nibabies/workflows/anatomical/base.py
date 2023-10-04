"""Base anatomical preprocessing."""
from __future__ import annotations

import typing as ty
from pathlib import Path

from nipype.interfaces import utility as niu
from nipype.interfaces.ants.base import Info as ANTsInfo
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow
from niworkflows.utils.spaces import Reference, SpatialReferences
from smriprep.workflows.norm import init_anat_norm_wf

from nibabies import config
from nibabies.utils.misc import fix_multi_source_name

# Relative imports to avoid verbosity
from .brain_extraction import init_infant_brain_extraction_wf
from .outputs import (
    init_anat_derivatives_wf,
    init_anat_reports_wf,
    init_coreg_report_wf,
)
from .preproc import init_anat_preproc_wf
from .registration import init_coregister_derivatives_wf, init_coregistration_wf
from .segmentation import init_anat_segmentations_wf
from .surfaces import init_anat_ribbon_wf
from .template import init_anat_template_wf

if ty.TYPE_CHECKING:
    from nibabies.utils.bids import Derivatives


ANAT_OUT_FIELDS = [
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
    "morphometrics",
    "anat_aseg",
    "anat_mcrib",
    "anat_aparc",
    "anat_ribbon",
    "template",
    # registration sphere space is dependent on surface recon method
    "sphere_reg",
    "sphere_reg_fsLR",
    "midthickness_fsLR",
]


def init_infant_anat_wf(
    *,
    age_months: int,
    ants_affine_init: bool,
    t1w: list,
    t2w: list,
    contrast: ty.Literal['T1w', 'T2w'],
    bids_root: str | Path,
    derivatives: Derivatives,
    freesurfer: bool,
    hires: bool | None,
    longitudinal: bool,
    omp_nthreads: int,
    output_dir: str | Path,
    segmentation_atlases: str | Path | None,
    skull_strip_mode: str,
    skull_strip_template: Reference,
    sloppy: bool,
    spaces: SpatialReferences,
    cifti_output: ty.Literal['91k', '170k'] | None,
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
    anat2fsnative_xfm
        LTA-style affine matrix translating from T1w to
        FreeSurfer-conformed subject space
    fsnative2anat_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed
        subject space to T1w
    surfaces
        GIFTI surfaces (gray/white boundary, midthickness, pial, inflated)
    """
    if not t1w or not t2w:
        # Error type?
        raise RuntimeError("Both T1w and T2w images are required to run this workflow.")

    num_t1w = len(t1w)
    num_t2w = len(t2w)
    wf = LiterateWorkflow(name=name)

    # Derivatives used based on the following truth table:
    # |--------|--------|---------------------------------|------------------|
    # | Has T1 | Has T2 | M-CRIB-S surface reconstruction | Derivatives Used |
    # |--------|--------|---------------------------------|------------------|
    # |   Yes  |   No   |             No                  |         T1       |
    # |   Yes  |   Yes  |             No                  |         T1       |
    # |   No   |   Yes  |             No                  |         T2       |
    # |   Yes  |   Yes  |             Yes                 |         T2       |

    recon_method = config.workflow.surface_recon_method
    t1w_mask = bool(derivatives.t1w_mask)
    t1w_aseg = bool(derivatives.t1w_aseg)
    t2w_mask = bool(derivatives.t2w_mask)
    t2w_aseg = bool(derivatives.t2w_aseg)

    # The T2 derivatives are only prioritized first if MCRIBS reconstruction is to be used.
    if recon_method == "mcribs":
        if t2w_aseg:
            t1w_aseg = False
        if t2w_mask:
            t1w_mask = False
    # Otherwise, prioritize T1 derivatives
    if t1w_mask:
        t2w_mask = False
    if t1w_aseg:
        t2w_aseg = False

    config.loggers.workflow.info(
        "Derivatives used:\n\t<T1 mask %s>\n\t<T1 aseg %s>\n\t<T2 mask %s>\n\t<T2 aseg %s>",
        t1w_mask,
        t1w_aseg,
        t2w_mask,
        t2w_aseg,
    )

    desc = _gen_anat_wf_desc(
        t1w=t1w,
        t2w=t2w,
        mask=t1w_mask or t2w_mask,
    )

    wf.__desc__ = desc.format(
        ants_ver=ANTsInfo.version() or "(version unknown)",
        skullstrip_tpl=skull_strip_template.fullname,
    )
    wf.__postdesc__ = ""

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["t1w", "t2w", "subject_id", "subjects_dir"]),  # FLAIR / ROI?
        name="inputnode",
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=ANAT_OUT_FIELDS), name="outputnode")

    # Define output workflows
    anat_reports_wf = init_anat_reports_wf(
        surface_recon=freesurfer, output_dir=output_dir, sloppy=sloppy
    )

    anat_derivatives_wf = init_anat_derivatives_wf(
        bids_root=bids_root,
        surface_recon=freesurfer,
        num_t1w=num_t1w,
        num_t2w=num_t2w,
        output_dir=output_dir,
        spaces=spaces,
        cifti_output=cifti_output,
    )

    t1w_template_wf = init_anat_template_wf(
        contrast="T1w",
        num_files=num_t1w,
        longitudinal=longitudinal,
        omp_nthreads=omp_nthreads,
        sloppy=sloppy,
        has_mask=t1w_mask,
        has_aseg=t1w_aseg,
        name="t1w_template_wf",
    )

    t2w_template_wf = init_anat_template_wf(
        contrast="T2w",
        num_files=num_t2w,
        longitudinal=longitudinal,
        omp_nthreads=omp_nthreads,
        sloppy=sloppy,
        has_mask=t2w_mask,
        has_aseg=t2w_aseg,
        name="t2w_template_wf",
    )

    # Clean up each anatomical template
    # Denoise, INU, + Clipping
    t1w_preproc_wf = init_anat_preproc_wf(name="t1w_preproc_wf")
    t2w_preproc_wf = init_anat_preproc_wf(name="t2w_preproc_wf")

    if skull_strip_mode != "force":
        raise NotImplementedError("Skull stripping is currently required.")

    coregistration_wf = init_coregistration_wf(
        omp_nthreads=omp_nthreads,
        sloppy=sloppy,
        debug="registration" in config.execution.debug,
        t1w_mask=t1w_mask,
        probmap=not t2w_mask,
    )
    coreg_report_wf = init_coreg_report_wf(
        output_dir=output_dir,
    )

    # Segmentation - initial implementation should be simple: JLF
    anat_seg_wf = init_anat_segmentations_wf(
        anat_modality=contrast.capitalize(),  # TODO: Revisit this option
        template_dir=segmentation_atlases,
        sloppy=sloppy,
        omp_nthreads=omp_nthreads,
        precomp_aseg=bool(derivatives.aseg),
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
            ("outputnode.anat_valid_list", "inputnode.t1w_source_files"),
            ("outputnode.anat_realign_xfm", "inputnode.t1w_ref_xfms")]),
        (t1w_template_wf, anat_reports_wf, [
            ("outputnode.out_report", "inputnode.anat_conform_report")]),

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
            ("outputnode.t1w_mask", "inputnode.anat_mask"),
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
            ("outputnode.anat_dseg", "inputnode.anat_dseg"),
            ("outputnode.anat_tpms", "inputnode.anat_tpms"),
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
            ("anat_preproc", "inputnode.anat_preproc"),
            ("anat_mask", "inputnode.anat_mask"),
            ("anat_dseg", "inputnode.anat_dseg"),
            ("std_preproc", "inputnode.std_t1w"),
            ("std_mask", "inputnode.std_mask"),
        ]),
    ])

    # Workflow to move derivatives between T1w/T2w spaces
    # May not be used, but define in case necessary.
    coreg_deriv_wf = init_coregister_derivatives_wf(
        t1w_mask=t1w_mask, t1w_aseg=t1w_aseg, t2w_aseg=t2w_aseg
    )
    deriv_buffer = pe.Node(
        niu.IdentityInterface(fields=['t2w_mask', 't1w_aseg', 't2w_aseg']),
        name='deriv_buffer',
    )
    if derivatives:
        wf.connect([
            (coregistration_wf, coreg_deriv_wf, [('outputnode.t1w2t2w_xfm', 'inputnode.t1w2t2w_xfm')]),
            (t1w_preproc_wf, coreg_deriv_wf, [('outputnode.anat_preproc', 'inputnode.t1w_ref')]),
            (t2w_preproc_wf, coreg_deriv_wf, [('outputnode.anat_preproc', 'inputnode.t2w_ref')]),
        ])

    # Derivative mask is present
    if t1w_mask:
        t1w_template_wf.inputs.inputnode.anat_mask = derivatives.t1w_mask
        t1w_template_wf.inputs.inputnode.mask_reference = derivatives.references['t1w_mask']
        # fmt:off
        wf.connect([
            (t1w_template_wf, coregistration_wf, [('outputnode.anat_mask', 'inputnode.in_mask')]),
            (t2w_preproc_wf, coregistration_wf, [('outputnode.anat_preproc', 'inputnode.in_t2w')]),
            (t1w_template_wf, coreg_deriv_wf, [('outputnode.anat_mask', 'inputnode.t1w_mask')]),
            (coreg_deriv_wf, deriv_buffer, [('outputnode.t2w_mask', 't2w_mask')])
        ])
        # fmt:on
    elif t2w_mask:
        t2w_template_wf.inputs.inputnode.anat_mask = derivatives.t2w_mask
        t2w_template_wf.inputs.inputnode.mask_reference = derivatives.references['t2w_mask']
        # fmt:on
        wf.connect([
            (t2w_template_wf, coregistration_wf, [('outputnode.anat_mask', 'inputnode.in_mask')]),
            (t2w_preproc_wf, coregistration_wf, [('outputnode.anat_preproc', 'inputnode.in_t2w')]),
            (t2w_template_wf, deriv_buffer, [('outputnode.anat_mask', 't2w_mask')]),
        ])
        # fmt:off
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

    # Derivative segmentation is present
    if derivatives.aseg:
        wf.connect(deriv_buffer, 't1w_aseg', anat_seg_wf, 'inputnode.anat_aseg')

        if t1w_aseg:
            t1w_template_wf.inputs.inputnode.anat_aseg = derivatives.t1w_aseg
            t1w_template_wf.inputs.inputnode.aseg_reference = derivatives.references['t1w_aseg']
            # fmt:off
            wf.connect([
                (t1w_template_wf, deriv_buffer, [('outputnode.anat_aseg', 't1w_aseg')]),
                (t1w_template_wf, coreg_deriv_wf, [('outputnode.anat_aseg', 'inputnode.t1w_aseg')]),
                (coreg_deriv_wf, deriv_buffer, [('outputnode.t2w_aseg', 't2w_aseg')]),
            ])
            # fmt:on
        elif t2w_aseg:
            t2w_template_wf.inputs.inputnode.anat_aseg = derivatives.t2w_aseg
            t2w_template_wf.inputs.inputnode.aseg_reference = derivatives.references['t2w_aseg']
            # fmt:off
            wf.connect([
                (t2w_template_wf, deriv_buffer, [('outputnode.anat_aseg', 't2w_aseg')]),
                (t2w_template_wf, coreg_deriv_wf, [('outputnode.anat_aseg', 'inputnode.t2w_aseg')]),
                (coreg_deriv_wf, deriv_buffer, [('outputnode.t1w_aseg', 't1w_aseg')]),
            ])
            # fmt:on

    if not freesurfer:
        return wf

    if recon_method == 'freesurfer':
        from smriprep.workflows.surfaces import init_surface_recon_wf

        surface_recon_wf = init_surface_recon_wf(omp_nthreads=omp_nthreads, hires=hires)
    elif recon_method == 'infantfs':
        from .surfaces import init_infantfs_surface_recon_wf

        # if running with precomputed aseg, or JLF, pass the aseg along to FreeSurfer
        use_aseg = bool(derivatives.aseg or segmentation_atlases)
        surface_recon_wf = init_infantfs_surface_recon_wf(
            age_months=age_months,
            use_aseg=use_aseg,
        )

    elif recon_method == 'mcribs':
        from nipype.interfaces.ants import DenoiseImage

        from .surfaces import init_mcribs_sphere_reg_wf, init_mcribs_surface_recon_wf

        # Denoise template T2w, since using the template / preproc resulted in intersection errors
        denoise_t2w = pe.Node(
            DenoiseImage(dimension=3, noise_model="Rician"), name='denoise_t2w'
        )
        # t2w mask, t2w aseg
        surface_recon_wf = init_mcribs_surface_recon_wf(
            omp_nthreads=omp_nthreads,
            use_aseg=bool(derivatives.aseg),  # TODO: Incorporate mcribs segmentation
            use_mask=bool(derivatives.mask),  # TODO: Pass in mask regardless of derivatives
            mcribs_dir=str(config.execution.mcribs_dir),  # Needed to preserve runs
        )
        # M-CRIB-S to dHCP42week (32k)
        sphere_reg_wf = init_mcribs_sphere_reg_wf()

        # fmt:off
        wf.connect([
            (t2w_template_wf, denoise_t2w, [('outputnode.anat_ref', 'input_image')]),
            (denoise_t2w, surface_recon_wf, [('output_image', 'inputnode.t2w')]),
        ])
        # fmt:on
        if derivatives.aseg:
            wf.connect(deriv_buffer, 't2w_aseg', surface_recon_wf, 'inputnode.ants_segs')
        if derivatives.mask:
            wf.connect(deriv_buffer, 't2w_mask', surface_recon_wf, 'inputnode.anat_mask')
    else:
        raise NotImplementedError

    if recon_method in ('freesurfer', 'infantfs'):
        from smriprep.workflows.surfaces import init_sphere_reg_wf

        # fsaverage to fsLR
        sphere_reg_wf = init_sphere_reg_wf()

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
            ("outputnode.t1w2fsnative_xfm", "anat2fsnative_xfm"),
            ("outputnode.fsnative2t1w_xfm", "fsnative2anat_xfm"),
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
        (surface_recon_wf, sphere_reg_wf, [
            ('outputnode.subject_id', 'inputnode.subject_id'),
            ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
        ]),
        (surface_recon_wf, anat_reports_wf, [
            ("outputnode.subject_id", "inputnode.subject_id"),
            ("outputnode.subjects_dir", "inputnode.subjects_dir"),
        ]),
        (surface_recon_wf, anat_derivatives_wf, [
            ("outputnode.out_aseg", "inputnode.anat_fs_aseg"),
            ("outputnode.out_aparc", "inputnode.anat_fs_aparc"),
            ("outputnode.t1w2fsnative_xfm", "inputnode.anat2fsnative_xfm"),
            ("outputnode.fsnative2t1w_xfm", "inputnode.fsnative2anat_xfm"),
            ("outputnode.surfaces", "inputnode.surfaces"),
            ("outputnode.morphometrics", "inputnode.morphometrics"),
        ]),
        (sphere_reg_wf, outputnode, [
            ('outputnode.sphere_reg', 'sphere_reg'),
            ('outputnode.sphere_reg_fsLR', 'sphere_reg_fsLR')]),
        (sphere_reg_wf, anat_derivatives_wf, [
            ('outputnode.sphere_reg', 'inputnode.sphere_reg'),
            ('outputnode.sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR')]),
    ])
    # fmt: on

    if cifti_output:
        from nibabies.workflows.anatomical.resampling import (
            init_anat_fsLR_resampling_wf,
        )

        is_mcribs = recon_method == "mcribs"
        # handles morph_grayords_wf
        anat_fsLR_resampling_wf = init_anat_fsLR_resampling_wf(cifti_output, mcribs=is_mcribs)
        anat_derivatives_wf.get_node('inputnode').inputs.cifti_density = cifti_output
        # fmt:off
        wf.connect([
            (sphere_reg_wf, anat_fsLR_resampling_wf, [
                ('outputnode.sphere_reg', 'inputnode.sphere_reg'),
                ('outputnode.sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR')]),
            (surface_recon_wf, anat_fsLR_resampling_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                ('outputnode.surfaces', 'inputnode.surfaces'),
                ('outputnode.morphometrics', 'inputnode.morphometrics')]),
            (anat_fsLR_resampling_wf, anat_derivatives_wf, [
                ("outputnode.cifti_morph", "inputnode.cifti_morph"),
                ("outputnode.cifti_metadata", "inputnode.cifti_metadata")]),
            (anat_fsLR_resampling_wf, outputnode, [
                ("outputnode.midthickness_fsLR", "midthickness_fsLR")])
        ])
        # fmt:on

    return wf


def init_infant_single_anat_wf(
    *,
    age_months: int,
    ants_affine_init: bool,
    t1w: list | None,
    t2w: list | None,
    contrast: ty.Literal['T1w', 'T2w'],
    bids_root: str | Path,
    derivatives: Derivatives,
    freesurfer: bool,
    hires: bool | None,
    longitudinal: bool,
    omp_nthreads: int,
    output_dir: str | Path,
    segmentation_atlases: str | Path | None,
    skull_strip_mode: str,
    skull_strip_template: Reference,
    sloppy: bool,
    spaces: SpatialReferences,
    cifti_output: ty.Literal['91k', '170k'] | None,
    name: str = "infant_single_anat_wf",
) -> LiterateWorkflow:
    """"""
    if t1w and t2w:
        # Error type?
        raise RuntimeError(
            "This workflow uses only T1w or T2w inputs, but both contrasts are available."
        )

    anat_files = t1w or t2w
    num_files = len(anat_files)
    workflow = LiterateWorkflow(name=name)

    # Precomputed derivatives
    if contrast == 'T1w':
        mask = derivatives.t1w_mask
        aseg = derivatives.t1w_aseg
    elif contrast == 'T2w':
        mask = derivatives.t2w_mask
        aseg = derivatives.t2w_aseg

    config.loggers.workflow.info(
        f"Derivatives used (%s):\n\t\t<Mask: %s>\n\t\t<Aseg: %s>\n",
        contrast,
        bool(mask),
        bool(aseg),
    )

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["t1w", "t2w", "subject_id", "subjects_dir"]),  # FLAIR / ROI?
        name="inputnode",
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=ANAT_OUT_FIELDS), name="outputnode")

    desc = _gen_anat_wf_desc(
        t1w=t1w or None,
        t2w=t2w or None,
        mask=bool(mask),
    )
    workflow.__desc__ = desc.format(
        ants_ver=ANTsInfo.version() or "(version unknown)",
        skullstrip_tpl=skull_strip_template.fullname,
    )
    workflow.__postdesc__ = ""

    # outputs
    recon_method = config.workflow.surface_recon_method  # TODO: Make workflow parameter
    anat_reports_wf = init_anat_reports_wf(
        surface_recon=recon_method, output_dir=output_dir, sloppy=sloppy
    )

    # TODO: Update transforms TO-FROM to reflect contrast
    anat_derivatives_wf = init_anat_derivatives_wf(
        bids_root=bids_root,
        output_dir=output_dir,
        surface_recon=recon_method,
        num_t1w=num_files if contrast == 'T1w' else None,
        num_t2w=num_files if contrast == 'T2w' else None,
        spaces=spaces,
        cifti_output=bool(cifti_output),
    )

    # template
    anat_template_wf = init_anat_template_wf(
        contrast=contrast,
        num_files=num_files,
        longitudinal=longitudinal,
        omp_nthreads=omp_nthreads,
        sloppy=sloppy,
        has_mask=bool(mask),
        has_aseg=bool(aseg),
        name=f"{contrast.lower()}_template_wf",
    )
    # preproc
    anat_preproc_wf = init_anat_preproc_wf(name=f"{contrast.lower()}_preproc_wf")
    # T2-only brain extraction
    anat_seg_wf = init_anat_segmentations_wf(
        anat_modality=contrast,
        template_dir=segmentation_atlases,
        sloppy=sloppy,
        omp_nthreads=omp_nthreads,
        precomp_aseg=bool(aseg),
    )
    # T2-only segmentation
    anat_norm_wf = init_anat_norm_wf(
        sloppy=sloppy,
        omp_nthreads=omp_nthreads,
        templates=spaces.get_spaces(nonstandard=False, dim=(3,)),
    )
    # Aggregate mask, applied mask
    mask_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_mask', 'anat_brain']),
        name='mask_buffer',
    )

    if mask:
        from niworkflows.interfaces.nibabel import ApplyMask

        anat_template_wf.inputs.inputnode.anat_mask = mask
        mask_ref = derivatives.references[f'{contrast.lower()}_mask']
        anat_template_wf.inputs.inputnode.mask_reference = mask_ref

        apply_deriv_mask = pe.Node(ApplyMask(), name='apply_deriv_mask')
        # fmt:off
        workflow.connect([
            (anat_template_wf, mask_buffer, [
                ('outputnode.anat_mask', 'anat_mask')]),
            (anat_preproc_wf, apply_deriv_mask, [
                ('outputnode.anat_preproc', 'in_file')]),
            (anat_template_wf, apply_deriv_mask, [
                ('outputnode.anat_mask', 'in_mask')]),
            (apply_deriv_mask, mask_buffer, [
                ('out_file', 'anat_brain')]),
        ])
        # fmt:on

    else:
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
        workflow.connect([
            (anat_preproc_wf, brain_extraction_wf, [
                ('outputnode.anat_preproc', 'inputnode.t2w_preproc')]),
            (brain_extraction_wf, mask_buffer, [
                ('outputnode.t2w_brain', 'anat_brain'),
                ('outputnode.out_mask', 'anat_mask')]),
        ])
        # fmt:on

    if aseg:
        anat_template_wf.inputs.inputnode.anat_aseg = aseg
        aseg_ref = derivatives.references[f'{contrast.lower()}_aseg']
        anat_template_wf.inputs.inputnode.aseg_reference = aseg_ref

        workflow.connect(
            anat_template_wf, 'outputnode.anat_aseg', anat_seg_wf, 'inputnode.anat_aseg'
        )

    # fmt:off
    workflow.connect([
        (inputnode, anat_template_wf, [(contrast.lower(), "inputnode.anat_files")]),
        (inputnode, anat_reports_wf, [(contrast.lower(), "inputnode.source_file")]),
        (inputnode, anat_norm_wf, [((contrast.lower(), fix_multi_source_name), "inputnode.orig_t1w")]),

        (anat_template_wf, outputnode, [
            ("outputnode.anat_realign_xfm", "anat_ref_xfms")]),
        (anat_template_wf, anat_preproc_wf, [
            ("outputnode.anat_ref", "inputnode.in_anat")]),
        (anat_template_wf, anat_derivatives_wf, [
            ("outputnode.anat_valid_list", f"inputnode.{contrast.lower()}_source_files"),
            ("outputnode.anat_realign_xfm", f"inputnode.{contrast.lower()}_ref_xfms")]),
        (anat_template_wf, anat_reports_wf, [
            ("outputnode.out_report", "inputnode.anat_conform_report")]),
        (anat_preproc_wf, anat_norm_wf, [
            ('outputnode.anat_preproc', 'inputnode.moving_image')]),
        (anat_preproc_wf, outputnode, [
            ('outputnode.anat_preproc', 'anat_preproc')]),
        (anat_preproc_wf, anat_derivatives_wf, [
            ('outputnode.anat_preproc', f'inputnode.{contrast.lower()}_preproc')]),
        (mask_buffer, anat_derivatives_wf, [
            ('anat_mask', 'inputnode.anat_mask')]),
        (mask_buffer, outputnode, [
            ('anat_mask', 'anat_mask')]),
        (mask_buffer, anat_seg_wf, [('anat_brain', 'inputnode.anat_brain')]),
        (anat_seg_wf, outputnode, [
            ("outputnode.anat_dseg", "anat_dseg"),
            ("outputnode.anat_tpms", "anat_tpms")]),
        (anat_seg_wf, anat_derivatives_wf, [
            ("outputnode.anat_dseg", "inputnode.anat_dseg"),
            ("outputnode.anat_tpms", "inputnode.anat_tpms")]),
        (mask_buffer, anat_norm_wf, [
            ('anat_mask', 'inputnode.moving_mask')]),
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
            ("anat_preproc", "inputnode.anat_preproc"),
            ("anat_mask", "inputnode.anat_mask"),
            ("anat_dseg", "inputnode.anat_dseg"),
            ("std_preproc", "inputnode.std_t1w"),
            ("std_mask", "inputnode.std_mask"),
        ]),
    ])
    # fmt:on
    # TODO: Remove `freesurfer` option
    if not recon_method:
        return workflow

    elif recon_method == 'freesurfer':
        from smriprep.workflows.surfaces import init_surface_recon_wf

        surface_recon_wf = init_surface_recon_wf(omp_nthreads=omp_nthreads, hires=hires)
    elif recon_method == 'infantfs':
        from .surfaces import init_infantfs_surface_recon_wf

        # if running with precomputed aseg, or JLF, pass the aseg along to FreeSurfer
        use_aseg = bool(derivatives.aseg or segmentation_atlases)
        surface_recon_wf = init_infantfs_surface_recon_wf(
            age_months=age_months,
            use_aseg=use_aseg,
        )

    elif recon_method == 'mcribs':
        from .surfaces import init_mcribs_sphere_reg_wf, init_mcribs_surface_recon_wf

        # t2w mask, t2w aseg
        surface_recon_wf = init_mcribs_surface_recon_wf(
            omp_nthreads=omp_nthreads,
            use_aseg=bool(aseg),  # TODO: Incorporate mcribs segmentation
            use_mask=bool(mask),  # TODO: Pass in mask regardless of derivatives
            mcribs_dir=str(config.execution.mcribs_dir),  # Needed to preserve runs
        )
        # M-CRIB-S to dHCP42week (32k)
        sphere_reg_wf = init_mcribs_sphere_reg_wf()

        # fmt:off
        workflow.connect([
            (anat_preproc_wf, surface_recon_wf, [('outputnode.anat_preproc', 'inputnode.t2w')]),
        ])
        # fmt:on
        if aseg:
            workflow.connect(
                anat_template_wf, 'outputnode.anat_aseg', surface_recon_wf, 'inputnode.ants_segs'
            )
        else:
            # TODO: Use MCRIBS segmentation
            ...
    else:
        raise NotImplementedError

    # Anatomical ribbon file using HCP signed-distance volume method
    anat_ribbon_wf = init_anat_ribbon_wf()

    # fmt:off
    workflow.connect([
        (inputnode, surface_recon_wf, [
            ("subject_id", "inputnode.subject_id"),
            ("subjects_dir", "inputnode.subjects_dir")]),
        (anat_template_wf, surface_recon_wf, [
            ("outputnode.anat_ref", "inputnode.t1w"),
        ]),
        (mask_buffer, surface_recon_wf, [
            ("anat_brain", "inputnode.skullstripped_t1"),
            ("anat_mask", "inputnode.anat_mask")]),
        (anat_preproc_wf, surface_recon_wf, [
            ("outputnode.anat_preproc", "inputnode.corrected_t1")]),
        (surface_recon_wf, outputnode, [
            ("outputnode.subjects_dir", "subjects_dir"),
            ("outputnode.subject_id", "subject_id"),
            ("outputnode.t1w2fsnative_xfm", "anat2fsnative_xfm"),
            ("outputnode.fsnative2t1w_xfm", "fsnative2anat_xfm"),
            ("outputnode.surfaces", "surfaces"),
            ("outputnode.morphometrics", "morphometrics"),
            ("outputnode.out_aparc", "anat_aparc"),
            ("outputnode.out_aseg", "anat_aseg"),
        ]),
        (mask_buffer, anat_ribbon_wf, [
            ("anat_mask", "inputnode.t1w_mask"),
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
        (surface_recon_wf, sphere_reg_wf, [
            ('outputnode.subject_id', 'inputnode.subject_id'),
            ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
        ]),
        (surface_recon_wf, anat_reports_wf, [
            ("outputnode.subject_id", "inputnode.subject_id"),
            ("outputnode.subjects_dir", "inputnode.subjects_dir"),
        ]),
        (surface_recon_wf, anat_derivatives_wf, [
            ("outputnode.out_aseg", "inputnode.anat_fs_aseg"),
            ("outputnode.out_aparc", "inputnode.anat_fs_aparc"),
            ("outputnode.t1w2fsnative_xfm", "inputnode.anat2fsnative_xfm"),
            ("outputnode.fsnative2t1w_xfm", "inputnode.fsnative2anat_xfm"),
            ("outputnode.surfaces", "inputnode.surfaces"),
            ("outputnode.morphometrics", "inputnode.morphometrics"),
        ]),
        (sphere_reg_wf, outputnode, [
            ('outputnode.sphere_reg', 'sphere_reg'),
            ('outputnode.sphere_reg_fsLR', 'sphere_reg_fsLR')]),
        (sphere_reg_wf, anat_derivatives_wf, [
            ('outputnode.sphere_reg', 'inputnode.sphere_reg'),
            ('outputnode.sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR')]),
    ])
    # fmt: on

    if cifti_output:
        from nibabies.workflows.anatomical.resampling import (
            init_anat_fsLR_resampling_wf,
        )

        is_mcribs = recon_method == "mcribs"
        # handles morph_grayords_wf
        anat_fsLR_resampling_wf = init_anat_fsLR_resampling_wf(cifti_output, mcribs=is_mcribs)
        anat_derivatives_wf.get_node('inputnode').inputs.cifti_density = cifti_output
        # fmt:off
        workflow.connect([
            (sphere_reg_wf, anat_fsLR_resampling_wf, [
                ('outputnode.sphere_reg', 'inputnode.sphere_reg'),
                ('outputnode.sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR')]),
            (surface_recon_wf, anat_fsLR_resampling_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                ('outputnode.surfaces', 'inputnode.surfaces'),
                ('outputnode.morphometrics', 'inputnode.morphometrics')]),
            (anat_fsLR_resampling_wf, anat_derivatives_wf, [
                ("outputnode.cifti_morph", "inputnode.cifti_morph"),
                ("outputnode.cifti_metadata", "inputnode.cifti_metadata")]),
            (anat_fsLR_resampling_wf, outputnode, [
                ("outputnode.midthickness_fsLR", "midthickness_fsLR")])
        ])
        # fmt:on
    return workflow


def _gen_anat_wf_desc(t1w: list | None, t2w: list | None, mask: bool) -> str:
    """Generate the anatomical workflow description."""
    if not t1w and not t2w:
        return ''

    # If only a single anatomical modality is provided
    modality = None
    anat = None
    if not t1w or not t2w:
        anat = t1w or t2w
        modality = 'T1w' if t1w else 'T2w'

    desc = """\n\nAnatomical data preprocessing\n:"""

    # Anatomicals found
    if anat is not None:
        desc += (
            f"A total of {len(anat)} {modality} images were found "
            "within the input BIDS dataset.\n"
        )
    else:
        desc += (
            f"A total of {len(t1w)} T1w and {len(t2w)} T2w images "
            "were found within the input BIDS dataset.\n"
        )

    # Template + Preproc workflows
    if t1w:
        if len(t1w) == 1:
            desc += (
                f"The T1-weighted (T1w) image was denoised "
                "and corrected for intensity non-uniformity (INU)"
            )
        else:
            desc += (
                "All of the T1-weighted images were corrected for intensity "
                "non-uniformity (INU)"
            )
        desc += (
            "with `N4BiasFieldCorrection` [@n4], distributed with ANTs {ants_ver} "
            "[@ants, RRID:SCR_004757]"
        )
        desc += ".\n" if len(t1w) > 1 else ", and used as T1w-reference throughout the workflow.\n"
    if t2w:
        if len(t2w) == 1:
            desc += (
                "The T2-weighted (T2w) image was denoised and corrected for intensity "
                "non-uniformity (INU)"
            )
        else:
            desc += (
                "All of the T2-weighted images were corrected for intensity "
                "non-uniformity (INU)"
            )
        desc += (
            "with `N4BiasFieldCorrection` [@n4], distributed with ANTs {ants_ver} "
            "[@ants, RRID:SCR_004757]"
        )
        desc += ".\n" if len(t2w) > 1 else ", and used as T2w-reference throughout the workflow.\n"

    # Precomputed derivatives
    if mask:
        desc += "A previously computed mask was used to skull-strip the anatomical image."

    else:
        desc += (
            f"The {modality or 'T2w'}-reference was then skull-stripped with a modified "
            "implementation of the `antsBrainExtraction.sh` workflow (from ANTs), using "
            "{skullstrip_tpl} as target template."
        )

    return desc
