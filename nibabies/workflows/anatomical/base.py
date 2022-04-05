"""Base anatomical preprocessing."""
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from ... import config


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
    from nipype.interfaces.ants.base import Info as ANTsInfo
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from ...utils.misc import fix_multi_source_name
    from .brain_extraction import (
        init_infant_brain_extraction_wf,
        init_precomputed_mask_wf,
    )
    from .norm import init_anat_norm_wf
    from .outputs import (
        init_anat_derivatives_wf,
        init_anat_reports_wf,
        init_coreg_report_wf,
    )
    from .preproc import init_anat_average_wf
    from .registration import init_coregistration_wf
    from .segmentation import init_anat_segmentations_wf
    from .surfaces import init_infant_surface_recon_wf

    # for now, T1w only
    num_t1w = len(t1w) if t1w else 0
    num_t2w = len(t2w) if t2w else 0

    precomp_mask = existing_derivatives.get("anat_mask")
    precomp_aseg = existing_derivatives.get("anat_aseg")

    # verify derivatives are relatively similar to T1w
    if precomp_mask or precomp_aseg:
        from ...utils.validation import validate_t1w_derivatives

        validated_derivatives = validate_t1w_derivatives(  # compare derivatives to the first T1w
            t1w[0], anat_mask=precomp_mask, anat_aseg=precomp_aseg
        )
        precomp_mask = validated_derivatives.get("anat_mask")
        precomp_aseg = validated_derivatives.get("anat_aseg")

    wf = Workflow(name=name)
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
                "anat2fsnative_xfm",
                "fsnative2anat_xfm",
                "surfaces",
                "anat_aseg",
                "anat_aparc",
                "template",
            ]
        ),
        name="outputnode",
    )

    desc += (
        """\
All of the T1-weighted images were corrected for intensity non-uniformity (INU)"""
        if num_t1w > 1
        else """\
The T1-weighted (T1w) image was corrected for intensity non-uniformity (INU)"""
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
    )

    # Multiple T1w files -> generate average reference
    t1w_template_wf = init_anat_average_wf(
        longitudinal=longitudinal,
        omp_nthreads=omp_nthreads,
        num_maps=num_t1w,
        name="t1w_template_wf",
    )

    t2w_template_wf = init_anat_average_wf(
        longitudinal=longitudinal,
        omp_nthreads=omp_nthreads,
        num_maps=num_t2w,
        name="t2w_template_wf",
    )

    # INU + Brain Extraction
    if skull_strip_mode != "force":
        raise NotImplementedError("Skull stripping is currently required.")

    if precomp_mask:
        precomp_mask_wf = init_precomputed_mask_wf(omp_nthreads=omp_nthreads)
        precomp_mask_wf.inputs.inputnode.t1w_mask = precomp_mask

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
        coregistration_wf = init_coregistration_wf(
            omp_nthreads=omp_nthreads,
            sloppy=sloppy,
            debug="registration" in config.execution.debug,
        )
        coreg_report_wf = init_coreg_report_wf(
            output_dir=output_dir,
        )
    t1w_preproc_wf = precomp_mask_wf if precomp_mask else coregistration_wf

    # Segmentation - initial implementation should be simple: JLF
    anat_seg_wf = init_anat_segmentations_wf(
        anat_modality=anat_modality.capitalize(),
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
        (inputnode, t1w_template_wf, [("t1w", "inputnode.in_files")]),
        (t1w_template_wf, outputnode, [
            ("outputnode.realign_xfms", "anat_ref_xfms"),
        ]),
        (anat_seg_wf, outputnode, [
            ("outputnode.anat_dseg", "anat_dseg"),
            ("outputnode.anat_tpms", "anat_tpms"),
        ]),
        (anat_seg_wf, anat_norm_wf, [
            ("outputnode.anat_dseg", "inputnode.moving_segmentation"),
            ("outputnode.anat_tpms", "inputnode.moving_tpms"),
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
        (t1w_preproc_wf, anat_norm_wf, [
            ("outputnode.t1w_preproc", "inputnode.moving_image"),
            ("outputnode.t1w_mask", "inputnode.moving_mask"),
        ]),
        (t1w_preproc_wf, anat_derivatives_wf, [
            ("outputnode.t1w_mask", "inputnode.t1w_mask"),
            ("outputnode.t1w_preproc", "inputnode.t1w_preproc"),
        ]),
        (t1w_preproc_wf, outputnode, [
            ("outputnode.t1w_preproc", "anat_preproc"),
            ("outputnode.t1w_brain", "anat_brain"),
            ("outputnode.t1w_mask", "anat_mask"),
        ]),
    ])

    if not precomp_aseg:
        wf.connect([
            (t1w_preproc_wf, anat_seg_wf, [("outputnode.t1w_brain", "inputnode.anat_brain")]),
        ])

    if precomp_mask:
        wf.connect([
            (t1w_template_wf, precomp_mask_wf, [
                ("outputnode.out_file", "inputnode.t1w"),
            ]),
        ])
    else:
        wf.connect([
            (inputnode, t2w_template_wf, [("t2w", "inputnode.in_files")]),
            (t2w_template_wf, brain_extraction_wf, [
                ("outputnode.out_file", "inputnode.in_t2w"),
            ]),
            (t1w_template_wf, coregistration_wf, [
                ("outputnode.out_file", "inputnode.in_t1w"),
            ]),
            (brain_extraction_wf, coregistration_wf, [
                ("outputnode.t2w_preproc", "inputnode.in_t2w_preproc"),
                ("outputnode.out_mask", "inputnode.in_mask"),
                ("outputnode.out_probmap", "inputnode.in_probmap"),
            ]),
            (inputnode, coreg_report_wf, [
                ("t1w", "inputnode.source_file"),
            ]),
            (t1w_preproc_wf, coreg_report_wf, [
                ("outputnode.t1w_preproc", "inputnode.t1w_preproc"),
                ("outputnode.t2w_preproc", "inputnode.t2w_preproc"),
                ("outputnode.t1w_mask", "inputnode.in_mask"),
            ]),
        ])

    wf.connect([
        # reports
        (inputnode, anat_reports_wf, [
            ("t1w", "inputnode.source_file"),
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
            ("outputnode.valid_list", "inputnode.source_files"),
            ("outputnode.realign_xfms", "inputnode.t1w_ref_xfms"),
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

    # if running with precomputed aseg, or JLF, pass the aseg along to FreeSurfer
    use_aseg = bool(precomp_aseg) or bool(segmentation_atlases)
    # FreeSurfer surfaces
    surface_recon_wf = init_infant_surface_recon_wf(
        age_months=age_months,
        use_aseg=use_aseg,
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
        (t1w_template_wf, surface_recon_wf, [
            ("outputnode.out_file", "inputnode.anat_orig"),
        ]),
        (t1w_preproc_wf, surface_recon_wf, [
            ("outputnode.t1w_brain", "inputnode.anat_skullstripped"),
            ("outputnode.t1w_preproc", "inputnode.anat_preproc"),
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
