from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl


def init_infant_anat_wf(
    *,
    age_months,
    anatomicals,
    anat_modality,
    bids_root,
    existing_derivatives,
    freesurfer,
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
    from niworkflows.interfaces.images import ValidateImage
    from smriprep.workflows.anatomical import init_anat_template_wf, _probseg_fast2bids, _pop
    from smriprep.workflows.norm import init_anat_norm_wf
    from smriprep.workflows.outputs import (
        init_anat_reports_wf,
        init_anat_derivatives_wf,
    )

    from ...utils.misc import fix_multi_source_name
    from .brain_extraction import init_infant_brain_extraction_wf
    from .segmentation import init_anat_seg_wf
    from .surfaces import init_infant_surface_recon_wf

    # for now, T1w only
    num_anats = len(anatomicals)
    wf = pe.Workflow(name=name)
    desc = """Anatomical data preprocessing

: """
    desc += f"""\
A total of {num_anats} anatomical images were found within the input
BIDS dataset."""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=["t1w", "t2w", "subject_id", "subjects_dir"]
        ),  # FLAIR / ROI?
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
            ]
        ),
        name="outputnode",
    )

    # Connect reportlets workflows
    anat_reports_wf = init_anat_reports_wf(
        freesurfer=freesurfer,
        output_dir=output_dir,
    )

    if existing_derivatives:
        raise NotImplementedError("Reusing derivatives is not yet supported.")

    desc += """
All of them were corrected for intensity non-uniformity (INU)
""" if num_anats > 1 else """\
The T1-weighted (T1w) image was corrected for intensity non-uniformity (INU)
"""
    desc += """\
with `N4BiasFieldCorrection` [@n4], distributed with ANTs {ants_ver} \
[@ants, RRID:SCR_004757]"""
    desc += '.\n' if num_anats > 1 else ", and used as T1w-reference throughout the workflow.\n"

    desc += """\
The T1w-reference was then skull-stripped with a modified implementation of
the `antsBrainExtraction.sh` workflow (from ANTs), using {skullstrip_tpl}
as target template.
Brain tissue segmentation of cerebrospinal fluid (CSF),
white-matter (WM) and gray-matter (GM) was performed on
the brain-extracted T1w using ANTs JointFusion, distributed with ANTs {ants_ver}.
"""

    workflow.__desc__ = desc.format(
        ants_ver=ANTsInfo.version() or '(version unknown)',
        skullstrip_tpl=skull_strip_template.fullname,
    )
    # Define output workflows
    anat_reports_wf = init_anat_reports_wf(freesurfer=freesurfer, output_dir=output_dir)
    anat_derivatives_wf = init_anat_derivatives_wf(
        bids_root=bids_root,
        freesurfer=freesurfer,
        num_t1w=num_anats,
        output_dir=output_dir,
        spaces=spaces,
    )

    # Multiple T1w files -> generate average reference
    # TODO: Add path for T2w
    anat_template_wf = init_anat_template_wf(
        longitudial=longitudial,
        omp_nthreads=omp_nthreads,
        num_t1w=num_anat,
    )

    anat_validate = pe.Node(
        ValidateImage(),
        name='anat_validate',
        run_without_submitting=True,
    )

    # INU + Brain Extraction
    if skull_strip_mode != 'force':
        raise NotImplementedError("Skull stripping is currently required.")

    brain_extraction_wf = init_infant_brain_extraction_wf(
        age_months=age_months,
        anat_modality=anat_modality,
        ants_affine_init=True,
        skull_strip_template=skull_strip_template,
        template_specs=template_specs,
        omp_nthreads=omp_nthreads,
        output_dir=output_dir,
        sloppy=sloppy,
    )
    # Ensure single outputs
    be_buffer = pe.Node(
        niu.IdentityInterface(
            fields=["anat_preproc", "anat_brain"]
        ),
        name='be_buffer'
    )

    # Segmentation - initial implementation should be simple: JLF
    anat_seg_wf = init_anat_seg_wf(
        age_months=age_months,
        anat_modality=anat_modality,
        template_dir=segmentation_atlases,
        sloppy=sloppy,
        omp_nthreads=omp_nthreads,
    )

    # Spatial normalization (requires segmentation)
    anat_norm_wf = init_anat_norm_wf(
        debug=debug,
        omp_nthreads=omp_nthreads,
        templates=spaces.get_spaces(nonstandard=False, dim=(3,)),
    )

    # fmt: off
    wf.connect([
        (inputnode, anat_template_wf, [
            ('t1w', 'inputnode.t1w'),
        ]),
        (anat_template_wf, outputnode, [
            ('outputnode.t1w_realign_xfm', 'anat_ref_xfms'),
        ]),
        (anat_template_wf, anat_validate, [
            ('outputnode.t1w_ref', 'in_file'),
        ]),
        (anat_validate, brain_extraction_wf, [
            ('out_file', 'inputnode.in_file'),
        ]),
        (brain_extraction_wf, be_buffer, [
            (('outputnode.out_corrected', _pop), 'anat_preproc'),
            (('outputnode.out_brain', _pop), 'anat_brain'),
            (('outputnode.out_mask', _pop), 'anat_mask'),
        ]),
        (be_buffer, outputnode, [
            ('anat_preproc', 'anat_preproc'),
            ('anat_brain', 'anat_brain'),
            ('anat_mask', 'anat_mask'),
        ]),
        (be_buffer, anat_seg_wf, [
            ('anat_brain', 'inputnode.anat_brain'),
        ]),
        (anat_seg_wf, outputnode, [
            ('outputnode.anat_dseg', 'anat_dseg'),
        ]),
        (anat_seg_wf, anat_norm_wf, [
            ('outputnode.anat_dseg', 'inputnode.moving_segmentation'),
        ]),
        (anat_seg_wf, anat_derivatives_wf, [
            ('outputnode.anat_aseg', 'inputnode.t1w_fs_aseg'),
        ]),
        (be_buffer, anat_norm_wf, [
            ('anat_preproc', 'inputnode.moving_image'),
            ('anat_mask', 'inputnode.moving_mask'),
        ]),
        (anat_norm_wf, outputnode, [
            ('poutputnode.standardized', 'std_preproc'),
            ('poutputnode.std_mask', 'std_mask'),
            ('poutputnode.std_dseg', 'std_dseg'),
            ('poutputnode.std_tpms', 'std_tpms'),
            ('outputnode.template', 'template'),
            ('outputnode.anat2std_xfm', 'anat2std_xfm'),
            ('outputnode.std2anat_xfm', 'std2anat_xfm'),
        ]),
        (inputnode, anat_norm_wf, [
            (('t1w', fix_multi_source_name), 'inputnode.orig_t1w'),  # anat_validate? not used...
        ]),
    ])

    wf.connect([
        # reports
        (inputnode, anat_reports_wf, [
            ('t1w', 'inputnode.source_file'),
        ]),
        (outputnode, anat_reports_wf, [
            ('anat_preproc', 'inputnode.t1w_preproc'),
            ('anat_mask', 'inputnode.t1w_mask'),
            ('anat_dseg', 'inputnode.t1w_dseg'),
            ('std_preproc', 'inputnode.std_t1w'),
            ('std_mask', 'inputnode.std_mask'),
        ]),
        (anat_template_wf, anat_reports_wf, [
            ('outputnode.out_report', 'inputnode.t1w_conform_report'),
        ]),
        (anat_norm_wf, anat_reports_wf, [
            ('poutputnode.template', 'inputnode.template'),
        ]),
        # derivatives
        (anat_template_wf, anat_derivatives_wf, [
            ('outputnode.t1w_valid_list', 'inputnode.source_files'),
        ]),
        (anat_norm_wf, anat_derivatives_wf, [
            ('outputnode.template', 'inputnode.template'),
            ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
            ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
        ]),
        (outputnode, anat_derivatives_wf, [
            ('anat_ref_xfms', 'inputnode.t1w_ref_xfms'),
            ('anat_preproc', 'inputnode.t1w_preproc'),
            ('anat_mask', 'inputnode.t1w_mask'),
            ('anat_dseg', 'inputnode.t1w_dseg'),
            ('anat_tpms', 'inputnode.t1w_tpms'),
        ]),
    ])

    if not freesurfer:
        return wf

    # FreeSurfer surfaces
    surface_recon_wf = init_infant_surface_recon_wf(age_months=age_months)

    wf.connect([
        (inputnode, surface_recon_wf, [
            ('subject_id', 'inputnode.subject_id'),
            ('subject_dir', 'inputnode.subject_dir'),
            ('t2w', 'inputnode.t2w'),
        ]),
        (anat_validate, surface_recon_wf, [
            ('out_file', 'inputnode.anat_orig'),
        ]),
        (be_buffer, surface_recon_wf, [
            ('anat_brain', 'inputnode.anat_skullstripped'),
            ('anat_preproc', 'inputnode.anat_preproc'),
        ]),
        (surface_recon_wf, outputnode, [
            ('outputnode.subjects_dir', 'subjects_dir'),
            ('outputnode.subject_id', 'subject_id'),
            ('outputnode.t1w2fsnative_xfm', 't1w2fsnative_xfm'),
            ('outputnode.fsnative2t1w_xfm', 'fsnative2t1w_xfm'),
            ('outputnode.surfaces', 'surfaces'),
            ('outputnode.out_aseg', 't1w_aseg'),
            ('outputnode.out_aparc', 't1w_aparc'),
        ]),
        (surface_recon_wf, anat_reports_wf, [
            ('outputnode.subject_id', 'inputnode.subject_id'),
            ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
        ]),
        (surface_recon_wf, anat_derivatives_wf, [
            ('outputnode.out_aparc', 'inputnode.t1w_fs_aparc'),
        ]),
        (surface_recon_wf, anat_derivatives_wf, [
            ('outputnode.t1w2fsnative_xfm', 'inputnode.t1w2fsnative_xfm'),
            ('outputnode.fsnative2t1w_xfm', 'inputnode.fsnative2t1w_xfm'),
            ('outputnode.surfaces', 'inputnode.surfaces'),
        ]),
    ])
    # fmt: on
    return wf
