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

    from smriprep.workflows.anatomical import init_anat_template_wf, _probseg_fast2bids, _pop
    from smriprep.workflows.norm import init_anat_norm_wf
    from smriprep.workflows.outputs import (
        init_anat_reports_wf,
        init_anat_derivatives_wf,
    )

    from .brain_extraction import init_infant_brain_extraction_wf
    from .segmentation import init_anat_seg_wf
    from .surfaces import init_infant_surface_recon_wf

    # for now, T1w only
    num_anats = len(anatomicals)
    wf = pe.Workflow(name=name)

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
    anat_template_wf = init_anat_template_wf(
        longitudial=longitudial,
        omp_nthreads=omp_nthreads,
        num_t1w=num_anat,
    )

    # INU + Brain Extraction
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
    be_buffer = pe.Node(niu.IdentityInterface(fields=["anat_preproc", "anat_brain"]), name='be_buffer')

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

    # FreeSurfer surfaces
    surface_recon_wf = init_infant_surface_recon_wf(age_months=age_months)
    applyrefined = pe.Node(fsl.ApplyMask(), name="applyrefined")

    # fmt: off
    wf.connect([
        (inputnode, brain_extraction_wf, [
            ('t1w', 'inputnode.t1w'),
            ('t2w', 'inputnode.t2w')]),
        (brain_extraction_wf, be_buffer, [
            (('outputnode.out_corrected', _pop), 'anat_preproc'),
            (('outputnode.out_brain', _pop), 'anat_brain'),
            (('outputnode.out_mask', _pop), 'anat_mask')]),
        (be_buffer, outputnode, [
            ('anat_preproc', 'anat_preproc'),
            ('anat_brain', 'anat_brain'),
            ('anat_mask', 'anat_mask')]),
        (be_buffer, anat_seg_wf, [
            ('anat_brain', 'inputnode.anat_brain')]),
        (surface_recon_wf, outputnode, [
            ('outputnode.subjects_dir', 'subjects_dir')]),
        (surface_recon_wf, anat_reports_wf, [
            ('outputnode.subjects_dir', 'subjects_dir'),
            ('outputnode.subject_id', 'subject_id')]),
    ])
    if freesurfer:
        wf.connect([
            (inputnode, surface_recon_wf, [
                ('subject_id', 'inputnode.subject_id'),
                ('subject_dir', 'inputnode.subject_dir')]),
            (be_buffer, surface_recon_wf, [('anat_brain', 'inputnode.masked_file')]),
    # fmt: on
    return wf
