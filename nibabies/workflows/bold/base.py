# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Orchestrating the BOLD-preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_func_preproc_wf
.. autofunction:: init_func_derivatives_wf

"""
from ... import config

import os

import nibabel as nb
from nipype.interfaces.fsl import Split as FSLSplit
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from niworkflows.utils.connections import pop_file, listify
from niworkflows.interfaces.header import ValidateImage
from niworkflows.interfaces.utility import KeySelect
from sdcflows.interfaces.brainmask import BrainExtraction

from ...interfaces import DerivativesDataSink
from ...interfaces.reports import FunctionalSummary
from ...utils.bids import extract_entities
from ...utils.misc import combine_meepi_source

# BOLD workflows
from .confounds import init_bold_confs_wf, init_carpetplot_wf
from .hmc import init_bold_hmc_wf
from .stc import init_bold_stc_wf
from .t2s import init_bold_t2s_wf
from .registration import init_bold_t1_trans_wf, init_bold_reg_wf
from .resampling import (
    init_bold_surf_wf,
    init_bold_std_trans_wf,
    init_bold_preproc_trans_wf,
)
from .outputs import init_func_derivatives_wf


def init_func_preproc_wf(bold_file, has_fieldmap=False):
    """
    This workflow controls the functional preprocessing stages of *fMRIPrep*.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.tests import mock_config
            from fmriprep import config
            from fmriprep.workflows.bold.base import init_func_preproc_wf
            with mock_config():
                bold_file = config.execution.bids_dir / 'sub-01' / 'func' \
                    / 'sub-01_task-mixedgamblestask_run-01_bold.nii.gz'
                wf = init_func_preproc_wf(str(bold_file))

    Parameters
    ----------
    bold_file
        BOLD series NIfTI file
    has_fieldmap
        Signals the workflow to use inputnode fieldmap files

    Inputs
    ------
    bold_file
        BOLD series NIfTI file
    t1w_preproc
        Bias-corrected structural template image
    t1w_mask
        Mask of the skull-stripped template image
    t1w_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    t1w_asec
        Segmentation of structural image, done with FreeSurfer.
    t1w_aparc
        Parcellation of structural image, done with FreeSurfer.
    t1w_tpms
        List of tissue probability maps in T1w space
    template
        List of templates to target
    anat2std_xfm
        List of transform files, collated with templates
    std2anat_xfm
        List of inverse transform files, collated with templates
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    t1w2fsnative_xfm
        LTA-style affine matrix translating from T1w to FreeSurfer-conformed subject space
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed subject space to T1w
    bold_ref
        BOLD reference file
    bold_ref_xfm
        Transform file in LTA format from bold to reference
    n_dummy_scans
        Number of nonsteady states at the beginning of the BOLD run

    Outputs
    -------
    bold_t1
        BOLD series, resampled to T1w space
    bold_mask_t1
        BOLD series mask in T1w space
    bold_std
        BOLD series, resampled to template space
    bold_mask_std
        BOLD series mask in template space
    confounds
        TSV of confounds
    surfaces
        BOLD series, resampled to FreeSurfer surfaces
    aroma_noise_ics
        Noise components identified by ICA-AROMA
    melodic_mix
        FSL MELODIC mixing matrix
    bold_cifti
        BOLD CIFTI image
    cifti_variant
        combination of target spaces for `bold_cifti`

    See Also
    --------

    * :py:func:`~niworkflows.func.util.init_bold_reference_wf`
    * :py:func:`~fmriprep.workflows.bold.stc.init_bold_stc_wf`
    * :py:func:`~fmriprep.workflows.bold.hmc.init_bold_hmc_wf`
    * :py:func:`~fmriprep.workflows.bold.t2s.init_bold_t2s_wf`
    * :py:func:`~fmriprep.workflows.bold.registration.init_bold_t1_trans_wf`
    * :py:func:`~fmriprep.workflows.bold.registration.init_bold_reg_wf`
    * :py:func:`~fmriprep.workflows.bold.confounds.init_bold_confs_wf`
    * :py:func:`~fmriprep.workflows.bold.confounds.init_ica_aroma_wf`
    * :py:func:`~fmriprep.workflows.bold.resampling.init_bold_std_trans_wf`
    * :py:func:`~fmriprep.workflows.bold.resampling.init_bold_preproc_trans_wf`
    * :py:func:`~fmriprep.workflows.bold.resampling.init_bold_surf_wf`
    * :py:func:`~sdcflows.workflows.fmap.init_fmap_wf`
    * :py:func:`~sdcflows.workflows.pepolar.init_pepolar_unwarp_wf`
    * :py:func:`~sdcflows.workflows.phdiff.init_phdiff_wf`
    * :py:func:`~sdcflows.workflows.syn.init_syn_sdc_wf`
    * :py:func:`~sdcflows.workflows.unwarp.init_sdc_unwarp_wf`

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.utility import DictMerge
    from niworkflows.workflows.epi.refmap import init_epi_reference_wf

    mem_gb = {"filesize": 1, "resampled": 1, "largemem": 1}
    bold_tlen = 10

    # Have some options handy
    omp_nthreads = config.nipype.omp_nthreads
    freesurfer = config.workflow.run_reconall
    spaces = config.workflow.spaces
    nibabies_dir = str(config.execution.nibabies_dir)

    # Extract BIDS entities and metadata from BOLD file(s)
    entities = extract_entities(bold_file)
    layout = config.execution.layout

    # Take first file as reference
    ref_file = pop_file(bold_file)
    metadata = layout.get_metadata(ref_file)
    # get original image orientation
    ref_orientation = get_img_orientation(ref_file)

    echo_idxs = listify(entities.get("echo", []))
    multiecho = len(echo_idxs) > 2
    if len(echo_idxs) == 1:
        config.loggers.workflow.warning(
            f"Running a single echo <{ref_file}> from a seemingly multi-echo dataset."
        )
        bold_file = ref_file  # Just in case - drop the list

    if len(echo_idxs) == 2:
        raise RuntimeError(
            "Multi-echo processing requires at least three different echos (found two)."
        )

    if multiecho:
        # Drop echo entity for future queries, have a boolean shorthand
        entities.pop("echo", None)
        # reorder echoes from shortest to largest
        tes, bold_file = zip(
            *sorted([(layout.get_metadata(bf)["EchoTime"], bf) for bf in bold_file])
        )
        ref_file = bold_file[0]  # Reset reference to be the shortest TE

    if os.path.isfile(ref_file):
        bold_tlen, mem_gb = _create_mem_gb(ref_file)

    wf_name = _get_wf_name(ref_file)
    config.loggers.workflow.debug(
        f'Creating bold processing workflow for <{ref_file}> ({mem_gb["filesize"]:.2f} GB '
        f'/ {bold_tlen} TRs). Memory resampled/largemem={mem_gb["resampled"]:.2f}'
        f'/{mem_gb["largemem"]:.2f} GB.'
    )

    # Find associated sbref, if possible
    entities["suffix"] = "sbref"
    entities["extension"] = [".nii", ".nii.gz"]  # Overwrite extensions
    sbref_files = layout.get(return_type="file", **entities)

    sbref_msg = f"No single-band-reference found for {os.path.basename(ref_file)}."
    if sbref_files and "sbref" in config.workflow.ignore:
        sbref_msg = "Single-band reference file(s) found and ignored."
    elif sbref_files:
        sbref_msg = "Using single-band reference file(s) {}.".format(
            ",".join([os.path.basename(sbf) for sbf in sbref_files])
        )
    config.loggers.workflow.info(sbref_msg)

    if has_fieldmap:
        # First check if specified via B0FieldSource
        estimator_key = listify(metadata.get("B0FieldSource"))

        if not estimator_key:
            from pathlib import Path
            import re
            from sdcflows.fieldmaps import get_identifier

            # Fallback to IntendedFor
            bold_rel = re.sub(
                r"^sub-[a-zA-Z0-9]*/", "", str(Path(bold_file).relative_to(layout.root))
            )
            estimator_key = get_identifier(bold_rel)

        if not estimator_key:
            has_fieldmap = False
            config.loggers.workflow.critical(
                f"None of the available B0 fieldmaps are associated to <{bold_rel}>"
            )
        else:
            config.loggers.workflow.info(f"Found usable B0 fieldmap <{estimator_key}>")

    # Short circuits: (True and True and (False or 'TooShort')) == 'TooShort'
    run_stc = (
        bool(metadata.get("SliceTiming"))
        and "slicetiming" not in config.workflow.ignore
        and (_get_series_len(ref_file) > 4 or "TooShort")
    )

    # Build workflow
    workflow = Workflow(name=wf_name)
    workflow.__postdesc__ = """\
All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. head-motion
transform matrices, susceptibility distortion correction when available,
and co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `antsApplyTransforms` (ANTs),
configured with Lanczos interpolation to minimize the smoothing
effects of other kernels [@lanczos].
Non-gridded (surface) resamplings were performed using `mri_vol2surf`
(FreeSurfer).
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_file",
                # from smriprep
                "anat_preproc",
                "anat_brain",
                "anat_mask",
                "anat_dseg",
                "anat_tpms",
                "anat_aseg",
                "anat_aparc",
                "anat2std_xfm",
                "std2anat_xfm",
                "template",
                # from bold reference workflow
                "bold_ref",
                "bold_ref_xfm",
                "n_dummy_scans",
                # from sdcflows (optional)
                "fmap",
                "fmap_ref",
                "fmap_coeff",
                "fmap_mask",
                "fmap_id",
                "sdc_method",
                # if reconstructing with FreeSurfer (optional)
                "anat2fsnative_xfm",
                "fsnative2anat_xfm",
                "subject_id",
                "subjects_dir",
            ]
        ),
        name="inputnode",
    )
    inputnode.inputs.bold_file = bold_file

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "bold_anat",
                "bold_anat_ref",
                "bold2anat_xfm",
                "anat2bold_xfm",
                "bold_mask_anat",
                "bold_aseg_anat",
                "bold_aparc_anat",
                "bold_std",
                "bold_std_ref",
                "bold_mask_std",
                "bold_aseg_std",
                "bold_aparc_std",
                "bold_native",
                "bold_cifti",
                "cifti_variant",
                "cifti_metadata",
                "cifti_density",
                "surfaces",
                "confounds",
                "aroma_noise_ics",
                "melodic_mix",
                "nonaggr_denoised_file",
                "confounds_metadata",
            ]
        ),
        name="outputnode",
    )

    # BOLD buffer: an identity used as a pointer to either the original BOLD
    # or the STC'ed one for further use.
    boldbuffer = pe.Node(niu.IdentityInterface(fields=["bold_file"]), name="boldbuffer")

    summary = pe.Node(
        FunctionalSummary(
            slice_timing=run_stc,
            registration=("FSL", "FreeSurfer")[freesurfer],
            registration_dof=config.workflow.bold2t1w_dof,
            registration_init=config.workflow.bold2t1w_init,
            pe_direction=metadata.get("PhaseEncodingDirection"),
            echo_idx=echo_idxs,
            tr=metadata.get("RepetitionTime"),
            orientation=ref_orientation,
        ),
        name="summary",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )
    summary.inputs.dummy_scans = config.workflow.dummy_scans
    # TODO: SDC: make dynamic
    summary.inputs.distortion_correction = "None" if not has_fieldmap else "TOPUP"

    func_derivatives_wf = init_func_derivatives_wf(
        bids_root=layout.root,
        cifti_output=config.workflow.cifti_output,
        freesurfer=freesurfer,
        metadata=metadata,
        output_dir=nibabies_dir,
        spaces=spaces,
        use_aroma=config.workflow.use_aroma,
        debug=config.execution.debug,
    )

    # fmt: off
    workflow.connect([
        (outputnode, func_derivatives_wf, [
            ('bold_anat', 'inputnode.bold_t1'),
            ('bold_anat_ref', 'inputnode.bold_t1_ref'),
            ('bold2anat_xfm', 'inputnode.bold2anat_xfm'),
            ('anat2bold_xfm', 'inputnode.anat2bold_xfm'),
            ('bold_aseg_anat', 'inputnode.bold_aseg_t1'),
            ('bold_aparc_anat', 'inputnode.bold_aparc_t1'),
            ('bold_mask_anat', 'inputnode.bold_mask_t1'),
            ('bold_native', 'inputnode.bold_native'),
            ('confounds', 'inputnode.confounds'),
            ('surfaces', 'inputnode.surf_files'),
            ('aroma_noise_ics', 'inputnode.aroma_noise_ics'),
            ('melodic_mix', 'inputnode.melodic_mix'),
            ('nonaggr_denoised_file', 'inputnode.nonaggr_denoised_file'),
            ('bold_cifti', 'inputnode.bold_cifti'),
            ('cifti_variant', 'inputnode.cifti_variant'),
            ('cifti_metadata', 'inputnode.cifti_metadata'),
            ('cifti_density', 'inputnode.cifti_density'),
            ('confounds_metadata', 'inputnode.confounds_metadata'),
            ('acompcor_masks', 'inputnode.acompcor_masks'),
            ('tcompcor_mask', 'inputnode.tcompcor_mask'),
        ]),
    ])
    # fmt: on

    # Extract BOLD validation from init_bold_reference_wf
    val_bold = pe.MapNode(
        ValidateImage(),
        name="val_bold",
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
        iterfield=["in_file"],
    )
    val_bold.inputs.in_file = listify(bold_file)

    # Top-level BOLD splitter
    bold_split = pe.Node(FSLSplit(dimension="t"), name="bold_split", mem_gb=mem_gb["filesize"] * 3)

    # HMC on the BOLD
    bold_hmc_wf = init_bold_hmc_wf(
        name="bold_hmc_wf", mem_gb=mem_gb["filesize"], omp_nthreads=omp_nthreads
    )

    # calculate BOLD registration to T1w
    bold_reg_wf = init_bold_reg_wf(
        bold2t1w_dof=config.workflow.bold2t1w_dof,
        bold2t1w_init=config.workflow.bold2t1w_init,
        freesurfer=freesurfer,
        mem_gb=mem_gb["resampled"],
        name="bold_reg_wf",
        omp_nthreads=omp_nthreads,
        sloppy=config.execution.sloppy,
        use_bbr=config.workflow.use_bbr,
    )

    # apply BOLD registration to T1w
    bold_t1_trans_wf = init_bold_t1_trans_wf(
        name="bold_t1_trans_wf",
        freesurfer=freesurfer,
        mem_gb=mem_gb["resampled"],
        omp_nthreads=omp_nthreads,
        use_compression=False,
    )
    bold_t1_trans_wf.inputs.inputnode.fieldwarp = "identity"

    # get confounds
    bold_confounds_wf = init_bold_confs_wf(
        mem_gb=mem_gb["largemem"],
        metadata=metadata,
        freesurfer=freesurfer,
        regressors_all_comps=config.workflow.regressors_all_comps,
        regressors_fd_th=config.workflow.regressors_fd_th,
        regressors_dvars_th=config.workflow.regressors_dvars_th,
        fd_radius=config.workflow.fd_radius,
        name="bold_confounds_wf",
    )
    bold_confounds_wf.get_node("inputnode").inputs.t1_transform_flags = [False]

    # Apply transforms in 1 shot
    # Only use uncompressed output if AROMA is to be run
    bold_bold_trans_wf = init_bold_preproc_trans_wf(
        mem_gb=mem_gb["resampled"],
        omp_nthreads=omp_nthreads,
        use_compression=not config.execution.low_mem,
        use_fieldwarp=has_fieldmap,
        name="bold_bold_trans_wf",
    )
    bold_bold_trans_wf.inputs.inputnode.name_source = ref_file
    bold_bold_trans_wf.inputs.inputnode.fieldwarp = "identity"

    # SLICE-TIME CORRECTION (or bypass) #############################################
    if run_stc is True:  # bool('TooShort') == True, so check True explicitly
        bold_stc_wf = init_bold_stc_wf(name="bold_stc_wf", metadata=metadata)
        # fmt: off
        workflow.connect([
            (inputnode, bold_stc_wf, [('n_dummy_scans', 'inputnode.skip_vols')]),
            (bold_stc_wf, boldbuffer, [('outputnode.stc_file', 'bold_file')]),
        ])
        # fmt: on
        if not multiecho:
            # fmt: off
            workflow.connect([
                (val_bold, bold_stc_wf, [(("out_file", pop_file), 'inputnode.bold_file')])])
            # fmt: on
        else:  # for meepi, iterate through stc_wf for all workflows
            meepi_echos = boldbuffer.clone(name="meepi_echos")
            meepi_echos.iterables = ("bold_file", bold_file)
            workflow.connect([(meepi_echos, bold_stc_wf, [("bold_file", "inputnode.bold_file")])])
    elif not multiecho:  # STC is too short or False
        # bypass STC from original BOLD to the splitter through boldbuffer
        # fmt: off
        workflow.connect([
            (val_bold, boldbuffer, [(("out_file", pop_file), 'bold_file')])])
        # fmt: on
    else:
        # for meepi, iterate over all meepi echos to boldbuffer
        boldbuffer.iterables = ("bold_file", bold_file)

    # MULTI-ECHO EPI DATA #############################################
    if multiecho:  # instantiate relevant interfaces, imports
        from niworkflows.func.util import init_skullstrip_bold_wf

        skullstrip_bold_wf = init_skullstrip_bold_wf(name="skullstrip_bold_wf")

        split_opt_comb = bold_split.clone(name="split_opt_comb")

        inputnode.inputs.bold_file = ref_file  # Replace reference w first echo

        join_echos = pe.JoinNode(
            niu.IdentityInterface(fields=["bold_files", "skullstripped_bold_files"]),
            joinsource=("meepi_echos" if run_stc is True else "boldbuffer"),
            joinfield=["bold_files", "skullstripped_bold_files"],
            name="join_echos",
        )

        # create optimal combination, adaptive T2* map
        bold_t2s_wf = init_bold_t2s_wf(
            echo_times=tes,
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            name="bold_t2smap_wf",
        )

    # Mask input BOLD reference image
    initial_boldref_mask = pe.Node(BrainExtraction(), name="initial_boldref_mask")

    # This final boldref will be calculated after bold_bold_trans_wf, which includes one or more:
    # HMC (head motion correction)
    # STC (slice time correction)
    # SDC (susceptibility distortion correction)
    final_boldref_wf = init_epi_reference_wf(
        auto_bold_nss=True,
        omp_nthreads=omp_nthreads,
        name="final_boldref_wf",
    )
    final_boldref_mask = pe.Node(BrainExtraction(), name="final_boldref_mask")

    # MAIN WORKFLOW STRUCTURE #######################################################
    # fmt: off
    workflow.connect([
        # BOLD buffer has slice-time corrected if it was run, original otherwise
        (boldbuffer, bold_split, [('bold_file', 'in_file')]),
        # HMC
        (inputnode, bold_hmc_wf, [
            ('bold_ref', 'inputnode.raw_ref_image')]),
        (inputnode, initial_boldref_mask, [('bold_ref', 'in_file')]),
        (val_bold, bold_hmc_wf, [
            (("out_file", pop_file), 'inputnode.bold_file')]),
        (inputnode, summary, [
            ('n_dummy_scans', 'algo_dummy_scans')]),
        # EPI-T1 registration workflow
        (inputnode, bold_reg_wf, [
            ('anat_dseg', 'inputnode.t1w_dseg'),
            # Undefined if --fs-no-reconall, but this is safe
            ('subjects_dir', 'inputnode.subjects_dir'),
            ('subject_id', 'inputnode.subject_id'),
            ('fsnative2anat_xfm', 'inputnode.fsnative2t1w_xfm')]),
        (inputnode, bold_reg_wf, [
            ('anat_brain', 'inputnode.t1w_brain')]),
        (inputnode, bold_t1_trans_wf, [
            ('bold_file', 'inputnode.name_source'),
            ('anat_mask', 'inputnode.t1w_mask'),
            ('anat_brain', 'inputnode.t1w_brain'),
            ('anat_aseg', 'inputnode.t1w_aseg'),
            ('anat_aparc', 'inputnode.t1w_aparc')]),
        (bold_reg_wf, outputnode, [
            ('outputnode.itk_bold_to_t1', 'bold2anat_xfm'),
            ('outputnode.itk_t1_to_bold', 'anat2bold_xfm')]),
        (bold_reg_wf, bold_t1_trans_wf, [
            ('outputnode.itk_bold_to_t1', 'inputnode.itk_bold_to_t1')]),
        (bold_t1_trans_wf, outputnode, [('outputnode.bold_t1', 'bold_anat'),
                                        ('outputnode.bold_t1_ref', 'bold_anat_ref'),
                                        ('outputnode.bold_aseg_t1', 'bold_aseg_anat'),
                                        ('outputnode.bold_aparc_t1', 'bold_aparc_anat')]),
        (bold_reg_wf, summary, [('outputnode.fallback', 'fallback')]),
        # Connect bold_confounds_wf
        (inputnode, bold_confounds_wf, [('anat_tpms', 'inputnode.t1w_tpms'),
                                        ('anat_mask', 'inputnode.t1w_mask')]),
        (bold_hmc_wf, bold_confounds_wf, [
            ('outputnode.movpar_file', 'inputnode.movpar_file'),
            ('outputnode.rmsd_file', 'inputnode.rmsd_file')]),
        (bold_reg_wf, bold_confounds_wf, [
            ('outputnode.itk_t1_to_bold', 'inputnode.t1_bold_xform')]),
        (inputnode, bold_confounds_wf, [
            ('n_dummy_scans', 'inputnode.skip_vols')]),
        (bold_confounds_wf, outputnode, [
            ('outputnode.confounds_file', 'confounds'),
            ('outputnode.confounds_metadata', 'confounds_metadata'),
            ('outputnode.acompcor_masks', 'acompcor_masks'),
            ('outputnode.tcompcor_mask', 'tcompcor_mask'),
        ]),
        # Connect bold_bold_trans_wf
        (bold_split, bold_bold_trans_wf, [
            ('out_files', 'inputnode.bold_file')]),
        (bold_hmc_wf, bold_bold_trans_wf, [
            ('outputnode.xforms', 'inputnode.hmc_xforms')]),
        # Summary
        (outputnode, summary, [('confounds', 'confounds_file')]),
        (final_boldref_wf, final_boldref_mask, [('outputnode.epi_ref_file', 'in_file')]),
    ])
    # fmt: on

    # for standard EPI data, pass along correct file
    if not multiecho:
        # TODO: Add SDC
        # fmt: off
        workflow.connect([
            (inputnode, func_derivatives_wf, [
                ('bold_file', 'inputnode.source_file')]),
            (bold_bold_trans_wf, final_boldref_wf, [
                ('outputnode.bold', 'inputnode.in_files')]),
            (bold_bold_trans_wf, bold_confounds_wf, [
                ('outputnode.bold', 'inputnode.bold')]),
            (bold_split, bold_t1_trans_wf, [
                ('out_files', 'inputnode.bold_split')]),
            (bold_hmc_wf, bold_t1_trans_wf, [
                ('outputnode.xforms', 'inputnode.hmc_xforms')]),
        ])
        # fmt: on
    else:  # for meepi, use optimal combination
        # fmt: off
        workflow.connect([
            # update name source for optimal combination
            (inputnode, func_derivatives_wf, [
                (('bold_file', combine_meepi_source), 'inputnode.source_file')]),
            (bold_bold_trans_wf, join_echos, [
                ('outputnode.bold', 'bold_files')]),
            (join_echos, final_boldref_wf, [
                ('bold_files', 'inputnode.in_files')]),
            # TODO: Check with multi-echo data
            (bold_bold_trans_wf, skullstrip_bold_wf, [
                ('outputnode.bold', 'inputnode.in_file')]),
            (skullstrip_bold_wf, join_echos, [
                ('outputnode.skull_stripped_file', 'skullstripped_bold_files')]),
            (join_echos, bold_t2s_wf, [
                ('skullstripped_bold_files', 'inputnode.bold_file')]),
            (bold_t2s_wf, bold_confounds_wf, [
                ('outputnode.bold', 'inputnode.bold')]),
            (bold_t2s_wf, split_opt_comb, [
                ('outputnode.bold', 'in_file')]),
            (split_opt_comb, bold_t1_trans_wf, [
                ('out_files', 'inputnode.bold_split')]),
        ])
        # fmt: on

        # Already applied in bold_bold_trans_wf, which inputs to bold_t2s_wf
        bold_t1_trans_wf.inputs.inputnode.hmc_xforms = "identity"

    # Map final BOLD mask into T1w space (if required)
    nonstd_spaces = set(spaces.get_nonstandard())
    if nonstd_spaces.intersection(("T1w", "anat")):
        from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms

        boldmask_to_t1w = pe.Node(
            ApplyTransforms(interpolation="MultiLabel"), name="boldmask_to_t1w", mem_gb=0.1
        )
        # fmt: off
        workflow.connect([
            (bold_reg_wf, boldmask_to_t1w, [
                ('outputnode.itk_bold_to_t1', 'transforms')]),
            (bold_t1_trans_wf, boldmask_to_t1w, [
                ('outputnode.bold_mask_t1', 'reference_image')]),
            (boldmask_to_t1w, outputnode, [
                ('output_image', 'bold_mask_anat')]),
        ])
        # fmt: on

    if nonstd_spaces.intersection(("func", "run", "bold", "boldref", "sbref")):
        workflow.connect(
            [
                (
                    bold_bold_trans_wf if not multiecho else bold_t2s_wf,
                    outputnode,
                    [("outputnode.bold", "bold_native")],
                )
            ]
        )

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # Apply transforms in 1 shot
        # Only use uncompressed output if AROMA is to be run
        bold_std_trans_wf = init_bold_std_trans_wf(
            freesurfer=freesurfer,
            mem_gb=mem_gb["resampled"],
            omp_nthreads=omp_nthreads,
            spaces=spaces,
            name="bold_std_trans_wf",
            use_compression=not config.execution.low_mem,
        )
        bold_std_trans_wf.inputs.inputnode.fieldwarp = "identity"

        # fmt: off
        workflow.connect([
            (inputnode, bold_std_trans_wf, [
                ('template', 'inputnode.templates'),
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('bold_file', 'inputnode.name_source'),
                ('anat_aseg', 'inputnode.bold_aseg'),
                ('anat_aparc', 'inputnode.bold_aparc')]),
            (bold_reg_wf, bold_std_trans_wf, [
                ('outputnode.itk_bold_to_t1', 'inputnode.itk_bold_to_t1')]),
            (bold_std_trans_wf, outputnode, [('outputnode.bold_std', 'bold_std'),
                                             ('outputnode.bold_std_ref', 'bold_std_ref'),
                                             ('outputnode.bold_mask_std', 'bold_mask_std')]),
        ])
        # fmt: on

        if freesurfer:
            # fmt: off
            workflow.connect([
                (bold_std_trans_wf, func_derivatives_wf, [
                    ('outputnode.bold_aseg_std', 'inputnode.bold_aseg_std'),
                    ('outputnode.bold_aparc_std', 'inputnode.bold_aparc_std'),
                ]),
                (bold_std_trans_wf, outputnode, [
                    ('outputnode.bold_aseg_std', 'bold_aseg_std'),
                    ('outputnode.bold_aparc_std', 'bold_aparc_std')]),
            ])
            # fmt: on

        if not multiecho:
            # fmt: off
            workflow.connect([
                (bold_split, bold_std_trans_wf, [
                    ('out_files', 'inputnode.bold_split')]),
                (bold_hmc_wf, bold_std_trans_wf, [
                    ('outputnode.xforms', 'inputnode.hmc_xforms')]),
            ])
            # fmt: on
        else:
            # fmt: off
            workflow.connect([
                (split_opt_comb, bold_std_trans_wf, [
                    ('out_files', 'inputnode.bold_split')])
            ])
            # fmt: on

            # Already applied in bold_bold_trans_wf, which inputs to bold_t2s_wf
            bold_std_trans_wf.inputs.inputnode.hmc_xforms = "identity"

        # func_derivatives_wf internally parametrizes over snapshotted spaces.
        # fmt: off
        workflow.connect([
            (bold_std_trans_wf, func_derivatives_wf, [
                ('outputnode.template', 'inputnode.template'),
                ('outputnode.spatial_reference', 'inputnode.spatial_reference'),
                ('outputnode.bold_std_ref', 'inputnode.bold_std_ref'),
                ('outputnode.bold_std', 'inputnode.bold_std'),
                ('outputnode.bold_mask_std', 'inputnode.bold_mask_std'),
            ]),
        ])
        # fmt: on

        if config.workflow.use_aroma:  # ICA-AROMA workflow
            from .confounds import init_ica_aroma_wf

            ica_aroma_wf = init_ica_aroma_wf(
                mem_gb=mem_gb["resampled"],
                metadata=metadata,
                omp_nthreads=omp_nthreads,
                err_on_aroma_warn=config.workflow.aroma_err_on_warn,
                aroma_melodic_dim=config.workflow.aroma_melodic_dim,
                name="ica_aroma_wf",
            )

            join = pe.Node(
                niu.Function(output_names=["out_file"], function=_to_join), name="aroma_confounds"
            )

            mrg_conf_metadata = pe.Node(
                niu.Merge(2), name="merge_confound_metadata", run_without_submitting=True
            )
            mrg_conf_metadata2 = pe.Node(
                DictMerge(), name="merge_confound_metadata2", run_without_submitting=True
            )
            # fmt: off
            workflow.disconnect([
                (bold_confounds_wf, outputnode, [
                    ('outputnode.confounds_file', 'confounds'),
                ]),
                (bold_confounds_wf, outputnode, [
                    ('outputnode.confounds_metadata', 'confounds_metadata'),
                ]),
            ])
            workflow.connect([
                (inputnode, ica_aroma_wf, [
                    ('bold_file', 'inputnode.name_source')]),
                (bold_hmc_wf, ica_aroma_wf, [
                    ('outputnode.movpar_file', 'inputnode.movpar_file')]),
                (inputnode, ica_aroma_wf, [
                    ('n_dummy_scans', 'inputnode.skip_vols')]),
                (bold_confounds_wf, join, [
                    ('outputnode.confounds_file', 'in_file')]),
                (bold_confounds_wf, mrg_conf_metadata,
                    [('outputnode.confounds_metadata', 'in1')]),
                (ica_aroma_wf, join,
                    [('outputnode.aroma_confounds', 'join_file')]),
                (ica_aroma_wf, mrg_conf_metadata,
                    [('outputnode.aroma_metadata', 'in2')]),
                (mrg_conf_metadata, mrg_conf_metadata2, [('out', 'in_dicts')]),
                (ica_aroma_wf, outputnode,
                    [('outputnode.aroma_noise_ics', 'aroma_noise_ics'),
                     ('outputnode.melodic_mix', 'melodic_mix'),
                     ('outputnode.nonaggr_denoised_file', 'nonaggr_denoised_file')]),
                (join, outputnode, [('out_file', 'confounds')]),
                (mrg_conf_metadata2, outputnode, [('out_dict', 'confounds_metadata')]),
                (bold_std_trans_wf, ica_aroma_wf, [
                    ('outputnode.bold_std', 'inputnode.bold_std'),
                    ('outputnode.bold_mask_std', 'inputnode.bold_mask_std'),
                    ('outputnode.spatial_reference', 'inputnode.spatial_reference')]),
            ])
            # fmt: on

    # SURFACES ##################################################################################
    # Freesurfer
    freesurfer_spaces = spaces.get_fs_spaces()
    if freesurfer and freesurfer_spaces:
        config.loggers.workflow.debug("Creating BOLD surface-sampling workflow.")
        bold_surf_wf = init_bold_surf_wf(
            mem_gb=mem_gb["resampled"],
            surface_spaces=freesurfer_spaces,
            medial_surface_nan=config.workflow.medial_surface_nan,
            name="bold_surf_wf",
        )

        # fmt: off
        workflow.connect([
            (inputnode, bold_surf_wf, [
                ('subjects_dir', 'inputnode.subjects_dir'),
                ('subject_id', 'inputnode.subject_id'),
                ('anat2fsnative_xfm', 'inputnode.t1w2fsnative_xfm')]),
            (bold_t1_trans_wf, bold_surf_wf, [('outputnode.bold_t1', 'inputnode.source_file')]),
            (bold_surf_wf, outputnode, [('outputnode.surfaces', 'surfaces')]),
            (bold_surf_wf, func_derivatives_wf, [
                ('outputnode.target', 'inputnode.surf_refs')]),
        ])
        # fmt: on

        # CIFTI output
        if config.workflow.cifti_output:
            from .alignment import init_subcortical_rois_wf, init_subcortical_mni_alignment_wf
            from .resampling import init_bold_grayords_wf

            key = None
            for space in spaces.get_spaces(nonstandard=False, dim=(3,)):
                if "MNIInfant" in space:
                    key = space.replace(":", "_")

            assert key is not None, f"MNIInfant not found in SpatialReferences: {spaces}"

            # BOLD/ROIs should be in MNIInfant space
            cifti_select_std = pe.Node(
                KeySelect(fields=["bold_std", "bold_aseg_std"], key=key),
                name="cifti_select_std",
                run_without_submitting=True,
            )

            subcortical_rois_wf = init_subcortical_rois_wf()
            subcortical_mni_alignment_wf = init_subcortical_mni_alignment_wf()
            bold_grayords_wf = init_bold_grayords_wf(
                grayord_density=config.workflow.cifti_output,
                mem_gb=mem_gb["resampled"],
                repetition_time=metadata["RepetitionTime"],
            )

            # fmt: off
            workflow.connect([
                (bold_std_trans_wf, cifti_select_std, [
                    ("outputnode.bold_std", "bold_std"),
                    ("outputnode.bold_aseg_std", "bold_aseg_std"),
                    ("outputnode.spatial_reference", "keys")]),
                (cifti_select_std, subcortical_rois_wf, [
                    ("bold_aseg_std", "inputnode.MNIInfant_aseg")]),
                (cifti_select_std, subcortical_mni_alignment_wf, [
                    ("bold_std", "inputnode.MNIInfant_bold")]),
                (subcortical_rois_wf, subcortical_mni_alignment_wf, [
                    ("outputnode.MNIInfant_rois", "inputnode.MNIInfant_rois"),
                    ("outputnode.MNI152_rois", "inputnode.MNI152_rois")]),
                (subcortical_mni_alignment_wf, bold_grayords_wf, [
                    ("outputnode.subcortical_volume", "inputnode.subcortical_volume"),
                    ("outputnode.subcortical_labels", "inputnode.subcortical_labels")]),
                (bold_surf_wf, bold_grayords_wf, [
                    ('outputnode.surfaces', 'inputnode.surf_files'),
                    ('outputnode.target', 'inputnode.surf_refs')]),
                (bold_grayords_wf, outputnode, [
                    ('outputnode.cifti_bold', 'bold_cifti'),
                    ('outputnode.cifti_variant', 'cifti_variant'),
                    ('outputnode.cifti_metadata', 'cifti_metadata'),
                    ('outputnode.cifti_density', 'cifti_density')]),
            ])
            # fmt: on

    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        if not config.workflow.cifti_output:
            config.loggers.workflow.critical("The carpetplot requires CIFTI outputs")
        else:
            carpetplot_wf = init_carpetplot_wf(
                mem_gb=mem_gb["resampled"],
                metadata=metadata,
                cifti_output=bool(config.workflow.cifti_output),
                name="carpetplot_wf",
            )

            # fmt: off
            workflow.connect([
                (bold_grayords_wf, carpetplot_wf, [
                    ('outputnode.cifti_bold', 'inputnode.cifti_bold')]),
                (bold_confounds_wf, carpetplot_wf, [
                    ('outputnode.confounds_file', 'inputnode.confounds_file')]),
            ])
            # fmt: on

    # REPORTING ############################################################
    ds_report_summary = pe.Node(
        DerivativesDataSink(desc="summary", datatype="figures", dismiss_entities=("echo",)),
        name="ds_report_summary",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    ds_report_validation = pe.Node(
        DerivativesDataSink(desc="validation", datatype="figures", dismiss_entities=("echo",)),
        name="ds_report_validation",
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # fmt: off
    workflow.connect([
        (summary, ds_report_summary, [('out_report', 'in_file')]),
        (val_bold, ds_report_validation, [
            (("out_report", pop_file), 'in_file')]),
    ])
    # fmt: on

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split(".")[-1].startswith("ds_report"):
            workflow.get_node(node).inputs.base_directory = nibabies_dir
            workflow.get_node(node).inputs.source_file = ref_file

    # fmt: off
    workflow.connect([
        (final_boldref_mask, bold_t1_trans_wf, [
            ('out_mask', 'inputnode.ref_bold_mask'),
            ('out_file', 'inputnode.ref_bold_brain'),
        ]),
        (final_boldref_mask, bold_reg_wf, [
            ('out_file', 'inputnode.ref_bold_brain'),
        ]),
        (final_boldref_mask, bold_confounds_wf, [
            ('out_mask', 'inputnode.bold_mask')
        ]),
    ])
    # fmt: on
    if nonstd_spaces.intersection(("T1w", "anat")):
        # fmt: off
        workflow.connect([
            (final_boldref_mask, boldmask_to_t1w, [
                ('out_mask', 'input_image')]),
        ])
        # fmt: on
    if nonstd_spaces.intersection(("func", "run", "bold", "boldref", "sbref")):
        # fmt: off
        workflow.connect([
            (final_boldref_mask, func_derivatives_wf, [
                ('out_file', 'inputnode.bold_native_ref'),
                ('out_mask', 'inputnode.bold_mask_native')]),
        ])
        # fmt: on
    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        # fmt: off
        workflow.connect([
            (final_boldref_mask, bold_std_trans_wf, [
                ('out_mask', 'inputnode.bold_mask')]),
        ])
        # fmt: on

    if not has_fieldmap:
        summary.inputs.distortion_correction = "None"
        return workflow

    # SDC
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from sdcflows.workflows.apply.registration import init_coeff2epi_wf
    from sdcflows.workflows.apply.correction import init_unwarp_wf

    coeff2epi_wf = init_coeff2epi_wf(
        debug="fieldmaps" in config.execution.debug,
        omp_nthreads=config.nipype.omp_nthreads,
        write_coeff=True,
    )
    unwarp_wf = init_unwarp_wf(
        debug="fieldmaps" in config.execution.debug, omp_nthreads=config.nipype.omp_nthreads
    )
    unwarp_wf.inputs.inputnode.metadata = layout.get_metadata(str(bold_file))

    output_select = pe.Node(
        KeySelect(fields=["fmap", "fmap_ref", "fmap_coeff", "fmap_mask", "sdc_method"]),
        name="output_select",
        run_without_submitting=True,
    )
    output_select.inputs.key = estimator_key[0]
    if len(estimator_key) > 1:
        config.loggers.workflow.warning(
            f"Several fieldmaps <{', '.join(estimator_key)}> are "
            f"'IntendedFor' <{bold_file}>, using {estimator_key[0]}"
        )

    sdc_report = pe.Node(
        SimpleBeforeAfter(before_label="Distorted", after_label="Corrected", dismiss_affine=True),
        name="sdc_report",
        mem_gb=0.1,
    )

    ds_report_sdc = pe.Node(
        DerivativesDataSink(
            base_directory=nibabies_dir,
            desc="sdc",
            suffix="bold",
            datatype="figures",
            dismiss_entities=("echo",),
        ),
        name="ds_report_sdc",
        run_without_submitting=True,
    )

    # fmt: off
    workflow.connect([
        (inputnode, output_select, [("fmap", "fmap"),
                                    ("fmap_ref", "fmap_ref"),
                                    ("fmap_coeff", "fmap_coeff"),
                                    ("fmap_mask", "fmap_mask"),
                                    ("sdc_method", "sdc_method"),
                                    ("fmap_id", "keys")]),
        (output_select, coeff2epi_wf, [
            ("fmap_ref", "inputnode.fmap_ref"),
            ("fmap_coeff", "inputnode.fmap_coeff"),
            ("fmap_mask", "inputnode.fmap_mask")]),
        (output_select, summary, [("sdc_method", "distortion_correction")]),
        (inputnode, coeff2epi_wf, [
            ("bold_ref", "inputnode.target_ref")]),
        (initial_boldref_mask, coeff2epi_wf, [
            ("out_file", "inputnode.target_mask")]),  # skull-stripped brain
        (inputnode, unwarp_wf, [("bold_ref", "inputnode.distorted")]),
        (coeff2epi_wf, unwarp_wf, [
            ("outputnode.fmap_coeff", "inputnode.fmap_coeff")]),
        (inputnode, sdc_report, [("bold_ref", "before")]),
        (unwarp_wf, sdc_report, [("outputnode.corrected", "after"),
                                 ("outputnode.corrected_mask", "wm_seg")]),
        (inputnode, ds_report_sdc, [("bold_file", "source_file")]),
        (sdc_report, ds_report_sdc, [("out_report", "in_file")]),
        (unwarp_wf, bold_bold_trans_wf, [('outputnode.fieldwarp', 'inputnode.fieldwarp')]),
    ])
    # fmt: on

    if not multiecho:
        # fmt: off
        workflow.connect([
            (unwarp_wf, bold_t1_trans_wf, [
                ('outputnode.fieldwarp', 'inputnode.fieldwarp')]),
            (unwarp_wf, bold_std_trans_wf, [
                ('outputnode.fieldwarp', 'inputnode.fieldwarp')]),
        ])
        # fmt: on

    return workflow


def _get_series_len(bold_fname):
    from nipype.algorithms.confounds import is_outlier
    import nibabel as nb
    import numpy as np

    img = nb.load(bold_fname)
    if len(img.shape) < 4:
        return 1

    data = img.get_fdata(dtype="float32")
    # Data can come with outliers showing very high numbers - preemptively prune
    data = np.clip(
        data,
        a_min=0.0,
        a_max=np.percentile(data, 99.8),
    )
    outliers = is_outlier(np.mean(data, axis=(0, 1, 2)))
    return img.shape[3] - outliers


def _create_mem_gb(bold_fname):
    bold_size_gb = os.path.getsize(bold_fname) / (1024 ** 3)
    bold_tlen = nb.load(bold_fname).shape[-1]
    mem_gb = {
        "filesize": bold_size_gb,
        "resampled": bold_size_gb * 4,
        "largemem": bold_size_gb * (max(bold_tlen / 100, 1.0) + 4),
    }

    return bold_tlen, mem_gb


def _get_wf_name(bold_fname):
    """
    Derive the workflow name for supplied BOLD file.
    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_bold.nii.gz')
    'func_preproc_task_nback_wf'
    >>> _get_wf_name('/completely/made/up/path/sub-01_task-nback_run-01_echo-1_bold.nii.gz')
    'func_preproc_task_nback_run_01_echo_1_wf'
    """
    from nipype.utils.filemanip import split_filename

    fname = split_filename(bold_fname)[1]
    fname_nosub = "_".join(fname.split("_")[1:])
    # if 'echo' in fname_nosub:
    #     fname_nosub = '_'.join(fname_nosub.split("_echo-")[:1]) + "_bold"
    name = "func_preproc_" + fname_nosub.replace(".", "_").replace(" ", "").replace(
        "-", "_"
    ).replace("_bold", "_wf")

    return name


def _to_join(in_file, join_file):
    """Join two tsv files if the join_file is not ``None``."""
    from niworkflows.interfaces.utils import JoinTSVColumns

    if join_file is None:
        return in_file
    res = JoinTSVColumns(in_file=in_file, join_file=join_file).run()
    return res.outputs.out_file


def get_img_orientation(imgf):
    """Return the image orientation as a string"""
    img = nb.load(imgf)
    return "".join(nb.aff2axcodes(img.affine))
