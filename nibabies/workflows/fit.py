import logging

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.utils.connections import pop_file

from nibabies import config

LOGGER = logging.getLogger('nipype.workflow')

def init_infant_anat_fit_wf(
    age_months,
    t1w,
    t2w,
    bids_root,
    precomputed,
    hires,
    longitudinal,
    omp_nthreads,
    output_dir,
    segmentation_atlases,
    skull_strip_mode,
    skull_strip_template,
    sloppy,
    spaces,
    cifti_output,
    name='infant_anat_fit_wf',
):
    """
    Stage the anatomical preprocessing steps:
    - T1w reference
    - T2w reference
    - Brain extraction and INU (bias field) correction
    - Brain tissue segmentation
    - Spatial normalization to standard spaces.
    - Surface reconstruction (MCRIBS / infant_recon_all / recon-all)
    """

    workflow = Workflow(name=name)
    num_t1w = len(t1w)
    num_t2w = len(t2w)

    if not num_t1w and not num_t2w:
        raise FileNotFoundError('No anatomical scans provided!')

    if not num_t1w or not num_t2w:
        modality = 'T1w' if num_t1w else 'T2w'
        anatomicals = t1w or t2w

        workflow = init_infant_single_anat_fit_wf(
            modality,
            age_months=age_months,
            anatomicals=anatomicals,
            bids_root=bids_root,
            precomputed=precomputed,
            hires=hires,
            longitudinal=longitudinal,
            omp_nthreads=omp_nthreads,
            output_dir=output_dir,
            segmentation_atlases=segmentation_atlases,
            skull_strip_mode=skull_strip_mode,
            skull_strip_template=skull_strip_mode,
            sloppy=sloppy,
            spaces=spaces,
            cifti_output=cifti_output,
        )

        return workflow

    # Organization
    # ------------
    # This workflow takes the usual (inputnode -> graph -> outputnode) format
    # The graph consists of (input -> compute -> datasink -> buffer) units,
    # and all inputs to outputnode are buffer.
    # If precomputed inputs are found, then these units are replaced with (buffer)
    #     At the time of writing, t1w_mask is an exception, which takes the form
    #     (t1w_buffer -> refined_buffer -> datasink -> outputnode)
    # All outputnode components should therefore point to files in the input or
    # output directories.
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['t1w', 't2w', 'roi', 'flair', 'subjects_dir', 'subject_id']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                # Primary derivatives
                't1w_preproc',
                't2w_preproc',
                't1w2t2w_xfm',
                't1w_mask',
                't1w_dseg',
                't1w_tpms',
                'anat2std_xfm',
                'fsnative2t1w_xfm',
                # Surface and metric derivatives for fsLR resampling
                'white',
                'pial',
                'midthickness',
                'sphere',
                'thickness',
                'sulc',
                'sphere_reg',
                'sphere_reg_fsLR',
                'sphere_reg_msm',
                'anat_ribbon',
                # Reverse transform; not computable from forward transform
                'std2anat_xfm',
                # Metadata
                'template',
                'subjects_dir',
                'subject_id',
                't1w_valid_list',
            ]
        ),
        name='outputnode',
    )

    # If all derivatives exist, inputnode could go unconnected, so add explicitly
    workflow.add_nodes([inputnode])

    # Stage 1 inputs (filtered)
    sourcefile_buffer = pe.Node(
        niu.IdentityInterface(fields=['source_files']),
        name='sourcefile_buffer',
    )

    # Stage 2 - Anatomicals
    t1w_buffer = pe.Node(
        niu.IdentityInterface(fields=['t1w_preproc', 't1w_mask', 't1w_brain']),
        name='t1w_buffer',
    )
    t2w_buffer = pe.Node(
        niu.IdentityInterface(fields=['t2w_preproc', 't2w_mask', 't2w_brain'])
    )

    # Stage 3 - Coregistration
    t1w2t2w_buffer = pe.Node(niu.Merge(2), name='t1w2t2w_buffer')
    t2w2t1w_buffer = pe.Node(niu.Merge(2), name='t2w2t1w_buffer')

    # Stage 4 - Segmentation
    seg_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_dseg', 'anat_tpms']),
        name='seg_buffer',
    )
    # Stage 5 - collated template names, forward and reverse transforms
    template_buffer = pe.Node(niu.Merge(2), name='template_buffer')
    anat2std_buffer = pe.Node(niu.Merge(2), name='anat2std_buffer')
    std2anat_buffer = pe.Node(niu.Merge(2), name='std2anat_buffer')

    # Stage 6 results: Refined stage 2 results; may be direct copy if no refinement
    refined_buffer = pe.Node(
        niu.IdentityInterface(fields=['t1w_mask', 't1w_brain']),
        name='refined_buffer',
    )

    # Stage 8 results: GIFTI surfaces
    surfaces_buffer = pe.Node(
        niu.IdentityInterface(
            fields=['white', 'pial', 'midthickness', 'sphere', 'sphere_reg', 'thickness', 'sulc']
        ),
        name='surfaces_buffer',
    )

    # Stage 9 and 10 results: fsLR sphere registration
    fsLR_buffer = pe.Node(niu.IdentityInterface(fields=['sphere_reg_fsLR']), name='fsLR_buffer')
    msm_buffer = pe.Node(niu.IdentityInterface(fields=['sphere_reg_msm']), name='msm_buffer')

    workflow.connect([
        (seg_buffer, outputnode, [
            ('t1w_dseg', 't1w_dseg'),
            ('t1w_tpms', 't1w_tpms'),
        ]),
        (anat2std_buffer, outputnode, [('out', 'anat2std_xfm')]),
        (std2anat_buffer, outputnode, [('out', 'std2anat_xfm')]),
        (template_buffer, outputnode, [('out', 'template')]),
        (sourcefile_buffer, outputnode, [('source_files', 't1w_valid_list')]),
        (surfaces_buffer, outputnode, [
            ('white', 'white'),
            ('pial', 'pial'),
            ('midthickness', 'midthickness'),
            ('sphere', 'sphere'),
            ('sphere_reg', 'sphere_reg'),
            ('thickness', 'thickness'),
            ('sulc', 'sulc'),
        ]),
        (fsLR_buffer, outputnode, [('sphere_reg_fsLR', 'sphere_reg_fsLR')]),
        (msm_buffer, outputnode, [('sphere_reg_msm', 'sphere_reg_msm')]),
    ])  # fmt:skip

    # Reporting
    recon_method = config.workflow.surface_recon_method
    anat_reports_wf = init_anat_reports_wf(
        surface_recon=recon_method,
        output_dir=output_dir,
        sloppy=sloppy,
    )

    workflow.connect([
        (outputnode, anat_reports_wf, [
            ('t1w_valid_list', 'inputnode.source_file'),
            ('t1w_preproc', 'inputnode.t1w_preproc'),
            ('t1w_mask', 'inputnode.t1w_mask'),
            ('t1w_dseg', 'inputnode.t1w_dseg'),
            ('template', 'inputnode.template'),
            ('anat2std_xfm', 'inputnode.anat2std_xfm'),
            ('subjects_dir', 'inputnode.subjects_dir'),
            ('subject_id', 'inputnode.subject_id'),
        ]),
    ])  # fmt:skip

    desc = (
        '\nAnatomical data preprocessing\n\n: ',
        f'A total of {len(t1w)} T1w and {len(t2w)} T2w images '
        'were found within the input BIDS dataset.'
    )

    # Stage 1: Conform & valid T1w/T2w images
    t1w_validate = pe.Node(ValidateImage(), name='anat_validate', run_without_submitting=True)
    t2w_validate = t1w_validate.clone('t2w_validate')

    if not precomputed.t1w_preproc:
        LOGGER.info('ANAT Stage 1: Adding T1w template workflow')
        desc += (
            'The T1-weighted (T1w) image was denoised and corrected for intensity '
            'non-uniformity (INU)'
        )

        t1w_template_wf = init_anat_template_wf(
            contrast='T1w',
            num_files=num_t1w,
            longitudinal=longitudinal,
            omp_nthreads=omp_nthreads,
            sloppy=sloppy,
            name='t1w_template_wf',
        )
        ds_t1w_template_wf = init_ds_template_wf(
            modality='T1w',
            output_dir=output_dir,
            num_anat=num_t1w,
            name='ds_t1w_template_wf',
        )
    else:
        LOGGER.info('ANAT Found preprocessed T1w - skipping Stage 1')
        desc += (
            ' A preprocessed T1w image was provided as a precomputed input and used as '
            'T1w-reference through the workflow.'
        )

        t1w_validate.inputs.in_file = precomputed.t1w_preproc
        sourcefile_buffer.inputs.source_files = [precomputed.t1w_preproc]

        workflow.connect([
            (t1w_validate, t1w_buffer, [('out_file', 't1w_preproc')]),
            (t1w_buffer, outputnode, [('t1w_preproc', 't1w_preproc')]),
        ])  # fmt:skip

    if not precomputed.t2w_preproc:
        LOGGER.info('ANAT Stage 1: Adding T2w template workflow')
        desc += (
            'The T1-weighted (T2w) image was denoised and corrected for intensity '
            'non-uniformity (INU)'
        )

        t2w_template_wf = init_anat_template_wf(
            contrast='T2w',
            num_files=num_t1w,
            longitudinal=longitudinal,
            omp_nthreads=omp_nthreads,
            sloppy=sloppy,
            name='t2w_template_wf',
        )
        ds_t2w_template_wf = init_ds_template_wf(
            modality='T2w',
            output_dir=output_dir,
            num_anat=num_t2w,
            name='ds_t2w_template_wf',
        )
    else:
        LOGGER.info('ANAT Found preprocessed T2w - skipping Stage 1')
        desc += (
            ' A preprocessed T2w image was provided as a precomputed input and used as '
            'T2w-reference through the workflow.'
        )

        t2w_validate.inputs.in_file = precomputed.t2w_preproc
        sourcefile_buffer.inputs.source_files = [precomputed.t2w_preproc]

        workflow.connect([
            (t2w_validate, t2w_buffer, [('out_file', 't2w_preproc')]),
            (t2w_buffer, outputnode, [('t2w_preproc', 't2w_preproc')]),
        ])  # fmt:skip


    # Stage 2: Use previously computed mask or calculate
    # If we only have one mask (could be either T1w/T2w),
    # just apply transform to get it in the other space
    only_t1w_mask = precomputed.t1w_mask and not precomputed.t2w_mask
    only_t2w_mask = precomputed.t2w_mask and not precomputed.t1w_mask

    if precomputed.t1w_mask or only_t2w_mask:
        desc += (
            ' A pre-computed T1w brain mask was provided as input and '
            'used throughout the workflow.'
        )
        # A mask is available and will be applied
        apply_t1w_mask = pe.Node(ApplyMask(), name='apply_t1w_mask')
        workflow.connect([
                (t1w_validate, apply_t1w_mask, [('out_file', 'in_file')]),
                (refined_buffer, outputnode, [('t1w_mask', 't1w_mask')])
            ])  # fmt:skip
        if precomputed.t1w_mask:
            LOGGER.info('ANAT Found T1w brain mask')

            t1w_buffer.inputs.t1w_mask = precomputed.t1w_mask
            # If we have a mask, always apply it
            apply_t1w_mask.inputs.in_mask = precomputed.t1w_mask
        elif only_t2w_mask:
            LOGGER.info('ANAT No T1w brain mask but a T2w mask is available')

            transform_t2w_mask = pe.Node(
                ApplyTransforms(interpolation='MultiLabel'),
                name='transform_t2w_mask'
            )
            workflow.connect([
                (refined_buffer, transform_t2w_mask, [('t2w_mask', 'input_image')]),
                (t2w_buffer, transform_t2w_mask, [('t2w_preproc', 'reference_image')]),
                (t2w2t1w_buffer, transform_t2w_mask, [('t2w2t1w_xfm', 'transforms')]),
                (transform_t2w_mask, apply_t1w_mask, [('output_image', 'in_file')]),
            ])  # fmt:skip

        if not precomputed.t1w_preproc:
            LOGGER.info('ANAT Skipping T1w skull-strip, INU-correction only')
            n4_only_wf = init_n4_only_wf(
                omp_nthreads=omp_nthreads,
                atropos_use_random_seed=not skull_strip_fixed_seed,
            )
            workflow.connect([
                (apply_t1w_mask, n4_only_wf, [('out_file', 'inputnode.in_files')]),
                (n4_only_wf, t1w_buffer, [
                    (('outputnode.bias_corrected', pop_file), 't1w_preproc'),
                    (('outputnode.out_file', pop_file), 't1w_brain'),
                ]),
            ])  # fmt:skip
        else:
            LOGGER.info('ANAT Applying T1w mask to precomputed T1w')
            workflow.connect(apply_t1w_mask, 'out_file', t1w_buffer, 't1w_brain')
    else:
        # T2w will be used for brain extraction
        # so just use the one from the coregistration workflow
        workflow.connect([
            ()
        ])

    if precomputed.t2w_mask:
        LOGGER.info('ANAT Found T2w brain mask')

        t2w_buffer.inputs.t2w_mask = precomputed.t2w_mask
        # If we have a mask, always apply it
        apply_t2w_mask = pe.Node(ApplyMask(in_mask=precomputed.t2w_mask), name='apply_t2w_mask')
        workflow.connect([
            (t2w_validate, apply_t2w_mask, [('out_file', 'in_file')]),
            (refined_buffer, outputnode, [('t2w_mask', 't2w_mask')])
        ])  # fmt:skip

        if not precomputed.t2w_preproc:
            LOGGER.info('ANAT Skipping skull-strip, INU-correction only')
            n4_only_wf = init_n4_only_wf(
                omp_nthreads=omp_nthreads,
                atropos_use_random_seed=not skull_strip_fixed_seed,
            )
            workflow.connect([
                (apply_t2w_mask, n4_only_wf, [('out_file', 'inputnode.in_files')]),
                (n4_only_wf, t2w_buffer, [
                    (('outputnode.bias_corrected', pop_file), 't2w_preproc'),
                    (('outputnode.out_file', pop_file), 't2w_brain'),
                ]),
            ])  # fmt:skip
        else:
            LOGGER.info('ANAT Applying T2w mask to precomputed T2w')
            workflow.connect(apply_t2w_mask, 'out_file', t2w_buffer, 't2w_brain')
    elif only_t1w_mask:
        workflow.connect([

        ])
    else:
        LOGGER.info('ANAT Stage 2: Preparing brain extraction workflow')
        if skull_strip_mode == 'auto':
            run_skull_strip = not all(_is_skull_stripped(img) for img in t1w)
        else:
            run_skull_strip = {'force': True, 'skip': False}[skull_strip_mode]

    # Stage 3: Coregistration
    # To use the found xfm, requires both precomputed anatomicals to be found as well
    if precomputed.t1w_preproc and precomputed.t2w_preproc and precomputed.t1w2t2w_xfm:
        LOGGER.info('ANAT Found T1w-T2w xfm')
        desc += (
            ' A T1w-T2w coregistration transform was provided as input and used throughout the workflow.'
        )

        t1w2t2w_buffer.inputs.t1w2t2w_xfm = precomputed.t1w2t2w_xfm
    else:
        LOGGER.info('ANAT Coregistering anatomicals')
        desc += (
            ' The T1w and T2w reference volumes were co-registered using ANTs.'
        )

        coregistration_wf = init_coregistration_wf(
            omp_nthreads=omp_nthreads,
            sloppy=sloppy,
            debug='registration' in config.execution.debug,
            t1w_mask=False,
            probmap=not precomputed.t2w_mask,
        )

        workflow.connect([
            (t1w_buffer, coregistration_wf, [('t1w_preproc', 'inputnode.in_t1w')]),
            (t2w_buffer, coregistration_wf, [
                ('t2w_preproc', 'inputnode.in_t2w'),
                ('t2w_mask', 'inputnode.in_mask'),
            ]),
            (coregistration_wf, t1w2t2w_buffer, [('outputnode.t1w2t2w_xfm', 't1w2t2w_xfm')]),
        ])  # fmt:skip

    # Stage 4: Segmentation
    if precomputed.t1w_dseg:
        ...


    workflow.__desc__ = desc
    return workflow


def init_infant_single_anat_fit_wf(
    modality,
    *,
    age_months: int,
    anatomicals: list,
    bids_root: str,
    precomputed,
    hires,
    longitudinal,
    omp_nthreads,
    output_dir,
    segmentation_atlases,
    skull_strip_mode,
    skull_strip_template,
    sloppy,
    spaces,
    cifti_output,
    name='infant_single_anat_fit_wf',
):
    desc = (
        '\nAnatomical data preprocessing\n\n: ',
        f'A total of {len(anatomicals)} {modality} images were found '
        'within the input BIDS dataset.\n'
    )
