import logging
import typing as ty
from pathlib import Path

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from niworkflows.interfaces.header import ValidateImage
from niworkflows.interfaces.nibabel import ApplyMask, Binarize
from niworkflows.utils.connections import pop_file
from smriprep.workflows.anatomical import init_anat_ribbon_wf, init_anat_template_wf
from smriprep.workflows.surfaces import (
    init_fsLR_reg_wf,
    init_gifti_surfaces_wf,
    init_gifti_morphometrics_wf,
    init_refinement_wf,
)
from smriprep.workflows.fit.registration import init_register_template_wf
from smriprep.workflows.outputs import (
    init_ds_dseg_wf,
    init_ds_fs_registration_wf,
    init_ds_mask_wf,
    init_ds_surface_metrics_wf,
    init_ds_surfaces_wf,
    init_ds_template_registration_wf,
    init_ds_template_wf,
    init_ds_tpms_wf,
)

from nibabies import config
from nibabies.workflows.anatomical.registration import init_coregistration_wf
from nibabies.workflows.anatomical.segmentation import init_segmentation_wf
from nibabies.workflows.anatomical.surfaces import init_mcribs_dhcp_wf

if ty.TYPE_CHECKING:
    from niworkflows.utils.spaces import Reference, SpatialReferences

    from nibabies.utils.bids import Derivatives

LOGGER = logging.getLogger('nipype.workflow')


def init_infant_anat_fit_wf(
    age_months: int,
    t1w: list,
    t2w: list,
    flair: list,
    bids_root: Path,
    precomputed: Derivatives,
    hires: bool,
    longitudinal: bool,
    omp_nthreads: int,
    output_dir: Path,
    segmentation_atlases: Path | None,
    skull_strip_mode: ty.Literal['auto', 'skip', 'force'],
    skull_strip_template: 'Reference',
    sloppy: bool,
    spaces: 'SpatialReferences',
    recon_method: ty.Literal['freesurfer', 'infantfs', 'mcribs'] | None,
    cifti_output: ty.Literal['91k', '170k'] | None,
    name: str = 'infant_anat_fit_wf',
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
    #     At the time of writing, t1w_mask / t2w_mask are an exception, which takes the form
    #     (t{1,2}w_buffer -> refined_buffer -> datasink -> outputnode)
    # All outputnode components should therefore point to files in the input or
    # output directories.
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['t1w', 't2w', 'roi', 'flair', 'subjects_dir', 'subject_id'],
        ),
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

    # Duplicates will be added to the sourcefile / anat buffers as `anat_...` to facilitate
    # usage of the preferred anatomical

    # Stage 1 inputs (filtered)
    sourcefile_buffer = pe.Node(
        niu.IdentityInterface(
            fields=['t1w_source_files', 't2w_source_files', 'anat_source_files'],
        ),
        name='sourcefile_buffer',
    )

    # Stage 2 - Anatomicals
    t1w_buffer = pe.Node(
        niu.IdentityInterface(fields=['t1w_preproc', 't1w_mask' 't1w_brain']),
        name='t1w_buffer',
    )
    t2w_buffer = pe.Node(
        niu.IdentityInterface(fields=['t2w_preproc', 't2w_mask' 't2w_brain']),
        name='t2w_buffer',
    )
    anat_buffer = pe.Node(
        niu.IdentityInterface(
            fields=[
                'anat_preproc',
                'anat_mask',
                'anat_brain',
            ]
        ),
        name='anat_buffer',
    )

    # At this point, we should decide which anatomical we will use as the reference space.
    # This will depend on the age of the participant, as myelination should be somewhat complete
    # by 9+ months
    reference_anat = 't2w' if age_months <= 8 else 't1w'
    image_type = reference_anat.capitalize()

    if reference_anat == 't1w':
        LOGGER.info('ANAT: Using T1w as the reference anatomical')
        workflow.connect([
            (t1w_buffer, anat_buffer, [
                ('t1w_preproc', 'anat_preproc'),
                ('t1w_mask', 'anat_mask'),
                ('t1w_brain', 'anat_brain'),
            ]),
        ])  # fmt:skip
    elif reference_anat == 't2w':
        LOGGER.info('ANAT: Using T2w as the reference anatomical')
        workflow.connect([
            (t2w_buffer, anat_buffer, [
                ('t2w_preproc', 'anat_preproc'),
                ('t2w_mask', 'anat_mask'),
                ('t2w_brain', 'anat_brain'),
            ]),
        ])  # fmt:skip

    # Stage 3 - Coregistration transforms
    coreg_buffer = pe.Node(
        niu.IdentityInterface(fields=['t1w2t2w_xfm', 't2w2t1w_xfm']),
        name='coreg_buffer',
    )

    # Stage 4 - Segmentation
    seg_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_dseg', 'anat_tpms', 'ants_segs']),
        name='seg_buffer',
    )
    # Stage 5 - collated template names, forward and reverse transforms
    template_buffer = pe.Node(niu.Merge(2), name='template_buffer')
    anat2std_buffer = pe.Node(niu.Merge(2), name='anat2std_buffer')
    std2anat_buffer = pe.Node(niu.Merge(2), name='std2anat_buffer')

    # Stage 6 results: Refined stage 2 results; may be direct copy if no refinement
    refined_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_mask', 'anat_brain']),
        name='refined_buffer',
    )

    fsnative_buffer = pe.Node(
        niu.IdentityInterface(fields=['fsnative2anat_xfm', 'anat2fsnative_xfm']),
        name='fsnative_buffer',
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
            ('anat_dseg', 'anat_dseg'),
            ('anat_tpms', 'anat_tpms'),
        ]),
        (anat2std_buffer, outputnode, [('out', 'anat2std_xfm')]),
        (std2anat_buffer, outputnode, [('out', 'std2anat_xfm')]),
        (template_buffer, outputnode, [('out', 'template')]),
        (sourcefile_buffer, outputnode, [
            ('t1w_source_files', 't1w_valid_list'),
            ('t2w_source_files', 't2w_valid_list'),
        ]),
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
    anat_reports_wf = init_anat_reports_wf(
        surface_recon=recon_method,
        output_dir=output_dir,
        sloppy=sloppy,
    )

    workflow.connect([
        (outputnode, anat_reports_wf, [
            ('anat_valid_list', 'inputnode.source_file'),
            ('anat_preproc', 'inputnode.anat_preproc'),
            ('anat_mask', 'inputnode.anat_mask'),
            ('anat_dseg', 'inputnode.anat_dseg'),
            ('template', 'inputnode.template'),
            ('anat2std_xfm', 'inputnode.anat2std_xfm'),
            ('subjects_dir', 'inputnode.subjects_dir'),
            ('subject_id', 'inputnode.subject_id'),
        ]),
    ])  # fmt:skip

    desc = (
        '\nAnatomical data preprocessing\n\n: '
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
            image_type='T1w',
            output_dir=output_dir,
            num_anat=num_t1w,
            name='ds_t1w_template_wf',
        )

        if reference_anat == 't1w':
            workflow.connect([
                (t1w_template_wf, sourcefile_buffer, [('outputnode.anat_valid_list', 'anat_source_files')]),
            ])  # fmt:skip

        workflow.connect([
            (inputnode, t1w_template_wf, [('t1w', 'inputnode.anat_files')]),
            (t1w_template_wf, t1w_validate, [('outputnode.anat_ref', 'in_file')]),
            (t1w_template_wf, sourcefile_buffer, [
                ('outputnode.anat_valid_list', 't1w_source_files'),
            ]),
            (t1w_template_wf, anat_reports_wf, [
                ('outputnode.out_report', 'inputnode.anat_conform_report'),
            ]),
            (t1w_template_wf, ds_t1w_template_wf, [
                ('outputnode.anat_realign_xfm', 'inputnode.anat_ref_xfms'),
            ]),
            (sourcefile_buffer, ds_t1w_template_wf, [
                ('t1w_source_files', 'inputnode.source_files'),
            ]),
            (anat_buffer, ds_t1w_template_wf, [('t1w_preproc', 'inputnode.anat_preproc')]),
            (ds_t1w_template_wf, outputnode, [('outputnode.t1w_preproc', 't1w_preproc')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Found preprocessed T1w - skipping Stage 1')
        desc += ' A preprocessed T1w image was provided as input.'

        t1w_validate.inputs.in_file = precomputed.t1w_preproc
        sourcefile_buffer.inputs.t1w_source_files = [precomputed.t1w_preproc]
        if reference_anat == 't1w':
            sourcefile_buffer.inputs.anat_source_files = [precomputed.t1w_preproc]

        workflow.connect([
            (t1w_validate, t1w_buffer, [('out_file', 't1w_preproc')]),
        ])  # fmt:skip

    if not precomputed.t2w_preproc:
        LOGGER.info('ANAT Stage 1: Adding T2w template workflow')
        desc += (
            'The T2-weighted (T2w) image was denoised and corrected for intensity '
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
            image_type='T2w',
            output_dir=output_dir,
            num_anat=num_t2w,
            name='ds_t2w_template_wf',
        )

        if reference_anat == 't2w':
            workflow.connect(
                t2w_template_wf, 'outputnode.anat_valid_list',
                sourcefile_buffer, 'anat_source_files',
            )  # fmt:skip

        workflow.connect([
            (inputnode, t2w_template_wf, [('t2w', 'inputnode.anat_files')]),
            (t2w_template_wf, t2w_validate, [('outputnode.anat_ref', 'in_file')]),
            (t2w_template_wf, sourcefile_buffer, [
                ('outputnode.anat_valid_list', 't2w_source_files'),
            ]),
            (t2w_template_wf, anat_reports_wf, [
                ('outputnode.out_report', 'inputnode.anat_conform_report'),
            ]),
            (t2w_template_wf, ds_t2w_template_wf, [
                ('outputnode.anat_realign_xfm', 'inputnode.anat_ref_xfms'),
            ]),
            (sourcefile_buffer, ds_t2w_template_wf, [
                ('t2w_source_files', 'inputnode.source_files'),
            ]),
            (anat_buffer, ds_t2w_template_wf, [('t2w_preproc', 'inputnode.anat_preproc')]),
            (ds_t2w_template_wf, outputnode, [('outputnode.t2w_preproc', 't2w_preproc')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Found preprocessed T2w - skipping Stage 1')
        desc += ' A preprocessed T2w image was provided as input.'

        t2w_validate.inputs.in_file = precomputed.t2w_preproc
        sourcefile_buffer.inputs.t2w_source_files = [precomputed.t2w_preproc]
        if precomputed.t2w_preproc:
            sourcefile_buffer.inputs.anat_source_files = [precomputed.t2w_preproc]

        workflow.connect([
            (t2w_validate, t2w_buffer, [('out_file', 't2w_preproc')]),
        ])  # fmt:skip

    # Stage 2: Use previously computed mask or calculate
    # If we only have one mask (could be either T1w/T2w),
    # just apply transform to get it in the other space
    # only_t1w_mask = precomputed.t1w_mask and not precomputed.t2w_mask
    # only_t2w_mask = precomputed.t2w_mask and not precomputed.t1w_mask

    t1w_mask = precomputed.t1w_mask
    t2w_mask = precomputed.t2w_mask
    anat_mask = None
    # T1w masking - define pre-emptively
    apply_t1w_mask = pe.Node(ApplyMask(), name='apply_t1w_mask')

    if not t1w_mask:
        if skull_strip_mode == 'auto':
            run_t1w_skull_strip = not all(_is_skull_stripped(img) for img in t1w)
        else:
            run_t1w_skull_strip = {'force': True, 'skip': False}[skull_strip_mode]

        if not run_t1w_skull_strip:  # Image is masked
            desc += (
                'The provided T1w image was previously skull-stripped; '
                'a brain mask was derived from the input image.'
            )

            if not precomputed.t1w_preproc:
                LOGGER.info('ANAT Stage 2: Skipping skull-strip, INU-correction only')

                n4_only_wf = init_n4_only_wf(
                    omp_nthreads=omp_nthreads,
                    atropos_use_random_seed=not skull_strip_fixed_seed,
                )
                workflow.connect([
                    (t1w_validate, n4_only_wf, [('out_file', 'inputnode.in_files')]),
                    (n4_only_wf, t1w_buffer, [
                        (('outputnode.bias_corrected', pop_file), 't1w_preproc'),
                        ('outputnode.out_mask', 't1w_mask'),
                        (('outputnode.out_file', pop_file), 't1w_brain'),
                        ('outputnode.out_segm', 'ants_seg'),
                    ]),
                ])  # fmt:skip
            else:
                LOGGER.info('ANAT Stage 2: Skipping skull-strip, generating mask from input')
                binarize_t1w = pe.Node(Binarize(thresh_low=2), name='binarize_t1w')
                workflow.connect([
                    (t1w_validate, binarize_t1w, [('out_file', 'in_file')]),
                    (t1w_validate, t1w_buffer, [('out_file', 't1w_brain')]),
                    (binarize_t1w, t1w_buffer, [('out_file', 't1w_mask')]),
                ])  # fmt:skip
        else:
            # T2w -> T1w transformation of the mask will occur if either
            # 1) reusing a precomputed T2w mask
            # 2) calculating with T2w template brain extraction
            LOGGER.info('ANAT T2w mask will be transformed into T1w space')
            transform_t2w_mask = pe.Node(
                ApplyTransforms(interpolation='MultiLabel'), name='transform_t2w_mask'
            )
            workflow.connect([
                (t2w_buffer, transform_t2w_mask, [('t2w_mask', 'input_image')]),
                (coreg_buffer, transform_t2w_mask, [('t2w2t1w_xfm', 'transforms')]),
                (transform_t2w_mask, apply_t1w_mask, [('output_image', 'in_mask')]),
                (t1w_buffer, apply_t1w_mask, [('t1w_preproc', 'in_file')]),  # TODO: Unsure about this connection
            ])  # fmt:skip

        # Save T1w mask
        ds_t1w_mask_wf = init_ds_mask_wf(
            bids_root=bids_root,
            output_dir=output_dir,
            mask_type='brain',
            name='ds_t1w_mask_wf',
        )
        workflow.connect([
            (sourcefile_buffer, ds_t1w_mask_wf, [('t1w_source_files', 'inputnode.source_files')]),
        ])  # fmt:skip

        if reference_anat == 't1w':
            workflow.connect([
                (refined_buffer, ds_t1w_mask_wf, [('anat_mask', 'inputnode.mask_file')]),
                (ds_t1w_mask_wf, outputnode, [('outputnode.mask_file', 'anat_mask')]),
            ])  # fmt:skip
        else:
            workflow.connect([
                (t1w_buffer, ds_t1w_mask_wf, [('t1w_mask', 'inputnode.mask_file')]),
            ])  # fmt:skip
    else:
        LOGGER.info('ANAT Found T1w brain mask')
        if reference_anat == 't1w':
            desc += (
                'A pre-computed T1w brain mask was provided as input and used throughout the '
                'workflow.'
            )
        t1w_buffer.inputs.t1w_mask = precomputed.t1w_mask
        apply_t1w_mask.inputs.in_mask = precomputed.t1w_mask
        workflow.connect(t1w_validate, 'out_file', apply_t1w_mask, 'in_file')

        if not precomputed.t1w:
            LOGGER.info('ANAT Skipping skull-strip, INU-correction only')
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
            LOGGER.info('ANAT Skipping T1w masking')
            workflow.connect(apply_t1w_mask, 'out_file', t1w_buffer, 't1w_brain')

    save_t2w_mask = True
    if precomputed.t2w_mask or only_t1w_mask:
        desc += (
            ' A pre-computed T2w brain mask was provided as input and '
            'used throughout the workflow.'
        )
        # A mask is available and will be applied
        apply_t2w_mask = pe.Node(ApplyMask(), name='apply_t2w_mask')
        workflow.connect(t2w_validate, 'out_file', apply_t2w_mask, 'in_file')
        if precomputed.t2w_mask:
            LOGGER.info('ANAT Found T2w brain mask')

            save_t2w_mask = False
            anat_buffer.inputs.t2w_mask = precomputed.t2w_mask
            apply_t1w_mask.inputs.in_mask = precomputed.t2w_mask
            workflow.connect(refined_buffer, 't2w_mask', outputnode, 't1w_mask')
        elif only_t2w_mask:
            LOGGER.info('ANAT No T2w brain mask but a T1w mask is available')

            transform_t2w_mask = pe.Node(
                ApplyTransforms(interpolation='MultiLabel'), name='transform_t1w_mask'
            )
            workflow.connect([
                (refined_buffer, transform_t2w_mask, [('t1w_mask', 'input_image')]),
                (anat_buffer, transform_t2w_mask, [('t2w_preproc', 'reference_image')]),
                (coreg_buffer, transform_t2w_mask, [('t2w2t1w_xfm', 'transforms')]),
                (transform_t2w_mask, apply_t2w_mask, [('output_image', 'in_file')]),
                (apply_t2w_mask, refined_buffer, [('out_file', 't2w_mask')]),
            ])  # fmt:skip

        if not precomputed.t2w_preproc:
            LOGGER.info('ANAT Skipping T1w skull-strip, INU-correction only')
            n4_only_wf = init_n4_only_wf(
                omp_nthreads=omp_nthreads,
                atropos_use_random_seed=not skull_strip_fixed_seed,
            )
            workflow.connect([
                (apply_t1w_mask, n4_only_wf, [('out_file', 'inputnode.in_files')]),
                (n4_only_wf, anat_buffer, [
                    (('outputnode.bias_corrected', pop_file), 't2w_preproc'),
                    (('outputnode.out_file', pop_file), 't2w_brain'),
                ]),
            ])  # fmt:skip
        else:
            LOGGER.info('ANAT Applying T2w mask to precomputed T2w')
            workflow.connect(apply_t1w_mask, 'out_file', anat_buffer, 't2w_brain')

    else:
        LOGGER.info('ANAT Stage 2: Preparing brain extraction workflow')
        if skull_strip_mode == 'auto':
            run_skull_strip = not all(_is_skull_stripped(img) for img in t1w)
        else:
            run_skull_strip = {'force': True, 'skip': False}[skull_strip_mode]
        ...

    if save_t2w_mask:
        ds_t2w_mask_wf = init_ds_mask_wf(
            bids_root=bids_root,
            output_dir=output_dir,
            mask_type='brain',
            name='ds_t2w_mask_wf',
        )
        workflow.connect([
            (sourcefile_buffer, ds_t2w_mask_wf, [('t2w_source_files', 'inputnode.source_files')]),
            (refined_buffer, ds_t2w_mask_wf, [('t2w_mask', 'inputnode.mask_file')]),
            (ds_t2w_mask_wf, outputnode, [('outputnode.mask_file', 't2w_mask')]),
        ])  # fmt:skip

    # Stage 3: Coregistration
    # To use the found xfm, requires both precomputed anatomicals to be found as well
    if precomputed.t1w_preproc and precomputed.t2w_preproc:
        if precomputed.t1w2t2w_xfm:
            LOGGER.info('ANAT Found T1w-T2w xfm')
            desc += ' A T1w-T2w coregistration transform was provided as input and used throughout the workflow.'
            coreg_buffer.inputs.t1w2t2w_xfm = precomputed.t1w2t2w_xfm
        if precomputed.t2w2t1w_xfm:
            LOGGER.info('ANAT Found T2w-T1w xfm')
            coreg_buffer.inputs.t2w2t1w_xfm = precomputed.t2w2t1w_xfm
    else:
        LOGGER.info('ANAT Coregistering anatomical references')
        desc += ' The T1w and T2w reference volumes were co-registered using ANTs.'

        coregistration_wf = init_coregistration_wf(
            omp_nthreads=omp_nthreads,
            sloppy=sloppy,
            debug='registration' in config.execution.debug,
            t1w_mask=False,
            probmap=not precomputed.t2w_mask,
        )
        workflow.connect([
            (anat_buffer, coregistration_wf, [
                ('t1w_preproc', 'inputnode.in_t1w'),
                ('t2w_preproc', 'inputnode.in_t2w'),
                ('t2w_mask', 'inputnode.in_mask'),
            ]),
            (coregistration_wf, coreg_buffer, [
                ('outputnode.t1w2t2w_xfm', 't1w2t2w_xfm'),
                ('outputnode.t2w2t1w_xfm', 't2w2t1w_xfm'),
            ]),
        ])  # fmt:skip

    # Stage 4: Segmentation
    anat_dseg = getattr(precomputed, f'{image_type}_dseg', None)
    anat_tpms = getattr(precomputed, f'{image_type}_tpms', None)
    anat_aseg = getattr(precomputed, f'{image_type}_aseg', False)
    seg_method = 'jlf' if config.execution.segmentation_atlases_dir else 'fast'

    if not (anat_dseg and anat_tpms):
        LOGGER.info('ANAT Stage 4: Tissue segmentation')
        segmentation_wf = init_segmentation_wf(
            sloppy=sloppy,
            method=seg_method,
            image_type=image_type.capitalize(),
            omp_nthreads=omp_nthreads,
            has_aseg=bool(anat_aseg),
        )

        workflow.connect([
            (anat_buffer, segmentation_wf, [(f'{image_type}_brain', 'anat_brain')]),
            (segmentation_wf, seg_buffer, [
                ('outputnode.anat_dseg', 'anat_dseg'),
                ('outputnode.anat_tpms', 'anat_tpms'),
            ]),
        ])  # fmt:skip
        if anat_aseg or seg_method == 'jlf':
            workflow.connect(segmentation_wf, 'outputnode.anat_aseg', seg_buffer, 'anat_aseg')
            if anat_aseg:
                LOGGER.info('ANAT Found precomputed anatomical segmentation')
                segmentation_wf.inputs.inputnode.anat_aseg = anat_aseg

        # TODO: datasink
        if not anat_dseg:
            ds_dseg_wf = init_ds_dseg_wf(output_dir=output_dir)
            workflow.connect([
                (sourcefile_buffer, ds_dseg_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (segmentation_wf, ds_dseg_wf, [
                    ('outputnode.anat_dseg', 'inputnode.anat_dseg'),
                ]),
                (ds_dseg_wf, seg_buffer, [('outputnode.anat_dseg', 'anat_dseg')]),
            ])  # fmt:skip

        if not anat_tpms:
            ds_tpms_wf = init_ds_tpms_wf(output_dir=output_dir)
            workflow.connect([
                (sourcefile_buffer, ds_dseg_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (segmentation_wf, ds_tpms_wf, [
                    ('outputnode.anat_tpms', 'inputnode.anat_tpms'),
                ]),
                (ds_tpms_wf, seg_buffer, [('outputnode.anat_tpms', 'anat_tpms')]),
            ])  # fmt:skip
    else:
        LOGGER.info('ANAT Stage 4: Skipping segmentation workflow')
    if anat_dseg:
        LOGGER.info('ANAT Found discrete segmentation')
        desc += 'Precomputed discrete tissue segmentations were provided as inputs.\n'
        seg_buffer.inputs.anat_dseg = anat_dseg
    if anat_tpms:
        LOGGER.info('ANAT Found tissue probability maps')
        desc += 'Precomputed tissue probabiilty maps were provided as inputs.\n'
        seg_buffer.inputs.anat_tpms = anat_tpms

    # Stage 5: Normalization
    templates = []
    found_xfms = {}
    for template in spaces.get_spaces(nonstandard=False, dim=(3,)):
        xfms = precomputed.get('transforms', {}).get(template, {})
        if set(xfms) != {'forward', 'reverse'}:
            templates.append(template)
        else:
            found_xfms[template] = xfms

    template_buffer.inputs.in1 = list(found_xfms)
    anat2std_buffer.inputs.in1 = [xfm['forward'] for xfm in found_xfms.values()]
    std2anat_buffer.inputs.in1 = [xfm['reverse'] for xfm in found_xfms.values()]

    if templates:
        LOGGER.info(f'ANAT Stage 5: Preparing normalization workflow for {templates}')
        register_template_wf = init_register_template_wf(
            sloppy=sloppy,
            omp_nthreads=omp_nthreads,
            templates=templates,
        )
        ds_template_registration_wf = init_ds_template_registration_wf(
            output_dir=str(output_dir),
            image_type=image_type.capitalize(),
        )

        workflow.connect([
            (inputnode, register_template_wf, [('roi', 'inputnode.lesion_mask')]),
            (anat_buffer, register_template_wf, [(f'{image_type}_preproc', 'inputnode.moving_image')]),
            (refined_buffer, register_template_wf, [(f'{image_type}_mask', 'inputnode.moving_mask')]),
            (sourcefile_buffer, ds_template_registration_wf, [
                (f'{image_type}_source_files', 'inputnode.source_files')
            ]),
            (register_template_wf, ds_template_registration_wf, [
                ('outputnode.template', 'inputnode.template'),
                ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
            ]),
            (register_template_wf, template_buffer, [('outputnode.template', 'in2')]),
            (ds_template_registration_wf, std2anat_buffer, [('outputnode.std2anat_xfm', 'in2')]),
            (ds_template_registration_wf, anat2std_buffer, [('outputnode.anat2std_xfm', 'in2')]),
        ])  # fmt:skip
    if found_xfms:
        LOGGER.info(f'ANAT Stage 5: Found pre-computed registrations for {found_xfms}')

    # Only refine mask if necessary
    if anat_mask or recon_method == None:
        workflow.connect([
            (anat_buffer, refined_buffer, [
                ('anat_mask', 'anat_mask'),
                ('anat_brain', 'anat_brain'),
            ]),
        ])  # fmt:skip

    workflow.__desc__ = desc

    if not recon_method:
        LOGGER.info('ANAT Skipping Stages 6+')
        return workflow

    # Stage 6: Surface reconstruction

    if recon_method == 'mcribs':
        from nibabies.workflows.anatomical.surfaces import init_mcribs_surface_recon_wf

        LOGGER.info('ANAT Stage 6: Preparing M-CRIB-S reconstruction workflow')
        surface_recon_wf = init_mcribs_surface_recon_wf(
            omp_nthreads=omp_nthreads,
            use_aseg=bool(anat_aseg),
            use_mask=True,
            precomputed=precomputed,
            mcribs_dir=str(config.execution.mcribs_dir),
        )

        workflow.connect([
            (t2w_buffer, surface_recon_wf, [
                ('t2w_preproc', 'inputnode.t2w'),
                ('t2w_mask', 'inputnode.t2w_mask'),
            ]),
            (anat_buffer, surface_recon_wf, [
                ('anat_aseg', 'inputnode.t2w_aseg'),
            ]),
        ])  # fmt:skip

    else:
        from smriprep.utils.misc import fs_isRunning

        fs_isrunning = pe.Node(
            niu.Function(function=fs_isRunning), overwrite=True, name='fs_isrunning'
        )
        fs_isrunning.inputs.logger = LOGGER

        if recon_method == 'freesurfer':
            from smriprep.workflows.surfaces import init_surface_recon_wf

            LOGGER.info('ANAT Stage 6: Preparing FreeSurfer recon-all workflow')
            fs_isrunning = pe.Node(
                niu.Function(function=fs_isRunning), overwrite=True, name='fs_isrunning'
            )
            fs_isrunning.inputs.logger = LOGGER

            surface_recon_wf = init_surface_recon_wf(
                name='surface_recon_wf',
                omp_nthreads=omp_nthreads,
                hires=hires,
                fs_no_resume=fs_no_resume,
                precomputed=precomputed,
            )

            if t2w or flair:
                t2w_or_flair = 'T2-weighted' if t2w else 'FLAIR'
                if surface_recon_wf.__desc__:
                    surface_recon_wf.__desc__ += (
                        f'A {t2w_or_flair} image was used to improve pial surface refinement.'
                    )
            workflow.connect([
                (inputnode, surface_recon_wf, [
                    ('t2w', 'inputnode.t2w'),
                    ('flair', 'inputnode.flair'),
                ]),
            ])  # fmt:skip

        elif recon_method == 'infantfs':
            from nibabies.workflows.anatomical.surfaces import init_infantfs_surface_recon_wf

            LOGGER.info('ANAT Stage 6: Preparing Infant FreeSurfer workflow')
            surface_recon_wf = init_infantfs_surface_recon_wf(
                age_months=age_months,
                precomputed=precomputed,
                use_aseg=bool(anat_mask),
            )

        # Force use of the T1w image
        workflow.connect([
            (inputnode, fs_isrunning, [
                ('subjects_dir', 'subjects_dir'),
                ('subject_id', 'subject_id'),
            ]),
            (inputnode, surface_recon_wf, [
                ('subject_id', 'inputnode.subject_id'),
            ]),
            (fs_isrunning, surface_recon_wf, [('out', 'inputnode.subjects_dir')]),
            (t1w_validate, surface_recon_wf, [('out_file', 'inputnode.t1w')]),
            (t1w_buffer, surface_recon_wf, [('t1w_brain', 'inputnode.skullstripped_t1')]),
            (surface_recon_wf, outputnode, [
                ('outputnode.subjects_dir', 'subjects_dir'),
                ('outputnode.subject_id', 'subject_id'),
            ]),
        ])  # fmt:skip

    fsnative_xfms = precomputed.get('transforms', {}).get('fsnative')
    if not fsnative_xfms:
        ds_fs_registration_wf = init_ds_fs_registration_wf(
            image_type=image_type,
            output_dir=output_dir
        )
        workflow.connect([
            (sourcefile_buffer, ds_fs_registration_wf, [
                ('anat_source_files', 'inputnode.source_files'),
            ]),
            (surface_recon_wf, ds_fs_registration_wf, [
                ('outputnode.fsnative2t1w_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (ds_fs_registration_wf, outputnode, [
                ('outputnode.fsnative2anat_xfm', 'fsnative2anat_xfm'),
            ]),
        ])  # fmt:skip
    elif 'reverse' in fsnative_xfms:
        LOGGER.info('ANAT Found fsnative-to-anatomical transform - skipping registration')
        outputnode.inputs.fsnative2anat_xfm = fsnative_xfms['reverse']
    else:
        raise RuntimeError(
            'Found an anatomical-to-fsnative transform without the reverse. Time to handle this.'
        )

    if not have_mask:
        LOGGER.info('ANAT Stage 7: Preparing mask refinement workflow')
        # Stage 6: Refine ANTs mask with FreeSurfer segmentation
        refinement_wf = init_refinement_wf()
        applyrefined = pe.Node(ApplyMask(), name='applyrefined')

        workflow.connect([
            (surface_recon_wf, refinement_wf, [
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                ('outputnode.subject_id', 'inputnode.subject_id'),
            ]),
            (fsnative_buffer, refinement_wf, [
                ('fsnative2anat_xfm', 'inputnode.fsnative2t1w_xfm'),
            ]),
            (anat_buffer, refinement_wf, [
                ('anat_preproc', 'inputnode.reference_image'),
                ('ants_seg', 'inputnode.ants_segs'),  # TODO: Verify this is the same as dseg
            ]),
            (anat_buffer, applyrefined, [('anat_preproc', 'in_file')]),
            (refinement_wf, applyrefined, [('outputnode.out_brainmask', 'in_mask')]),
            (refinement_wf, refined_buffer, [('outputnode.out_brainmask', 'anat_mask')]),
            (applyrefined, refined_buffer, [('out_file', 'anat_brain')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Found brain mask - skipping Stage 7')

    # Stages 8-10: Surface conversion and registration
    # sphere_reg is needed to generate sphere_reg_fsLR
    # sphere and sulc are needed to generate sphere_reg_msm
    # white, pial, midthickness and thickness are needed to resample in the cortical ribbon
    # TODO: Consider paring down or splitting into a subworkflow that can be called on-demand
    # A subworkflow would still need to check for precomputed outputs
    needed_anat_surfs = ['white', 'pial', 'midthickness']
    needed_metrics = ['thickness', 'sulc']
    needed_spheres = ['sphere_reg', 'sphere']

    # Detect pre-computed surfaces
    found_surfs = {
        surf: sorted(precomputed[surf])
        for surf in needed_anat_surfs + needed_metrics + needed_spheres
        if len(precomputed.get(surf, [])) == 2
    }
    if found_surfs:
        LOGGER.info(f'ANAT Stage 8: Found pre-converted surfaces for {list(found_surfs)}')
        surfaces_buffer.inputs.trait_set(**found_surfs)

    # Stage 8: Surface conversion
    surfs = [surf for surf in needed_anat_surfs if surf not in found_surfs]
    spheres = [sphere for sphere in needed_spheres if sphere not in found_surfs]
    if surfs or spheres:
        LOGGER.info(f'ANAT Stage 8: Creating GIFTI surfaces for {surfs + spheres}')
    if surfs:
        gifti_surfaces_wf = init_gifti_surfaces_wf(surfaces=surfs)
        ds_surfaces_wf = init_ds_surfaces_wf(
            bids_root=str(bids_root),
            output_dir=str(output_dir),
            surfaces=surfs,
        )

        workflow.connect([
            (surface_recon_wf, gifti_surfaces_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
            ]),
            (fsnative_buffer, gifti_surfaces_wf, [
                ('outputnode.fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (gifti_surfaces_wf, surfaces_buffer, [
                (f'outputnode.{surf}', surf) for surf in surfs
            ]),
            (sourcefile_buffer, ds_surfaces_wf, [('source_files', 'inputnode.source_files')]),
            (gifti_surfaces_wf, ds_surfaces_wf, [
                (f'outputnode.{surf}', f'inputnode.{surf}') for surf in surfs
            ]),
        ])  # fmt:skip
    if spheres:
        gifti_spheres_wf = init_gifti_surfaces_wf(
            surfaces=spheres, to_scanner=False, name='gifti_spheres_wf'
        )
        ds_spheres_wf = init_ds_surfaces_wf(
            bids_root=str(bids_root),
            output_dir=str(output_dir),
            surfaces=spheres,
            name='ds_spheres_wf',
        )

        workflow.connect([
            (surface_recon_wf, gifti_spheres_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                # No transform for spheres, following HCP pipelines' lead
            ]),
            (gifti_spheres_wf, surfaces_buffer, [
                (f'outputnode.{sphere}', sphere) for sphere in spheres
            ]),
            (sourcefile_buffer, ds_spheres_wf, [('source_files', 'inputnode.source_files')]),
            (gifti_spheres_wf, ds_spheres_wf, [
                (f'outputnode.{sphere}', f'inputnode.{sphere}') for sphere in spheres
            ]),
        ])  # fmt:skip
    metrics = [metric for metric in needed_metrics if metric not in found_surfs]
    if metrics:
        LOGGER.info(f'ANAT Stage 8: Creating GIFTI metrics for {metrics}')
        gifti_morph_wf = init_gifti_morphometrics_wf(morphometrics=metrics)
        ds_morph_wf = init_ds_surface_metrics_wf(
            bids_root=str(bids_root),
            output_dir=str(output_dir),
            metrics=metrics,
            name='ds_morph_wf',
        )

        workflow.connect([
            (surface_recon_wf, gifti_morph_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
            ]),
            (gifti_morph_wf, surfaces_buffer, [
                (f'outputnode.{metric}', metric) for metric in metrics
            ]),
            (sourcefile_buffer, ds_morph_wf, [('source_files', 'inputnode.source_files')]),
            (gifti_morph_wf, ds_morph_wf, [
                (f'outputnode.{metric}', f'inputnode.{metric}') for metric in metrics
            ]),
        ])  # fmt:skip

    if 'anat_ribbon' not in precomputed:
        LOGGER.info('ANAT Stage 8a: Creating cortical ribbon mask')
        anat_ribbon_wf = init_anat_ribbon_wf()
        ds_ribbon_mask_wf = init_ds_mask_wf(
            bids_root=str(bids_root),
            output_dir=str(output_dir),
            mask_type='ribbon',
            name='ds_ribbon_mask_wf',
        )

        workflow.connect([
            (anat_buffer, anat_ribbon_wf, [
                ('anat_preproc', 'inputnode.ref_file'),
            ]),
            (surfaces_buffer, anat_ribbon_wf, [
                ('white', 'inputnode.white'),
                ('pial', 'inputnode.pial'),
            ]),
            (sourcefile_buffer, ds_ribbon_mask_wf, [('source_files', 'inputnode.source_files')]),
            (anat_ribbon_wf, ds_ribbon_mask_wf, [
                ('outputnode.anat_ribbon', 'inputnode.mask_file'),
            ]),
            (ds_ribbon_mask_wf, outputnode, [('outputnode.mask_file', 'anat_ribbon')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Stage 8a: Found pre-computed cortical ribbon mask')
        outputnode.inputs.anat_ribbon = precomputed['anat_ribbon']

    # Stage 9: Baseline fsLR registration
    if len(precomputed.get('sphere_reg_fsLR', [])) < 2:
        LOGGER.info('ANAT Stage 9: Creating fsLR registration sphere')
        if recon_method == 'mcribs':
            fsLR_reg_wf = init_mcribs_dhcp_wf()
        else:
            fsLR_reg_wf = init_fsLR_reg_wf()

        ds_fsLR_reg_wf = init_ds_surfaces_wf(
            bids_root=str(bids_root),
            output_dir=str(output_dir),
            surfaces=['sphere_reg_fsLR'],
            name='ds_fsLR_reg_wf',
        )

        workflow.connect([
            (surfaces_buffer, fsLR_reg_wf, [('sphere_reg', 'inputnode.sphere_reg')]),
            (sourcefile_buffer, ds_fsLR_reg_wf, [('source_files', 'inputnode.source_files')]),
            (fsLR_reg_wf, ds_fsLR_reg_wf, [
                ('outputnode.sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR')
            ]),
            (ds_fsLR_reg_wf, fsLR_buffer, [('outputnode.sphere_reg_fsLR', 'sphere_reg_fsLR')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Stage 9: Found pre-computed fsLR registration sphere')
        fsLR_buffer.inputs.sphere_reg_fsLR = sorted(precomputed['sphere_reg_fsLR'])
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
        'within the input BIDS dataset.\n',
    )
