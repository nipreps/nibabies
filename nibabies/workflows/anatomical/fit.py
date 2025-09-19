import logging
import typing as ty
from pathlib import Path

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.anat.ants import init_n4_only_wf
from niworkflows.engine import Workflow, tag
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from niworkflows.interfaces.header import ValidateImage
from niworkflows.interfaces.nibabel import ApplyMask, Binarize
from niworkflows.interfaces.utility import KeySelect
from niworkflows.utils.connections import pop_file
from smriprep.workflows.anatomical import (
    _is_skull_stripped,
    init_anat_template_wf,
)
from smriprep.workflows.fit.registration import init_register_template_wf
from smriprep.workflows.outputs import (
    init_ds_dseg_wf,
    init_ds_fs_registration_wf,
    init_ds_mask_wf,
    init_ds_surface_masks_wf,
    init_ds_surface_metrics_wf,
    init_ds_surfaces_wf,
    init_ds_template_registration_wf,
    init_ds_template_wf,
    init_ds_tpms_wf,
)
from smriprep.workflows.surfaces import (
    init_anat_ribbon_wf,
    init_cortex_masks_wf,
    init_fsLR_reg_wf,
    init_gifti_morphometrics_wf,
    init_gifti_surfaces_wf,
    init_refinement_wf,
)

from nibabies import config
from nibabies.interfaces import DerivativesDataSink
from nibabies.workflows.anatomical.brain_extraction import init_infant_brain_extraction_wf
from nibabies.workflows.anatomical.outputs import init_anat_reports_wf, init_coreg_report_wf
from nibabies.workflows.anatomical.preproc import (
    init_anat_preproc_wf,
    init_conform_derivative_wf,
    init_csf_norm_wf,
)
from nibabies.workflows.anatomical.registration import (
    init_concat_registrations_wf,
    init_coregistration_wf,
)
from nibabies.workflows.anatomical.segmentation import init_segmentation_wf
from nibabies.workflows.anatomical.surfaces import init_mcribs_dhcp_wf

if ty.TYPE_CHECKING:
    from niworkflows.utils.spaces import Reference, SpatialReferences


LOGGER = logging.getLogger('nipype.workflow')


@tag('anat.fit')
def init_infant_anat_fit_wf(
    *,
    age_months: int,
    t1w: list,
    t2w: list,
    flair: list,
    bids_root: str,
    precomputed: dict,
    longitudinal: bool,
    omp_nthreads: int,
    output_dir: str,
    segmentation_atlases: str | Path | None,
    skull_strip_mode: ty.Literal['auto', 'skip', 'force'],
    skull_strip_template: 'Reference',
    skull_strip_fixed_seed: bool,
    sloppy: bool,
    spaces: 'SpatialReferences',
    recon_method: ty.Literal['freesurfer', 'infantfs', 'mcribs'] | None,
    reference_anat: ty.Literal['T1w', 'T2w'],
    cifti_output: ty.Literal['91k', '170k', False],
    msm_sulc: bool = False,
    name: str = 'infant_anat_fit_wf',
) -> Workflow:
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

    anat = reference_anat.lower()

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
                # T1w
                't1w_preproc',
                't1w_mask',
                't1w_valid_list',
                # T2w
                't2w_preproc',
                't2w_mask',
                't2w_valid_list',
                # Anat specific
                'anat_preproc',
                'anat_mask',
                'anat_dseg',
                'anat_tpms',
                'anat2std_xfm',
                'fsnative2anat_xfm',
                't1w2t2w_xfm',
                't2w2t1w_xfm',
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
                'cortex_mask',
                'anat_ribbon',
                # Reverse transform; not computable from forward transform
                'std2anat_xfm',
                # Metadata
                'template',
                'subjects_dir',
                'subject_id',
                'anat_valid_list',
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
        niu.IdentityInterface(fields=['t1w_preproc', 't1w_mask', 't1w_brain']),
        name='t1w_buffer',
    )
    t2w_buffer = pe.Node(
        niu.IdentityInterface(fields=['t2w_preproc', 't2w_mask', 't2w_brain', 't2w_probmap']),
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

    # Additional buffer if CSF normalization is used
    anat_reg_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_preproc']),
        name='anat_reg_buffer',
    )
    if not config.workflow.norm_csf:
        workflow.connect(anat_buffer, 'anat_preproc', anat_reg_buffer, 'anat_preproc')

    if reference_anat == 'T1w':
        LOGGER.info('ANAT: Using T1w as the reference anatomical')
        workflow.connect([
            (t1w_buffer, anat_buffer, [
                ('t1w_preproc', 'anat_preproc'),
                ('t1w_mask', 'anat_mask'),
                ('t1w_brain', 'anat_brain'),
            ]),
        ])  # fmt:skip
    elif reference_anat == 'T2w':
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

    aseg_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_aseg']),
        name='aseg_buffer',
    )

    # Stage 4 - Segmentation
    seg_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_dseg', 'anat_tpms', 'ants_segs']),
        name='seg_buffer',
    )
    # Stage 5 - collated template names, forward and reverse transforms
    template_buffer = pe.Node(niu.Merge(3), name='template_buffer')
    anat2std_buffer = pe.Node(niu.Merge(3), name='anat2std_buffer')
    std2anat_buffer = pe.Node(niu.Merge(3), name='std2anat_buffer')

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
        (anat_buffer, outputnode, [
            ('anat_preproc', 'anat_preproc'),
        ]),
        (refined_buffer, outputnode, [
            ('anat_mask', 'anat_mask'),
            ('anat_brain', 'anat_brain'),
        ]),
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
            ('anat_source_files', 'anat_valid_list'),  # Alias for reference anat files
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
        spaces=spaces,
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

    t1w_preproc = precomputed.get('t1w_preproc')
    t2w_preproc = precomputed.get('t2w_preproc')
    anat_preproc = precomputed.get(f'{anat}_preproc')

    # Stage 1: Conform & valid T1w/T2w images
    # Note: Since stage 1 & 2 are tightly knit together, it may be more intuitive
    # to combine them in a separate workflow.
    t1w_validate = pe.Node(ValidateImage(), name='t1w_validate', run_without_submitting=True)
    t2w_validate = t1w_validate.clone('t2w_validate')

    if not t1w_preproc:
        LOGGER.info('ANAT Stage 1: Adding T1w template workflow')
        desc += (
            'The T1-weighted (T1w) image was denoised and corrected for intensity '
            'non-uniformity (INU)'
        )

        t1w_template_wf = init_anat_template_wf(
            image_type='T1w',
            num_files=num_t1w,
            longitudinal=longitudinal,
            omp_nthreads=omp_nthreads,
            name='t1w_template_wf',
        )
        ds_t1w_template_wf = init_ds_template_wf(
            image_type='T1w',
            output_dir=output_dir,
            num_anat=num_t1w,
            name='ds_t1w_template_wf',
        )

        if reference_anat == 'T1w':
            workflow.connect([
                (t1w_template_wf, sourcefile_buffer, [
                    ('outputnode.anat_valid_list', 'anat_source_files'),
                ]),
                (t1w_template_wf, anat_reports_wf, [
                    ('outputnode.out_report', 'inputnode.anat_conform_report'),
                ]),
            ])  # fmt:skip

        workflow.connect([
            (inputnode, t1w_template_wf, [('t1w', 'inputnode.anat_files')]),
            (t1w_template_wf, t1w_validate, [('outputnode.anat_ref', 'in_file')]),
            (t1w_template_wf, sourcefile_buffer, [
                ('outputnode.anat_valid_list', 't1w_source_files'),
            ]),
            (t1w_template_wf, ds_t1w_template_wf, [
                ('outputnode.anat_realign_xfm', 'inputnode.anat_ref_xfms'),
            ]),
            (sourcefile_buffer, ds_t1w_template_wf, [
                ('t1w_source_files', 'inputnode.source_files'),
            ]),
            (t1w_buffer, ds_t1w_template_wf, [('t1w_preproc', 'inputnode.anat_preproc')]),
            (ds_t1w_template_wf, outputnode, [('outputnode.anat_preproc', 't1w_preproc')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Found preprocessed T1w - skipping Stage 1')
        desc += ' A preprocessed T1w image was provided as input.'

        t1w_validate.inputs.in_file = t1w_preproc
        sourcefile_buffer.inputs.t1w_source_files = [t1w_preproc]
        if reference_anat == 'T1w':
            sourcefile_buffer.inputs.anat_source_files = [t1w_preproc]

        workflow.connect([
            (t1w_validate, t1w_buffer, [('out_file', 't1w_preproc')]),
        ])  # fmt:skip

    if not t2w_preproc:
        LOGGER.info('ANAT Stage 1: Adding T2w template workflow')
        desc += (
            'The T2-weighted (T2w) image was denoised and corrected for intensity '
            'non-uniformity (INU)'
        )

        t2w_template_wf = init_anat_template_wf(
            image_type='T2w',
            num_files=num_t1w,
            longitudinal=longitudinal,
            omp_nthreads=omp_nthreads,
            name='t2w_template_wf',
        )
        ds_t2w_template_wf = init_ds_template_wf(
            image_type='T2w',
            output_dir=output_dir,
            num_anat=num_t2w,
            name='ds_t2w_template_wf',
        )

        if reference_anat == 'T2w':
            workflow.connect([
                (t2w_template_wf, sourcefile_buffer, [
                    ('outputnode.anat_valid_list', 'anat_source_files'),
                ]),
                (t2w_template_wf, anat_reports_wf, [
                    ('outputnode.out_report', 'inputnode.anat_conform_report'),
                ]),
            ])  # fmt:skip

        workflow.connect([
            (inputnode, t2w_template_wf, [('t2w', 'inputnode.anat_files')]),
            (t2w_template_wf, t2w_validate, [('outputnode.anat_ref', 'in_file')]),
            (t2w_template_wf, sourcefile_buffer, [
                ('outputnode.anat_valid_list', 't2w_source_files'),
            ]),
            (t2w_template_wf, ds_t2w_template_wf, [
                ('outputnode.anat_realign_xfm', 'inputnode.anat_ref_xfms'),
            ]),
            (sourcefile_buffer, ds_t2w_template_wf, [
                ('t2w_source_files', 'inputnode.source_files'),
            ]),
            (t2w_buffer, ds_t2w_template_wf, [('t2w_preproc', 'inputnode.anat_preproc')]),
            (ds_t2w_template_wf, outputnode, [('outputnode.anat_preproc', 't2w_preproc')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Found preprocessed T2w - skipping Stage 1')
        desc += ' A preprocessed T2w image was provided as input.'

        t2w_validate.inputs.in_file = t2w_preproc
        sourcefile_buffer.inputs.t2w_source_files = [t2w_preproc]
        if reference_anat == 'T2w':
            sourcefile_buffer.inputs.anat_source_files = [t2w_preproc]

        workflow.connect([
            (t2w_validate, t2w_buffer, [('out_file', 't2w_preproc')]),
        ])  # fmt:skip

    # Stage 2: Use previously computed mask or calculate
    # If we only have one mask (could be either T1w/T2w),
    # just apply transform to get it in the other space
    t1w_mask = precomputed.get('t1w_mask')
    t2w_mask = precomputed.get('t2w_mask')
    anat_mask = precomputed.get(f'{anat}_mask')
    refine_mask = False
    # T1w masking - define preemptively
    apply_t1w_mask = pe.Node(ApplyMask(), name='apply_t1w_mask')
    apply_t2w_mask = apply_t1w_mask.clone(name='apply_t2w_mask')

    # T1w masking logic:
    #
    # PCM = Pre-computed mask
    # PCT = Pre-computed template
    # SS = Skull stripping required
    #
    # PCM, PCT, SS -> Apply PCM to PCT
    # PCM, PCT, !SS -> Apply PCM to PCT
    # PCM, !PCT, SS -> Apply PCM to template, then run INU
    # PCM, !PCT, !SS -> Apply PCM to template, then run INU
    # !PCM, PCT, SS -> Transform T2w brain mask with t2w2t1w_xfm
    # !PCM, PCT, !SS -> Binarize PCT
    # !PCM, !PCT, SS -> Transform T2w brain mask with t2w2t1w_xfm
    # !PCM, !PCT, !SS -> INU correct template

    if not t1w_mask:
        if skull_strip_mode == 'auto':
            run_t1w_skull_strip = not all(_is_skull_stripped(img) for img in t1w)
        else:
            run_t1w_skull_strip = {'force': True, 'skip': False}[skull_strip_mode]

        if not run_t1w_skull_strip:
            desc += (
                'The provided T1w image was previously skull-stripped; '
                'a brain mask was derived from the input image.'
            )

            if not t1w_preproc:
                LOGGER.info('ANAT Stage 2: Skipping skull-strip, INU-correction only')

                t1w_n4_only_wf = init_n4_only_wf(
                    omp_nthreads=omp_nthreads,
                    atropos_use_random_seed=not skull_strip_fixed_seed,
                )
                workflow.connect([
                    (t1w_validate, t1w_n4_only_wf, [('out_file', 'inputnode.in_files')]),
                    (t1w_n4_only_wf, t1w_buffer, [
                        (('outputnode.bias_corrected', pop_file), 't1w_preproc'),
                        ('outputnode.out_mask', 't1w_mask'),
                        (('outputnode.out_file', pop_file), 't1w_brain'),
                    ]),
                ])  # fmt:skip

                if reference_anat == 'T1w':
                    refine_mask = True
                    workflow.connect([
                        (t1w_n4_only_wf, seg_buffer, [
                            ('outputnode.out_segm', 'ants_segs'),
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
            # TODO: May need to differentiate presence of t1w_preproc
            workflow.connect([
                (t2w_buffer, transform_t2w_mask, [('t2w_mask', 'input_image')]),
                (coreg_buffer, transform_t2w_mask, [('t2w2t1w_xfm', 'transforms')]),
                (t1w_validate, transform_t2w_mask, [('out_file', 'reference_image')]),

                (transform_t2w_mask, t1w_buffer, [('output_image', 't1w_mask')]),
                (transform_t2w_mask, apply_t1w_mask, [('output_image', 'in_mask')]),
                (t1w_validate, apply_t1w_mask, [('out_file', 'in_file')]),
                (apply_t1w_mask, t1w_buffer, [('out_file', 't1w_brain')]),
            ])  # fmt:skip

            if not t1w_preproc:
                t1w_preproc_wf = init_anat_preproc_wf(name='t1w_preproc_wf')
                workflow.connect([
                    (t1w_validate, t1w_preproc_wf, [
                        ('out_file', 'inputnode.in_anat'),
                    ]),
                    (t1w_preproc_wf, t1w_buffer, [
                        ('outputnode.anat_preproc', 't1w_preproc'),
                    ]),
                ])  # fmt:skip

        # Save T1w mask
        ds_t1w_mask_wf = init_ds_mask_wf(
            bids_root=bids_root,
            output_dir=output_dir,
            mask_type='brain',
            extra_entities={'space': 'T1w'},
            name='ds_t1w_mask_wf',
        )
        workflow.connect([
            (sourcefile_buffer, ds_t1w_mask_wf, [('t1w_source_files', 'inputnode.source_files')]),
        ])  # fmt:skip

        if reference_anat == 'T1w':
            workflow.connect([
                (refined_buffer, ds_t1w_mask_wf, [('anat_mask', 'inputnode.mask_file')]),
            ])  # fmt:skip
        else:
            workflow.connect([
                (t1w_buffer, ds_t1w_mask_wf, [('t1w_mask', 'inputnode.mask_file')]),
            ])  # fmt:skip
    else:  # T1w mask derivative has been found
        LOGGER.info('ANAT Found T1w brain mask')
        if reference_anat == 'T1w':
            desc += (
                'A pre-computed T1w brain mask was provided as input and used throughout the '
                'workflow.'
            )
        workflow.connect(apply_t1w_mask, 'out_file', t1w_buffer, 't1w_brain')

        if not t1w_preproc:
            # Ensure compatibility with T1w template
            conform_t1w_mask_wf = init_conform_derivative_wf(
                in_file=t1w_mask, name='conform_t1w_mask_wf'
            )

            LOGGER.info('ANAT Skipping skull-strip, INU-correction only')
            t1w_n4_wf = init_anat_preproc_wf(name='t1w_n4_wf')
            workflow.connect([
                (t1w_validate, conform_t1w_mask_wf, [('out_file', 'inputnode.ref_file')]),
                (conform_t1w_mask_wf, t1w_buffer, [('outputnode.out_file', 't1w_mask')]),
                (conform_t1w_mask_wf, apply_t1w_mask, [('outputnode.out_file', 'in_mask')]),
                (t1w_validate, t1w_n4_wf, [('out_file', 'inputnode.in_anat')]),
                (t1w_n4_wf, t1w_buffer, [('outputnode.anat_preproc', 't1w_preproc')]),
                (t1w_n4_wf, apply_t1w_mask, [('outputnode.anat_preproc', 'in_file')]),
            ])  # fmt:skip
        else:
            LOGGER.info('ANAT Skipping T1w masking')
            workflow.connect(t1w_validate, 'out_file', apply_t1w_mask, 'in_file')
            t1w_buffer.inputs.t1w_mask = t1w_mask
            apply_t1w_mask.inputs.in_mask = t1w_mask

    # T2w masking logic:
    #
    # PCM-T2 = Pre-computed mask
    # PCM-T1 = Pre-computed mask for PCT-T1
    # PCT = Pre-computed template
    # SS = Skull stripping required
    #
    # PCM, PCT, SS -> Apply PCM to PCT
    # PCM, PCT, !SS -> Apply PCM to PCT
    # PCM, !PCT, SS -> Apply PCM to template, then run INU
    # PCM, !PCT, !SS -> Apply PCM to template, then run INU
    # !PCM, PCT, SS -> Run brain extraction
    # !PCM, PCT, !SS -> Binarize PCT
    # !PCM, !PCT, SS -> Run brain extraction
    # !PCM, !PCT, !SS -> INU correct template
    if not t2w_mask:
        if skull_strip_mode == 'auto':
            run_t2w_skull_strip = not all(_is_skull_stripped(img) for img in t2w)
        else:
            run_t2w_skull_strip = {'force': True, 'skip': False}[skull_strip_mode]

        if not run_t2w_skull_strip:
            desc += (
                'The T2w reference was previously skull-stripped; '
                'a brain mask was derived from the input image.'
            )

            if not t2w_preproc:
                LOGGER.info('ANAT Stage 2b: Skipping skull-strip, INU-correction only')

                t2w_n4_only_wf = init_n4_only_wf(
                    omp_nthreads=omp_nthreads,
                    bids_suffix=reference_anat,
                    atropos_use_random_seed=not skull_strip_fixed_seed,
                    name='t2w_n4_only_wf',
                )
                workflow.connect([
                    (t2w_validate, t2w_n4_only_wf, [('out_file', 'inputnode.in_files')]),
                    (t2w_n4_only_wf, t2w_buffer, [
                        (('outputnode.bias_corrected', pop_file), 't2w_preproc'),
                        ('outputnode.out_mask', 't2w_mask'),
                        (('outputnode.out_file', pop_file), 't2w_brain'),
                    ]),
                ])  # fmt:skip
                if reference_anat == 'T2w':
                    refine_mask = True
                    workflow.connect([
                        (t2w_n4_only_wf, seg_buffer, [
                            ('outputnode.out_segm', 'ants_segs'),
                        ]),
                    ])  # fmt:skip
            else:
                LOGGER.info('ANAT Stage 2b: Skipping skull-strip, generating mask from input')
                binarize_t2w = pe.Node(Binarize(thresh_low=2), name='binarize_t2w')
                workflow.connect([
                    (t2w_validate, binarize_t2w, [('out_file', 'in_file')]),
                    (t2w_validate, t2w_buffer, [('out_file', 't2w_brain')]),
                    (binarize_t2w, t2w_buffer, [('out_file', 't2w_mask')]),
                ])  # fmt:skip
        else:
            LOGGER.info('ANAT Atlas-based brain mask will be calculated on the T2w')
            brain_extraction_wf = init_infant_brain_extraction_wf(
                omp_nthreads=omp_nthreads,
                sloppy=sloppy,
                age_months=age_months,
                ants_affine_init=True,
                skull_strip_template=skull_strip_template.space,
                template_specs=skull_strip_template.spec,
                debug='registration' in config.execution.debug,
            )

            workflow.connect([
                (t2w_validate, brain_extraction_wf, [
                    ('out_file', 'inputnode.t2w_preproc'),
                ]),
                (brain_extraction_wf, t2w_buffer, [
                    ('outputnode.out_mask', 't2w_mask'),
                    ('outputnode.t2w_brain', 't2w_brain'),
                    ('outputnode.t2w_preproc', 't2w_preproc'),
                    ('outputnode.out_probmap', 't2w_probmap')
                ]),
            ])  # fmt:skip

        # Save T2w mask
        ds_t2w_mask_wf = init_ds_mask_wf(
            bids_root=bids_root,
            output_dir=output_dir,
            mask_type='brain',
            extra_entities={'space': 'T2w'},
            name='ds_t2w_mask_wf',
        )
        workflow.connect([
            (sourcefile_buffer, ds_t2w_mask_wf, [('t2w_source_files', 'inputnode.source_files')]),
        ])  # fmt:skip

        if reference_anat == 'T2w':
            workflow.connect([
                (refined_buffer, ds_t2w_mask_wf, [('anat_mask', 'inputnode.mask_file')]),
            ])  # fmt:skip
        else:
            workflow.connect([
                (t2w_buffer, ds_t2w_mask_wf, [('t2w_mask', 'inputnode.mask_file')]),
            ])  # fmt:skip
    else:
        LOGGER.info('ANAT Found T2w brain mask')
        if reference_anat == 'T2w':
            desc += (
                'A pre-computed T2w brain mask was provided as input and used throughout the '
                'workflow.'
            )
        workflow.connect(apply_t2w_mask, 'out_file', t2w_buffer, 't2w_brain')

        if not t2w_preproc:
            # Ensure compatibility with T2w template
            conform_t2w_mask_wf = init_conform_derivative_wf(
                in_file=t2w_mask,
                name='conform_t2w_mask_wf',
            )

            LOGGER.info('ANAT Skipping skull-strip, INU-correction only')
            t2w_n4_wf = init_anat_preproc_wf(name='t2w_n4_wf')
            workflow.connect([
                (t2w_validate, conform_t2w_mask_wf, [('out_file', 'inputnode.ref_file')]),
                (conform_t2w_mask_wf, t2w_buffer, [('outputnode.out_file', 't2w_mask')]),
                (conform_t2w_mask_wf, apply_t2w_mask, [('outputnode.out_file', 'in_mask')]),
                (t2w_validate, t2w_n4_wf, [('out_file', 'inputnode.in_anat')]),
                (t2w_n4_wf, t2w_buffer, [('outputnode.anat_preproc', 't2w_preproc')]),
                (t2w_n4_wf, apply_t2w_mask, [('outputnode.anat_preproc', 'in_file')]),
            ])  # fmt:skip
        else:
            LOGGER.info('ANAT Skipping T2w masking')
            workflow.connect(t2w_validate, 'out_file', apply_t2w_mask, 'in_file')
            t2w_buffer.inputs.t2w_mask = t2w_mask
            apply_t2w_mask.inputs.in_mask = t2w_mask

    # Stage 3: Coregistration
    t1w2t2w_xfm = precomputed.get('t1w2t2w_xfm')
    t2w2t1w_xfm = precomputed.get('t2w2t1w_xfm')

    # To use the found xfm, requires both precomputed anatomicals to be found as well
    if (t1w_preproc and t2w_preproc) and (t1w2t2w_xfm and t2w2t1w_xfm):
        LOGGER.info('ANAT Found T1w<->T2w xfms')
        desc += (
            ' A T1w-T2w coregistration transform was provided as input and used throughout '
            'the workflow.'
        )
        coreg_buffer.inputs.t1w2t2w_xfm = t1w2t2w_xfm
        coreg_buffer.inputs.t2w2t1w_xfm = t2w2t1w_xfm
    else:
        LOGGER.info('ANAT Coregistering anatomicals')
        desc += ' The T1w and T2w reference volumes were co-registered using ANTs.'

        probmap = not t2w_preproc and not t2w_mask
        coregistration_wf = init_coregistration_wf(
            omp_nthreads=omp_nthreads,
            sloppy=sloppy,
            debug='registration' in config.execution.debug,
            t1w_mask=False,
            probmap=probmap,
        )

        ds_t1w2t2w_xfm = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                to='T2w',
                mode='image',
                suffix='xfm',
                dismiss_entities=('desc', 'echo'),
                **{'from': 'T1w'},
            ),
            name='ds_t1w2t2w_xfm',
            run_without_submitting=True,
        )

        ds_t2w2t1w_xfm = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                to='T1w',
                mode='image',
                suffix='xfm',
                dismiss_entities=('desc', 'echo'),
                **{'from': 'T2w'},
            ),
            name='ds_t2w2t1w_xfm',
            run_without_submitting=True,
        )

        coreg_report_wf = init_coreg_report_wf(output_dir=output_dir)

        workflow.connect([
            (t1w_validate, coregistration_wf, [
                ('out_file', 'inputnode.in_t1w'),
            ]),
            (t2w_buffer, coregistration_wf, [
                ('t2w_preproc', 'inputnode.in_t2w'),
                ('t2w_mask', 'inputnode.in_mask'),
            ]),
            (coregistration_wf, ds_t1w2t2w_xfm, [
                ('outputnode.t1w2t2w_xfm', 'in_file'),
            ]),
            (sourcefile_buffer, ds_t1w2t2w_xfm, [
                ('t1w_source_files', 'source_file'),
            ]),
            (coregistration_wf, ds_t2w2t1w_xfm, [
                ('outputnode.t2w2t1w_xfm', 'in_file'),
            ]),
            (sourcefile_buffer, ds_t2w2t1w_xfm, [
                ('t2w_source_files', 'source_file'),
            ]),
            (coregistration_wf, coreg_buffer, [
                ('outputnode.t1w2t2w_xfm', 't1w2t2w_xfm'),
                ('outputnode.t2w2t1w_xfm', 't2w2t1w_xfm'),
            ]),
            (coregistration_wf, coreg_report_wf, [
                ('outputnode.t1w_preproc', 'inputnode.t1w_preproc'),
                ('outputnode.t2w_preproc', 'inputnode.t2w_preproc'),
                ('outputnode.t1w_mask', 'inputnode.in_mask'),
            ]),
            (sourcefile_buffer, coreg_report_wf, [('anat_source_files', 'inputnode.source_file')]),
        ])  # fmt:skip

        if probmap:
            workflow.connect([
                (t2w_buffer, coregistration_wf, [
                    ('t2w_probmap', 'inputnode.in_probmap'),
                ])
            ])  # fmt:skip

    # Stage 4: Segmentation
    anat_dseg = precomputed.get(f'{anat}_dseg')
    anat_tpms = precomputed.get(f'{anat}_tpms')
    anat_aseg = precomputed.get(f'{anat}_aseg')

    seg_method = 'jlf' if config.execution.segmentation_atlases_dir else 'fast'

    if anat_aseg:
        LOGGER.info('ANAT Found precomputed anatomical segmentation')
        # Ensure compatibility with anatomical template
        if not anat_preproc:
            conform_aseg_wf = init_conform_derivative_wf(
                in_file=anat_aseg,
                name='conform_aseg_wf',
            )

            workflow.connect([
                (anat_buffer, conform_aseg_wf, [('anat_preproc', 'inputnode.ref_file')]),
                (conform_aseg_wf, aseg_buffer, [('outputnode.out_file', 'anat_aseg')]),
            ])  # fmt:skip
        else:
            aseg_buffer.inputs.anat_aseg = anat_aseg

    if not (anat_dseg and anat_tpms):
        LOGGER.info('ANAT Stage 4: Tissue segmentation')
        segmentation_wf = init_segmentation_wf(
            sloppy=sloppy,
            method=seg_method,
            image_type=reference_anat,
            omp_nthreads=omp_nthreads,
            has_aseg=bool(anat_aseg),
        )

        workflow.connect([
            (anat_buffer, segmentation_wf, [('anat_brain', 'inputnode.anat_brain')]),
            (segmentation_wf, seg_buffer, [
                ('outputnode.anat_dseg', 'anat_dseg'),
                ('outputnode.anat_tpms', 'anat_tpms'),
            ]),
        ])  # fmt:skip

        if anat_aseg:
            workflow.connect(aseg_buffer, 'anat_aseg', segmentation_wf, 'inputnode.anat_aseg')
        elif seg_method == 'jlf':
            workflow.connect(segmentation_wf, 'outputnode.anat_aseg', aseg_buffer, 'anat_aseg')
            # TODO: datasink aseg

        if not anat_dseg:
            ds_dseg_wf = init_ds_dseg_wf(
                output_dir=str(output_dir),
                extra_entities={'space': reference_anat},
            )
            workflow.connect([
                (sourcefile_buffer, ds_dseg_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (segmentation_wf, ds_dseg_wf, [
                    ('outputnode.anat_dseg', 'inputnode.anat_dseg'),
                ]),
            ])  # fmt:skip

        if not anat_tpms:
            ds_tpms_wf = init_ds_tpms_wf(
                output_dir=str(output_dir),
                extra_entities={'space': reference_anat},
            )
            workflow.connect([
                (sourcefile_buffer, ds_tpms_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (segmentation_wf, ds_tpms_wf, [
                    ('outputnode.anat_tpms', 'inputnode.anat_tpms'),
                ]),
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

    # If selected adult templates are requested (MNI152 6th Gen or 2009)
    # opt to concatenate transforms first from native -> infant template (MNIInfant),
    # and then use a previously computed MNIInfant<cohort> -> MNI transform
    # this minimizes the chance of a bad registration.

    templates = []
    concat_xfms = []
    found_xfms = {}
    intermediate = None  # The intermediate space when concatenating xfms - includes cohort
    intermediate_targets = (
        {
            'MNI152NLin6Asym',  # TODO: 'MNI152NLin2009cAsym'
        }
        if config.workflow.multi_step_reg
        else set()
    )

    for template in spaces.get_spaces(nonstandard=False, dim=(3,)):
        # resolution / spec will not differentiate here
        if template.startswith('MNIInfant'):
            intermediate = template
        xfms = precomputed.get('transforms', {}).get(template, {})
        if set(xfms) != {'forward', 'reverse'}:
            if template in intermediate_targets:
                concat_xfms.append(template)
            else:
                templates.append(template)
        else:
            found_xfms[template] = xfms

    # Create another set of buffers to handle the case where we aggregate found and generated
    # xfms to be concatenated
    concat_template_buffer = pe.Node(niu.Merge(2), name='concat_template_buffer')
    concat_template_buffer.inputs.in1 = list(found_xfms)
    concat_anat2std_buffer = pe.Node(niu.Merge(2), name='concat_anat2std_buffer')
    concat_anat2std_buffer.inputs.in1 = [xfm['forward'] for xfm in found_xfms.values()]
    concat_std2anat_buffer = pe.Node(niu.Merge(2), name='concat_std2anat_buffer')
    concat_std2anat_buffer.inputs.in1 = [xfm['reverse'] for xfm in found_xfms.values()]

    template_buffer.inputs.in1 = list(found_xfms)
    anat2std_buffer.inputs.in1 = [xfm['forward'] for xfm in found_xfms.values()]
    std2anat_buffer.inputs.in1 = [xfm['reverse'] for xfm in found_xfms.values()]

    if config.workflow.norm_csf:
        csf_norm_wf = init_csf_norm_wf()

        workflow.connect([
            (anat_buffer, csf_norm_wf, [('anat_preproc', 'inputnode.anat_preproc')]),
            (seg_buffer, csf_norm_wf, [('anat_tpms', 'inputnode.anat_tpms')]),
            (csf_norm_wf, anat_reg_buffer, [('outputnode.anat_preproc', 'anat_preproc')]),
        ])  # fmt:skip

    if templates:
        LOGGER.info(f'ANAT Stage 5a: Preparing normalization workflow for {templates}')
        register_template_wf = init_register_template_wf(
            sloppy=sloppy,
            omp_nthreads=omp_nthreads,
            templates=templates,
            image_type=reference_anat,
        )
        ds_template_registration_wf = init_ds_template_registration_wf(
            output_dir=str(output_dir),
            image_type=reference_anat,
        )

        workflow.connect([
            (inputnode, register_template_wf, [('roi', 'inputnode.lesion_mask')]),
            (anat_reg_buffer, register_template_wf, [
                ('anat_preproc', 'inputnode.moving_image'),
            ]),
            (refined_buffer, register_template_wf, [('anat_mask', 'inputnode.moving_mask')]),
            (sourcefile_buffer, ds_template_registration_wf, [
                ('anat_source_files', 'inputnode.source_files')
            ]),
            (register_template_wf, ds_template_registration_wf, [
                ('outputnode.template', 'inputnode.template'),
                ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
            ]),
            (register_template_wf, template_buffer, [('outputnode.template', 'in2')]),
            (register_template_wf, concat_template_buffer, [('outputnode.template', 'in2')]),
            (register_template_wf, std2anat_buffer, [('outputnode.std2anat_xfm', 'in2')]),
            (register_template_wf, concat_std2anat_buffer, [('outputnode.std2anat_xfm', 'in2')]),
            (register_template_wf, anat2std_buffer, [('outputnode.anat2std_xfm', 'in2')]),
            (register_template_wf, concat_anat2std_buffer, [('outputnode.anat2std_xfm', 'in2')]),
        ])  # fmt:skip

    if concat_xfms and intermediate is None:
        LOGGER.error(
            'Intermediate is not set - skipping concatenation workflow. Spaces: %s',
            ' Intermediate targets: %s',
            spaces.get_spaces(nonstandard=False, dim=(3,)),
            intermediate_targets,
        )
    elif concat_xfms:
        LOGGER.info(f'ANAT Stage 5b: Concatenating normalization for {concat_xfms}')

        select_infant_mni = pe.Node(
            KeySelect(fields=['anat2std_xfm', 'std2anat_xfm'], key=intermediate),
            name='select_infant_mni',
            run_without_submitting=True,
        )
        concat_reg_wf = init_concat_registrations_wf(templates=concat_xfms)
        ds_concat_reg_wf = init_ds_template_registration_wf(
            output_dir=str(output_dir),
            image_type=reference_anat,
            name='ds_concat_registration_wf',
        )

        workflow.connect([
            (concat_template_buffer, select_infant_mni, [('out', 'keys')]),
            (concat_anat2std_buffer, select_infant_mni, [('out', 'anat2std_xfm')]),
            (concat_std2anat_buffer, select_infant_mni, [('out', 'std2anat_xfm')]),
            (select_infant_mni, concat_reg_wf, [
                ('key', 'inputnode.intermediate'),
                ('anat2std_xfm', 'inputnode.anat2int_xfm'),
                ('std2anat_xfm', 'inputnode.int2anat_xfm'),
            ]),
            (sourcefile_buffer, ds_concat_reg_wf, [
                ('anat_source_files', 'inputnode.source_files')
            ]),
            (concat_reg_wf, ds_concat_reg_wf, [
                ('outputnode.template', 'inputnode.template'),
                ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
            ]),
            (concat_reg_wf, template_buffer, [('outputnode.template', 'in3')]),
            (concat_reg_wf, anat2std_buffer, [('outputnode.anat2std_xfm', 'in3')]),
            (concat_reg_wf, std2anat_buffer, [('outputnode.std2anat_xfm', 'in3')]),
        ])  # fmt:skip

    if found_xfms:
        LOGGER.info(f'ANAT Stage 5c: Found pre-computed registrations for {found_xfms}')

    # Only refine mask if necessary
    if anat_mask or recon_method is None or not refine_mask:
        workflow.connect([
            (anat_buffer, refined_buffer, [
                ('anat_mask', 'anat_mask'),
                ('anat_brain', 'anat_brain'),
            ]),
        ])  # fmt:skip

    workflow.__desc__ = desc

    if recon_method is None:
        LOGGER.info('ANAT Skipping Stages 6+')
        return workflow

    # Stage 6: Surface reconstruction
    if recon_method == 'mcribs':
        from nibabies.workflows.anatomical.surfaces import init_mcribs_surface_recon_wf

        LOGGER.info('ANAT Stage 6: Preparing M-CRIB-S reconstruction workflow')
        surface_recon_wf = init_mcribs_surface_recon_wf(
            omp_nthreads=omp_nthreads,
            use_aseg=bool(anat_aseg),
            precomputed=precomputed,
            mcribs_dir=str(config.execution.mcribs_dir),
        )

        workflow.connect([
            (inputnode, surface_recon_wf, [
                ('subject_id', 'inputnode.subject_id'),
                ('subjects_dir', 'inputnode.subjects_dir'),
            ]),
            (t2w_validate, surface_recon_wf, [('out_file', 'inputnode.t2w')]),
            (t2w_buffer, surface_recon_wf, [('t2w_mask', 'inputnode.in_mask'),]),
            (aseg_buffer, surface_recon_wf, [
                ('anat_aseg', 'inputnode.in_aseg'),
            ]),
            (surface_recon_wf, outputnode, [
                ('outputnode.subjects_dir', 'subjects_dir'),
                ('outputnode.subject_id', 'subject_id'),
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
                hires=True,
                fs_no_resume=False,
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
                omp_nthreads=omp_nthreads,
                use_aseg=bool(anat_aseg),
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

        if anat_aseg:
            workflow.connect(aseg_buffer, 'anat_aseg', surface_recon_wf, 'inputnode.in_aseg')

    fsnative_xfms = precomputed.get('transforms', {}).get('fsnative')
    if not fsnative_xfms:
        ds_fs_registration_wf = init_ds_fs_registration_wf(
            image_type=reference_anat, output_dir=output_dir
        )

        if recon_method == 'freesurfer':
            workflow.connect([
                (surface_recon_wf, fsnative_buffer, [
                    ('outputnode.fsnative2t1w_xfm', 'fsnative2anat_xfm'),
                    ('outputnode.t1w2fsnative_xfm', 'anat2fsnative_xfm'),
                ]),
            ])  # fmt:skip
        else:
            workflow.connect([
                (surface_recon_wf, fsnative_buffer, [
                    ('outputnode.fsnative2anat_xfm', 'fsnative2anat_xfm'),
                    ('outputnode.anat2fsnative_xfm', 'anat2fsnative_xfm'),
                ]),
            ])  # fmt:skip

        workflow.connect([
            (sourcefile_buffer, ds_fs_registration_wf, [
                ('anat_source_files', 'inputnode.source_files'),
            ]),
            (fsnative_buffer, ds_fs_registration_wf, [
                ('fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (fsnative_buffer, outputnode, [
                ('fsnative2anat_xfm', 'fsnative2anat_xfm'),
            ]),
        ])  # fmt:skip
    elif 'reverse' in fsnative_xfms:
        LOGGER.info('ANAT Found fsnative-to-anatomical transform - skipping registration')
        outputnode.inputs.fsnative2anat_xfm = fsnative_xfms['reverse']
    else:
        raise RuntimeError(
            'Found an anatomical-to-fsnative transform without the reverse. Time to handle this.'
        )

    if not anat_mask and refine_mask:
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
                ('fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (anat_buffer, refinement_wf, [
                ('anat_preproc', 'inputnode.reference_image'),
            ]),
            (seg_buffer, refinement_wf, [
                ('ants_segs', 'inputnode.ants_segs'),  # TODO: Verify this is the same as dseg
            ]),
            (anat_buffer, applyrefined, [('anat_preproc', 'in_file')]),
            (refinement_wf, applyrefined, [('outputnode.out_brainmask', 'in_mask')]),
            (refinement_wf, refined_buffer, [('outputnode.out_brainmask', 'anat_mask')]),
            (applyrefined, refined_buffer, [('out_file', 'anat_brain')]),
        ])  # fmt:skip
    elif not refine_mask:
        LOGGER.info('ANAT Skipping mask refinement workflow')
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
        ds_surfaces_wf = init_ds_surfaces_wf(output_dir=output_dir, surfaces=surfs)

        workflow.connect([
            (surface_recon_wf, gifti_surfaces_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
            ]),
            (fsnative_buffer, gifti_surfaces_wf, [
                ('fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (sourcefile_buffer, ds_surfaces_wf, [('anat_source_files', 'inputnode.source_files')]),
            (gifti_surfaces_wf, ds_surfaces_wf, [
                (f'outputnode.{surf}', f'inputnode.{surf}') for surf in surfs
            ]),
            (ds_surfaces_wf, surfaces_buffer, [
                (f'outputnode.{surf}', surf) for surf in surfs
            ]),
        ])  # fmt:skip
    if spheres:
        gifti_spheres_wf = init_gifti_surfaces_wf(
            surfaces=spheres, to_scanner=False, name='gifti_spheres_wf'
        )
        ds_spheres_wf = init_ds_surfaces_wf(
            output_dir=output_dir,
            surfaces=spheres,
            name='ds_spheres_wf',
        )

        workflow.connect([
            (surface_recon_wf, gifti_spheres_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                # No transform for spheres, following HCP pipelines' lead
            ]),
            (sourcefile_buffer, ds_spheres_wf, [('anat_source_files', 'inputnode.source_files')]),
            (gifti_spheres_wf, ds_spheres_wf, [
                (f'outputnode.{sphere}', f'inputnode.{sphere}') for sphere in spheres
            ]),
            (ds_spheres_wf, surfaces_buffer, [
                (f'outputnode.{sphere}', sphere) for sphere in spheres
            ]),
        ])  # fmt:skip
    metrics = [metric for metric in needed_metrics if metric not in found_surfs]
    if metrics:
        LOGGER.info(f'ANAT Stage 8: Creating GIFTI metrics for {metrics}')
        gifti_morph_wf = init_gifti_morphometrics_wf(morphometrics=metrics)
        ds_morph_wf = init_ds_surface_metrics_wf(
            bids_root=bids_root,
            output_dir=output_dir,
            metrics=metrics,
            name='ds_morph_wf',
        )

        workflow.connect([
            (surface_recon_wf, gifti_morph_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
            ]),
            (sourcefile_buffer, ds_morph_wf, [('anat_source_files', 'inputnode.source_files')]),
            (gifti_morph_wf, ds_morph_wf, [
                (f'outputnode.{metric}', f'inputnode.{metric}') for metric in metrics
            ]),
            (ds_morph_wf, surfaces_buffer, [
                (f'outputnode.{metric}', metric) for metric in metrics
            ]),
        ])  # fmt:skip

    if 'anat_ribbon' not in precomputed:
        LOGGER.info('ANAT Stage 8a: Creating cortical ribbon mask')
        anat_ribbon_wf = init_anat_ribbon_wf()
        ds_ribbon_mask_wf = init_ds_mask_wf(
            bids_root=bids_root,
            output_dir=output_dir,
            mask_type='ribbon',
            extra_entities={'space': reference_anat},
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
            (sourcefile_buffer, ds_ribbon_mask_wf, [
                ('anat_source_files', 'inputnode.source_files'),
            ]),
            (anat_ribbon_wf, ds_ribbon_mask_wf, [
                ('outputnode.anat_ribbon', 'inputnode.mask_file'),
            ]),
            (ds_ribbon_mask_wf, outputnode, [('outputnode.mask_file', 'anat_ribbon')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Stage 8a: Found pre-computed cortical ribbon mask')
        outputnode.inputs.anat_ribbon = precomputed['anat_ribbon']

    # Stage 9: Baseline fsLR registration
    if recon_method == 'mcribs':
        if len(precomputed.get('sphere_reg_dhcpAsym', [])) < 2:
            LOGGER.info('ANAT Stage 9: Creating dhcp-fsLR registration sphere')
            fsLR_reg_wf = init_mcribs_dhcp_wf()

            ds_fsLR_reg_wf = init_ds_surfaces_wf(
                output_dir=output_dir,
                surfaces=['sphere_reg_dhcpAsym'],
                name='ds_fsLR_reg_wf',
            )

            workflow.connect([
                (surfaces_buffer, fsLR_reg_wf, [('sphere_reg', 'inputnode.sphere_reg')]),
                (sourcefile_buffer, ds_fsLR_reg_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (fsLR_reg_wf, ds_fsLR_reg_wf, [
                    ('outputnode.sphere_reg_dhcpAsym', 'inputnode.sphere_reg_dhcpAsym')
                ]),
                (ds_fsLR_reg_wf, fsLR_buffer, [
                    ('outputnode.sphere_reg_dhcpAsym', 'sphere_reg_fsLR'),
                ]),
            ])  # fmt:skip
        else:
            LOGGER.info('ANAT Stage 9: Found pre-computed dhcp-fsLR registration sphere')
            fsLR_buffer.inputs.sphere_reg_fsLR = sorted(precomputed['sphere_reg_dhcpAsym'])

    else:
        if len(precomputed.get('sphere_reg_fsLR', [])) < 2:
            LOGGER.info('ANAT Stage 9: Creating fsLR registration sphere')
            fsLR_reg_wf = init_fsLR_reg_wf()

            ds_fsLR_reg_wf = init_ds_surfaces_wf(
                output_dir=output_dir,
                surfaces=['sphere_reg_fsLR'],
                name='ds_fsLR_reg_wf',
            )

            workflow.connect([
                (surfaces_buffer, fsLR_reg_wf, [('sphere_reg', 'inputnode.sphere_reg')]),
                (sourcefile_buffer, ds_fsLR_reg_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (fsLR_reg_wf, ds_fsLR_reg_wf, [
                    ('outputnode.sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR')
                ]),
                (ds_fsLR_reg_wf, fsLR_buffer, [('outputnode.sphere_reg_fsLR', 'sphere_reg_fsLR')]),
            ])  # fmt:skip
        else:
            LOGGER.info('ANAT Stage 9: Found pre-computed fsLR registration sphere')
            fsLR_buffer.inputs.sphere_reg_fsLR = sorted(precomputed['sphere_reg_fsLR'])

    # Stage 10: Cortical surface mask
    if len(precomputed.get('cortex_mask', [])) < 2:
        LOGGER.info('ANAT Stage 11: Creating cortical surface mask')

        cortex_masks_wf = init_cortex_masks_wf()
        ds_cortex_masks_wf = init_ds_surface_masks_wf(
            output_dir=output_dir,
            mask_type='cortex',
            name='ds_cortex_masks_wf',
        )

        workflow.connect([
            (surfaces_buffer, cortex_masks_wf, [
                ('midthickness', 'inputnode.midthickness'),
                ('thickness', 'inputnode.thickness'),
            ]),
            (cortex_masks_wf, ds_cortex_masks_wf, [
                ('outputnode.cortex_masks', 'inputnode.mask_files'),
                ('outputnode.source_files', 'inputnode.source_files'),
            ]),
            (ds_cortex_masks_wf, outputnode, [('outputnode.mask_files', 'cortex_mask')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Stage 11: Found pre-computed cortical surface mask')
        outputnode.inputs.cortex_mask = sorted(precomputed['cortex_mask'])
    return workflow


@tag('anat.fit')
def init_infant_single_anat_fit_wf(
    *,
    age_months: int,
    t1w: list,
    t2w: list,
    flair: list,
    bids_root: str,
    precomputed: dict,
    longitudinal: bool,
    omp_nthreads: int,
    output_dir: str,
    segmentation_atlases: str | Path | None,
    skull_strip_mode: ty.Literal['force', 'skip', 'auto'],
    skull_strip_template: 'Reference',
    skull_strip_fixed_seed: bool,
    sloppy: bool,
    spaces: 'SpatialReferences',
    recon_method: ty.Literal['freesurfer', 'infantfs', 'mcribs'] | None,
    reference_anat: ty.Literal['T1w', 'T2w'],
    cifti_output: ty.Literal['91k', '170k', False],
    msm_sulc: bool = False,
    name: str = 'infant_single_anat_fit_wf',
) -> Workflow:
    """
    Create a fit workflow with just a single anatomical.

    Note: Treat this functionality as a workaround that enables the processing of incomplete data.
    For best results, especially in periods of transitioning myelination (usually 3-8 months),
    a combination of T1w and T2w images will produce more accurate results.
    """
    anatomicals = t1w or t2w

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['t1w', 't2w', 'anat', 'roi', 'flair', 'subjects_dir', 'subject_id'],
        ),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                # Primary derivatives
                'anat_preproc',
                'anat_mask',
                'anat_dseg',
                'anat_tpms',
                'anat2std_xfm',
                'fsnative2anat_xfm',
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
                'cortex_mask',
                'anat_ribbon',
                # Reverse transform; not computable from forward transform
                'std2anat_xfm',
                # Metadata
                'template',
                'subjects_dir',
                'subject_id',
                'anat_valid_list',
            ]
        ),
        name='outputnode',
    )

    anat = reference_anat.lower()
    workflow = Workflow(name=f'infant_single_{anat}_fit_wf')
    workflow.add_nodes([inputnode, outputnode])

    # Stage 1 inputs (filtered)
    sourcefile_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_source_files']),
        name='sourcefile_buffer',
    )

    # Stage 2
    anat_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_preproc', 'anat_mask', 'anat_brain']),
        name='anat_buffer',
    )

    # Additional buffer if CSF normalization is used
    anat_reg_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_preproc']),
        name='anat_reg_buffer',
    )
    if not config.workflow.norm_csf:
        workflow.connect(anat_buffer, 'anat_preproc', anat_reg_buffer, 'anat_preproc')

    aseg_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_aseg']),
        name='aseg_buffer',
    )

    # Stage 3 - Segmentation
    seg_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_dseg', 'anat_tpms', 'ants_segs']),
        name='seg_buffer',
    )
    # Stage 4 - collated template names, forward and reverse transforms
    template_buffer = pe.Node(niu.Merge(3), name='template_buffer')
    anat2std_buffer = pe.Node(niu.Merge(3), name='anat2std_buffer')
    std2anat_buffer = pe.Node(niu.Merge(3), name='std2anat_buffer')

    # Stage 5 results: Refined stage 2 results; may be direct copy if no refinement
    refined_buffer = pe.Node(
        niu.IdentityInterface(fields=['anat_mask', 'anat_brain']),
        name='refined_buffer',
    )

    fsnative_buffer = pe.Node(
        niu.IdentityInterface(fields=['fsnative2anat_xfm', 'anat2fsnative_xfm']),
        name='fsnative_buffer',
    )

    # Stage 6 results: GIFTI surfaces
    surfaces_buffer = pe.Node(
        niu.IdentityInterface(
            fields=['white', 'pial', 'midthickness', 'sphere', 'sphere_reg', 'thickness', 'sulc']
        ),
        name='surfaces_buffer',
    )

    # Stage 7 and 8 results: fsLR sphere registration
    fsLR_buffer = pe.Node(niu.IdentityInterface(fields=['sphere_reg_fsLR']), name='fsLR_buffer')
    msm_buffer = pe.Node(niu.IdentityInterface(fields=['sphere_reg_msm']), name='msm_buffer')

    workflow.connect([
        (anat_buffer, outputnode, [
            ('anat_preproc', 'anat_preproc'),
        ]),
        (refined_buffer, outputnode, [
            ('anat_mask', 'anat_mask'),
            ('anat_brain', 'anat_brain'),
        ]),
        (seg_buffer, outputnode, [
            ('anat_dseg', 'anat_dseg'),
            ('anat_tpms', 'anat_tpms'),
        ]),
        (anat2std_buffer, outputnode, [('out', 'anat2std_xfm')]),
        (std2anat_buffer, outputnode, [('out', 'std2anat_xfm')]),
        (template_buffer, outputnode, [('out', 'template')]),
        (sourcefile_buffer, outputnode, [('anat_source_files', 'anat_valid_list')]),
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
        spaces=spaces,
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
        f'A total of {len(anatomicals)} {anat} images were found '
        'within the input BIDS dataset.\n'
    )
    # Lowercase to match pattern
    anat = reference_anat.lower()

    # Derivatives
    anat_preproc = precomputed.get(f'{anat}_preproc')
    anat_mask = precomputed.get(f'{anat}_mask')
    anat_tpms = precomputed.get(f'{anat}_tpms')
    anat_dseg = precomputed.get(f'{anat}_dseg')
    anat_aseg = precomputed.get(f'{anat}_aseg')

    anat_validate = pe.Node(ValidateImage(), name='anat_validate', run_without_submitting=True)
    if not anat_preproc:
        LOGGER.info(f'ANAT Stage 1: Adding {reference_anat} template workflow')
        desc += (
            f'The {reference_anat} image was denoised and corrected for intensity '
            'non-uniformity (INU)'
        )

        anat_template_wf = init_anat_template_wf(
            image_type=reference_anat,
            num_files=len(anatomicals),
            longitudinal=longitudinal,
            omp_nthreads=omp_nthreads,
            name=f'{anat}_template_wf',
        )
        ds_anat_template_wf = init_ds_template_wf(
            image_type=reference_anat,
            output_dir=output_dir,
            num_anat=len(anatomicals),
            name='ds_anat_template_wf',
        )

        workflow.connect([
            (anat_template_wf, sourcefile_buffer, [
                ('outputnode.anat_valid_list', 'anat_source_files'),
            ]),
            (anat_template_wf, anat_reports_wf, [
                ('outputnode.out_report', 'inputnode.anat_conform_report'),
            ]),
        ])  # fmt:skip

        workflow.connect([
            (inputnode, anat_template_wf, [('anat', 'inputnode.anat_files')]),
            (anat_template_wf, anat_validate, [('outputnode.anat_ref', 'in_file')]),
            (anat_template_wf, ds_anat_template_wf, [
                ('outputnode.anat_realign_xfm', 'inputnode.anat_ref_xfms'),
            ]),
            (sourcefile_buffer, ds_anat_template_wf, [
                ('anat_source_files', 'inputnode.source_files'),
            ]),
            (anat_buffer, ds_anat_template_wf, [('anat_preproc', 'inputnode.anat_preproc')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Found preprocessed T1w - skipping Stage 1')
        desc += ' A preprocessed T1w image was provided as input.'

        anat_validate.inputs.in_file = anat_preproc
        sourcefile_buffer.inputs.anat_source_files = [anat_preproc]

        workflow.connect([
            (anat_validate, anat_buffer, [('out_file', 'anat_preproc')]),
        ])  # fmt:skip

    # Stage 2: Use previously computed mask or calculate
    # If we only have one mask (could be either T1w/T2w),
    # just apply transform to get it in the other space
    # T2w masking logic:
    # PCM = Pre-computed mask
    # PCM = Pre-computed mask for PCT
    # PCT = Pre-computed template
    # SS = Skull stripping required
    #
    # PCM, PCT, SS -> Apply PCM to PCT
    # PCM, PCT, !SS -> Apply PCM to PCT
    # PCM, !PCT, SS -> Apply PCM to template, then run INU
    # PCM, !PCT, !SS -> Apply PCM to template, then run INU
    # !PCM, PCT, SS -> Run brain extraction
    # !PCM, PCT, !SS -> Binarize PCT
    # !PCM, !PCT, SS -> Run brain extraction
    # !PCM, !PCT, !SS -> INU correct template
    refine_mask = False
    apply_mask = pe.Node(ApplyMask(), name='apply_mask')
    if not anat_mask:
        if skull_strip_mode == 'auto':
            run_skull_strip = not all(_is_skull_stripped(img) for img in anatomicals)
        else:
            run_skull_strip = {'force': True, 'skip': False}[skull_strip_mode]

        if not run_skull_strip:
            desc += (
                f'The {reference_anat} reference was previously skull-stripped; '
                'a brain mask was derived from the input image.'
            )

            if not anat_preproc:
                LOGGER.info('ANAT Stage 2b: Skipping skull-strip, INU-correction only')

                n4_only_wf = init_n4_only_wf(
                    omp_nthreads=omp_nthreads,
                    bids_suffix=reference_anat,
                    atropos_use_random_seed=not skull_strip_fixed_seed,
                    name='n4_only_wf',
                )

                refine_mask = True
                workflow.connect([
                    (anat_validate, n4_only_wf, [('out_file', 'inputnode.in_files')]),
                    (n4_only_wf, anat_buffer, [
                        (('outputnode.bias_corrected', pop_file), 'anat_preproc'),
                        ('outputnode.out_mask', 'anat_mask'),
                        (('outputnode.out_file', pop_file), 'anat_brain'),
                    ]),
                    (n4_only_wf, seg_buffer, [
                        ('outputnode.out_segm', 'ants_segs'),
                    ]),
                ])  # fmt:skip
            else:
                LOGGER.info('ANAT Stage 2b: Skipping skull-strip, generating mask from input')
                binarize = pe.Node(Binarize(thresh_low=2), name='binarize')
                workflow.connect([
                    (anat_validate, binarize, [('out_file', 'in_file')]),
                    (anat_validate, anat_buffer, [('out_file', 'anat_brain')]),
                    (binarize, anat_buffer, [('out_file', 'anat_mask')]),
                ])  # fmt:skip
        else:
            LOGGER.info(f'ANAT Atlas-based brain mask will be calculated on the {reference_anat}')
            brain_extraction_wf = init_infant_brain_extraction_wf(
                omp_nthreads=omp_nthreads,
                sloppy=sloppy,
                age_months=age_months,
                ants_affine_init=True,
                skull_strip_template=skull_strip_template.space,
                template_specs=skull_strip_template.spec,
                debug='registration' in config.execution.debug,
            )

            workflow.connect([
                (anat_validate, brain_extraction_wf, [
                    ('out_file', 'inputnode.t2w_preproc'),
                ]),
                (brain_extraction_wf, anat_buffer, [
                    ('outputnode.out_mask', 'anat_mask'),
                    ('outputnode.t2w_brain', 'anat_brain'),
                ]),
            ])  # fmt:skip

    else:
        LOGGER.info(f'ANAT Found {reference_anat} brain mask')
        desc += 'A pre-computed brain mask was provided as input and used throughout the workflow.'
        workflow.connect(apply_mask, 'out_file', anat_buffer, 'anat_brain')

        if not anat_preproc:
            conform_anat_mask_wf = init_conform_derivative_wf(
                in_file=anat_mask,
                name='conform_anat_mask_wf',
            )

            LOGGER.info('ANAT Skipping skull-strip, INU-correction only')
            anat_n4_wf = init_anat_preproc_wf(name='anat_n4_wf')
            workflow.connect([
                (anat_validate, conform_anat_mask_wf, [('out_file', 'inputnode.ref_file')]),
                (conform_anat_mask_wf, anat_buffer, [('outputnode.out_file', 'anat_mask')]),
                (conform_anat_mask_wf, apply_mask, [('outputnode.out_file', 'in_mask')]),
                (anat_validate, anat_n4_wf, [('out_file', 'inputnode.in_anat')]),
                (anat_n4_wf, anat_buffer, [('outputnode.anat_preproc', 'anat_preproc')]),
                (anat_n4_wf, apply_mask, [('outputnode.anat_preproc', 'in_file')]),
            ])  # fmt:skip
        else:
            LOGGER.info(f'ANAT Skipping {reference_anat} masking')
            workflow.connect(anat_validate, 'out_file', apply_mask, 'in_file')
            anat_buffer.inputs.anat_mask = anat_mask
            apply_mask.inputs.in_mask = anat_mask

    # Stage 3: Segmentation
    seg_method = 'jlf' if config.execution.segmentation_atlases_dir else 'fast'
    if anat_aseg:
        LOGGER.info('ANAT Found precomputed anatomical segmentation')
        # Ensure compatibility with anatomical template
        if not anat_preproc:
            conform_aseg_wf = init_conform_derivative_wf(
                in_file=anat_aseg,
                name='conform_aseg_wf',
            )

            workflow.connect([
                (anat_buffer, conform_aseg_wf, [('anat_preproc', 'inputnode.ref_file')]),
                (conform_aseg_wf, aseg_buffer, [('outputnode.out_file', 'anat_aseg')]),
            ])  # fmt:skip
        else:
            aseg_buffer.inputs.anat_aseg = anat_aseg

    if not (anat_dseg and anat_tpms):
        LOGGER.info('ANAT Stage 3: Tissue segmentation')
        segmentation_wf = init_segmentation_wf(
            sloppy=sloppy,
            method=seg_method,
            image_type=reference_anat,
            omp_nthreads=omp_nthreads,
            has_aseg=bool(anat_aseg),
        )

        workflow.connect([
            (anat_buffer, segmentation_wf, [('anat_brain', 'inputnode.anat_brain')]),
            (segmentation_wf, seg_buffer, [
                ('outputnode.anat_dseg', 'anat_dseg'),
                ('outputnode.anat_tpms', 'anat_tpms'),
            ]),
        ])  # fmt:skip

        if anat_aseg:
            workflow.connect(aseg_buffer, 'anat_aseg', segmentation_wf, 'inputnode.anat_aseg')
        elif seg_method == 'jlf':
            workflow.connect(segmentation_wf, 'outputnode.anat_aseg', aseg_buffer, 'anat_aseg')
            # TODO: datasink aseg

        if not anat_dseg:
            ds_dseg_wf = init_ds_dseg_wf(
                output_dir=str(output_dir),
                extra_entities={'space': reference_anat},
            )
            workflow.connect([
                (sourcefile_buffer, ds_dseg_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (segmentation_wf, ds_dseg_wf, [
                    ('outputnode.anat_dseg', 'inputnode.anat_dseg'),
                ]),
            ])  # fmt:skip

        if not anat_tpms:
            ds_tpms_wf = init_ds_tpms_wf(
                output_dir=str(output_dir),
                extra_entities={'space': reference_anat},
            )
            workflow.connect([
                (sourcefile_buffer, ds_tpms_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (segmentation_wf, ds_tpms_wf, [
                    ('outputnode.anat_tpms', 'inputnode.anat_tpms'),
                ]),
            ])  # fmt:skip
    else:
        LOGGER.info('ANAT Stage 3: Skipping segmentation workflow')
    if anat_dseg:
        LOGGER.info('ANAT Found discrete segmentation')
        desc += 'Precomputed discrete tissue segmentations were provided as inputs.\n'
        seg_buffer.inputs.anat_dseg = anat_dseg
    if anat_tpms:
        LOGGER.info('ANAT Found tissue probability maps')
        desc += 'Precomputed tissue probabiilty maps were provided as inputs.\n'
        seg_buffer.inputs.anat_tpms = anat_tpms

    # Stage 4: Normalization
    # If selected adult templates are requested (MNI152 6th Gen or 2009)
    # opt to concatenate transforms first from native -> infant template (MNIInfant),
    # and then use a previously computed MNIInfant<cohort> -> MNI transform
    # this minimizes the chance of a bad registration.

    templates = []
    concat_xfms = []
    found_xfms = {}
    intermediate = None  # The intermediate space when concatenating xfms - includes cohort
    intermediate_targets = (
        {
            'MNI152NLin6Asym',  # TODO: 'MNI152NLin2009cAsym'
        }
        if config.workflow.multi_step_reg
        else set()
    )

    for template in spaces.get_spaces(nonstandard=False, dim=(3,)):
        # resolution / spec will not differentiate here
        if template.startswith('MNIInfant'):
            intermediate = template
        xfms = precomputed.get('transforms', {}).get(template, {})
        if set(xfms) != {'forward', 'reverse'}:
            if template in intermediate_targets:
                concat_xfms.append(template)
            else:
                templates.append(template)
        else:
            found_xfms[template] = xfms

    # Create another set of buffers to handle the case where we aggregate found and generated
    # xfms to be concatenated
    concat_template_buffer = pe.Node(niu.Merge(2), name='concat_template_buffer')
    concat_template_buffer.inputs.in1 = list(found_xfms)
    concat_anat2std_buffer = pe.Node(niu.Merge(2), name='concat_anat2std_buffer')
    concat_anat2std_buffer.inputs.in1 = [xfm['forward'] for xfm in found_xfms.values()]
    concat_std2anat_buffer = pe.Node(niu.Merge(2), name='concat_std2anat_buffer')
    concat_std2anat_buffer.inputs.in1 = [xfm['reverse'] for xfm in found_xfms.values()]

    template_buffer.inputs.in1 = list(found_xfms)
    anat2std_buffer.inputs.in1 = [xfm['forward'] for xfm in found_xfms.values()]
    std2anat_buffer.inputs.in1 = [xfm['reverse'] for xfm in found_xfms.values()]

    if config.workflow.norm_csf:
        csf_norm_wf = init_csf_norm_wf()

        workflow.connect([
            (anat_buffer, csf_norm_wf, [('anat_preproc', 'inputnode.anat_preproc')]),
            (seg_buffer, csf_norm_wf, [('anat_tpms', 'inputnode.anat_tpms')]),
            (csf_norm_wf, anat_reg_buffer, [('outputnode.anat_preproc', 'anat_preproc')]),
        ])  # fmt:skip

    if templates:
        LOGGER.info(f'ANAT Stage 4: Preparing normalization workflow for {templates}')
        register_template_wf = init_register_template_wf(
            sloppy=sloppy,
            omp_nthreads=omp_nthreads,
            templates=templates,
            image_type=reference_anat,
        )
        ds_template_registration_wf = init_ds_template_registration_wf(
            output_dir=str(output_dir),
            image_type=reference_anat,
        )

        workflow.connect([
            (inputnode, register_template_wf, [('roi', 'inputnode.lesion_mask')]),
            (anat_reg_buffer, register_template_wf, [
                ('anat_preproc', 'inputnode.moving_image'),
            ]),
            (refined_buffer, register_template_wf, [('anat_mask', 'inputnode.moving_mask')]),
            (sourcefile_buffer, ds_template_registration_wf, [
                ('anat_source_files', 'inputnode.source_files')
            ]),
            (register_template_wf, ds_template_registration_wf, [
                ('outputnode.template', 'inputnode.template'),
                ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
            ]),
            (register_template_wf, template_buffer, [('outputnode.template', 'in2')]),
            (register_template_wf, concat_template_buffer, [('outputnode.template', 'in2')]),
            (register_template_wf, std2anat_buffer, [('outputnode.std2anat_xfm', 'in2')]),
            (register_template_wf, concat_std2anat_buffer, [('outputnode.std2anat_xfm', 'in2')]),
            (register_template_wf, anat2std_buffer, [('outputnode.anat2std_xfm', 'in2')]),
            (register_template_wf, concat_anat2std_buffer, [('outputnode.anat2std_xfm', 'in2')]),
        ])  # fmt:skip

    if concat_xfms and intermediate is None:
        LOGGER.error(
            'Intermediate is not set - skipping concatenation workflow.\nSpaces: %s\n',
            'Intermediate targets: %s',
            spaces.get_spaces(nonstandard=False, dim=(3,)),
            intermediate_targets,
        )
    elif concat_xfms:
        LOGGER.info(f'ANAT Stage 5b: Concatenating normalization for {concat_xfms}')

        select_infant_mni = pe.Node(
            KeySelect(fields=['anat2std_xfm', 'std2anat_xfm'], key=intermediate),
            name='select_infant_mni',
            run_without_submitting=True,
        )

        concat_reg_wf = init_concat_registrations_wf(templates=concat_xfms)
        ds_concat_reg_wf = init_ds_template_registration_wf(
            output_dir=str(output_dir),
            image_type=reference_anat,
            name='ds_concat_registration_wf',
        )

        workflow.connect([
            (concat_template_buffer, select_infant_mni, [('out', 'keys')]),
            (concat_anat2std_buffer, select_infant_mni, [('out', 'anat2std_xfm')]),
            (concat_std2anat_buffer, select_infant_mni, [('out', 'std2anat_xfm')]),
            (select_infant_mni, concat_reg_wf, [
                ('key', 'inputnode.intermediate'),
                ('anat2std_xfm', 'inputnode.anat2int_xfm'),
                ('std2anat_xfm', 'inputnode.int2anat_xfm'),
            ]),
            (sourcefile_buffer, ds_concat_reg_wf, [
                ('anat_source_files', 'inputnode.source_files')
            ]),
            (concat_reg_wf, ds_concat_reg_wf, [
                ('outputnode.template', 'inputnode.template'),
                ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
            ]),
            (concat_reg_wf, template_buffer, [('outputnode.template', 'in3')]),
            (concat_reg_wf, anat2std_buffer, [('outputnode.anat2std_xfm', 'in3')]),
            (concat_reg_wf, std2anat_buffer, [('outputnode.std2anat_xfm', 'in3')]),
        ])  # fmt:skip

    if found_xfms:
        LOGGER.info(f'ANAT Stage 4: Found pre-computed registrations for {found_xfms}')

    # Only refine mask if necessary
    if anat_mask or recon_method is None:
        workflow.connect([
            (anat_buffer, refined_buffer, [
                ('anat_mask', 'anat_mask'),
                ('anat_brain', 'anat_brain'),
            ]),
        ])  # fmt:skip

    workflow.__desc__ = desc

    if recon_method is None:
        LOGGER.info('ANAT Skipping Stages 5+')
        return workflow

    # Stage 5: Surface reconstruction
    if recon_method == 'mcribs':
        if reference_anat == 'T1w':
            LOGGER.warning('Attempting to use MCRIBS with a T1w file, good luck.')

        from nibabies.workflows.anatomical.surfaces import init_mcribs_surface_recon_wf

        LOGGER.info('ANAT Stage 5: Preparing M-CRIB-S reconstruction workflow')
        surface_recon_wf = init_mcribs_surface_recon_wf(
            omp_nthreads=omp_nthreads,
            use_aseg=bool(anat_aseg),
            precomputed=precomputed,
            mcribs_dir=str(config.execution.mcribs_dir),
        )

        workflow.connect([
            (inputnode, surface_recon_wf, [
                ('subject_id', 'inputnode.subject_id'),
                ('subjects_dir', 'inputnode.subjects_dir'),
            ]),
            (anat_validate, surface_recon_wf, [('out_file', 'inputnode.t2w')]),
            (anat_buffer, surface_recon_wf, [('anat_mask', 'inputnode.in_mask')]),
            (aseg_buffer, surface_recon_wf, [
                ('anat_aseg', 'inputnode.in_aseg'),
            ]),
            (surface_recon_wf, outputnode, [
                ('outputnode.subjects_dir', 'subjects_dir'),
                ('outputnode.subject_id', 'subject_id'),
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

            LOGGER.info('ANAT Stage 5: Preparing FreeSurfer recon-all workflow')
            fs_isrunning = pe.Node(
                niu.Function(function=fs_isRunning), overwrite=True, name='fs_isrunning'
            )
            fs_isrunning.inputs.logger = LOGGER

            surface_recon_wf = init_surface_recon_wf(
                name='surface_recon_wf',
                omp_nthreads=omp_nthreads,
                hires=True,
                fs_no_resume=False,
                precomputed=precomputed,
            )

        elif recon_method == 'infantfs':
            from nibabies.workflows.anatomical.surfaces import init_infantfs_surface_recon_wf

            LOGGER.info('ANAT Stage 5: Preparing Infant FreeSurfer workflow')
            surface_recon_wf = init_infantfs_surface_recon_wf(
                age_months=age_months,
                precomputed=precomputed,
                omp_nthreads=omp_nthreads,
                use_aseg=bool(anat_aseg),
            )

        workflow.connect([
            (inputnode, fs_isrunning, [
                ('subjects_dir', 'subjects_dir'),
                ('subject_id', 'subject_id'),
            ]),
            (inputnode, surface_recon_wf, [
                ('subject_id', 'inputnode.subject_id'),
            ]),
            (fs_isrunning, surface_recon_wf, [('out', 'inputnode.subjects_dir')]),
            (anat_validate, surface_recon_wf, [('out_file', 'inputnode.t1w')]),
            (anat_buffer, surface_recon_wf, [('anat_brain', 'inputnode.skullstripped_t1')]),
            (surface_recon_wf, outputnode, [
                ('outputnode.subjects_dir', 'subjects_dir'),
                ('outputnode.subject_id', 'subject_id'),
            ]),
        ])  # fmt:skip

        if anat_aseg:
            workflow.connect(aseg_buffer, 'anat_aseg', surface_recon_wf, 'inputnode.in_aseg')

    fsnative_xfms = precomputed.get('transforms', {}).get('fsnative')
    if not fsnative_xfms:
        ds_fs_registration_wf = init_ds_fs_registration_wf(
            image_type=reference_anat, output_dir=output_dir
        )

        if recon_method == 'freesurfer':
            workflow.connect([
                (surface_recon_wf, fsnative_buffer, [
                    ('outputnode.fsnative2t1w_xfm', 'fsnative2anat_xfm'),
                    ('outputnode.t1w2fsnative_xfm', 'anat2fsnative_xfm'),
                ]),
            ])  # fmt:skip
        else:
            workflow.connect([
                (surface_recon_wf, fsnative_buffer, [
                    ('outputnode.fsnative2anat_xfm', 'fsnative2anat_xfm'),
                    ('outputnode.anat2fsnative_xfm', 'anat2fsnative_xfm'),
                ]),
            ])  # fmt:skip

        workflow.connect([
            (sourcefile_buffer, ds_fs_registration_wf, [
                ('anat_source_files', 'inputnode.source_files'),
            ]),
            (fsnative_buffer, ds_fs_registration_wf, [
                ('fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (fsnative_buffer, outputnode, [
                ('fsnative2anat_xfm', 'fsnative2anat_xfm'),
            ]),
        ])  # fmt:skip
    elif 'reverse' in fsnative_xfms:
        LOGGER.info('ANAT Found fsnative-to-anatomical transform - skipping registration')
        outputnode.inputs.fsnative2anat_xfm = fsnative_xfms['reverse']
    else:
        raise RuntimeError(
            'Found an anatomical-to-fsnative transform without the reverse. Time to handle this.'
        )

    if not anat_mask and refine_mask:
        LOGGER.info('ANAT Stage 6: Preparing mask refinement workflow')
        # Stage 6: Refine ANTs mask with FreeSurfer segmentation
        refinement_wf = init_refinement_wf()
        applyrefined = pe.Node(ApplyMask(), name='applyrefined')

        workflow.connect([
            (surface_recon_wf, refinement_wf, [
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                ('outputnode.subject_id', 'inputnode.subject_id'),
            ]),
            (fsnative_buffer, refinement_wf, [
                ('fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (anat_buffer, refinement_wf, [
                ('anat_preproc', 'inputnode.reference_image'),
            ]),
            (seg_buffer, refinement_wf, [
                ('ants_segs', 'inputnode.ants_segs'),
            ]),
            (anat_buffer, applyrefined, [('anat_preproc', 'in_file')]),
            (refinement_wf, applyrefined, [('outputnode.out_brainmask', 'in_mask')]),
            (refinement_wf, refined_buffer, [('outputnode.out_brainmask', 'anat_mask')]),
            (applyrefined, refined_buffer, [('out_file', 'anat_brain')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Found brain mask - skipping Stage 6')

    # Stages 7-9: Surface conversion and registration
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
        LOGGER.info(f'ANAT Stage 7: Found pre-converted surfaces for {list(found_surfs)}')
        surfaces_buffer.inputs.trait_set(**found_surfs)

    # Stage 8: Surface conversion
    surfs = [surf for surf in needed_anat_surfs if surf not in found_surfs]
    spheres = [sphere for sphere in needed_spheres if sphere not in found_surfs]
    if surfs or spheres:
        LOGGER.info(f'ANAT Stage 7: Creating GIFTI surfaces for {surfs + spheres}')
    if surfs:
        gifti_surfaces_wf = init_gifti_surfaces_wf(surfaces=surfs)
        ds_surfaces_wf = init_ds_surfaces_wf(output_dir=output_dir, surfaces=surfs)

        workflow.connect([
            (surface_recon_wf, gifti_surfaces_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
            ]),
            (fsnative_buffer, gifti_surfaces_wf, [
                ('fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (sourcefile_buffer, ds_surfaces_wf, [('anat_source_files', 'inputnode.source_files')]),
            (gifti_surfaces_wf, ds_surfaces_wf, [
                (f'outputnode.{surf}', f'inputnode.{surf}') for surf in surfs
            ]),
            (ds_surfaces_wf, surfaces_buffer, [
                (f'outputnode.{surf}', surf) for surf in surfs
            ]),
        ])  # fmt:skip
    if spheres:
        gifti_spheres_wf = init_gifti_surfaces_wf(
            surfaces=spheres, to_scanner=False, name='gifti_spheres_wf'
        )
        ds_spheres_wf = init_ds_surfaces_wf(
            output_dir=output_dir,
            surfaces=spheres,
            name='ds_spheres_wf',
        )

        workflow.connect([
            (surface_recon_wf, gifti_spheres_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                # No transform for spheres, following HCP pipelines' lead
            ]),
            (sourcefile_buffer, ds_spheres_wf, [('anat_source_files', 'inputnode.source_files')]),
            (gifti_spheres_wf, ds_spheres_wf, [
                (f'outputnode.{sphere}', f'inputnode.{sphere}') for sphere in spheres
            ]),
            (ds_spheres_wf, surfaces_buffer, [
                (f'outputnode.{sphere}', sphere) for sphere in spheres
            ]),
        ])  # fmt:skip
    metrics = [metric for metric in needed_metrics if metric not in found_surfs]
    if metrics:
        LOGGER.info(f'ANAT Stage 8: Creating GIFTI metrics for {metrics}')
        gifti_morph_wf = init_gifti_morphometrics_wf(morphometrics=metrics)
        ds_morph_wf = init_ds_surface_metrics_wf(
            bids_root=bids_root,
            output_dir=output_dir,
            metrics=metrics,
            name='ds_morph_wf',
        )

        workflow.connect([
            (surface_recon_wf, gifti_morph_wf, [
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
            ]),
            (sourcefile_buffer, ds_morph_wf, [('anat_source_files', 'inputnode.source_files')]),
            (gifti_morph_wf, ds_morph_wf, [
                (f'outputnode.{metric}', f'inputnode.{metric}') for metric in metrics
            ]),
            (ds_morph_wf, surfaces_buffer, [
                (f'outputnode.{metric}', metric) for metric in metrics
            ]),
        ])  # fmt:skip

    if 'anat_ribbon' not in precomputed:
        LOGGER.info('ANAT Stage 8a: Creating cortical ribbon mask')
        anat_ribbon_wf = init_anat_ribbon_wf()
        ds_ribbon_mask_wf = init_ds_mask_wf(
            bids_root=bids_root,
            output_dir=output_dir,
            mask_type='ribbon',
            extra_entities={'space': reference_anat},
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
            (sourcefile_buffer, ds_ribbon_mask_wf, [
                ('anat_source_files', 'inputnode.source_files'),
            ]),
            (anat_ribbon_wf, ds_ribbon_mask_wf, [
                ('outputnode.anat_ribbon', 'inputnode.mask_file'),
            ]),
            (ds_ribbon_mask_wf, outputnode, [('outputnode.mask_file', 'anat_ribbon')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Stage 8a: Found pre-computed cortical ribbon mask')
        outputnode.inputs.anat_ribbon = precomputed['anat_ribbon']

    # Stage 9: Baseline fsLR registration
    if recon_method == 'mcribs':
        if len(precomputed.get('sphere_reg_dhcpAsym', [])) < 2:
            LOGGER.info('ANAT Stage 9: Creating dhcp-fsLR registration sphere')
            fsLR_reg_wf = init_mcribs_dhcp_wf()

            ds_fsLR_reg_wf = init_ds_surfaces_wf(
                output_dir=output_dir,
                surfaces=['sphere_reg_dhcpAsym'],
                name='ds_fsLR_reg_wf',
            )

            workflow.connect([
                (surfaces_buffer, fsLR_reg_wf, [('sphere_reg', 'inputnode.sphere_reg')]),
                (sourcefile_buffer, ds_fsLR_reg_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (fsLR_reg_wf, ds_fsLR_reg_wf, [
                    ('outputnode.sphere_reg_dhcpAsym', 'inputnode.sphere_reg_dhcpAsym')
                ]),
                (ds_fsLR_reg_wf, fsLR_buffer, [
                    ('outputnode.sphere_reg_dhcpAsym', 'sphere_reg_fsLR'),
                ]),
            ])  # fmt:skip
        else:
            LOGGER.info('ANAT Stage 9: Found pre-computed dhcp-fsLR registration sphere')
            fsLR_buffer.inputs.sphere_reg_fsLR = sorted(precomputed['sphere_reg_dhcpAsym'])

    else:
        if len(precomputed.get('sphere_reg_fsLR', [])) < 2:
            LOGGER.info('ANAT Stage 9: Creating fsLR registration sphere')
            fsLR_reg_wf = init_fsLR_reg_wf()

            ds_fsLR_reg_wf = init_ds_surfaces_wf(
                output_dir=output_dir,
                surfaces=['sphere_reg_fsLR'],
                name='ds_fsLR_reg_wf',
            )

            workflow.connect([
                (surfaces_buffer, fsLR_reg_wf, [('sphere_reg', 'inputnode.sphere_reg')]),
                (sourcefile_buffer, ds_fsLR_reg_wf, [
                    ('anat_source_files', 'inputnode.source_files'),
                ]),
                (fsLR_reg_wf, ds_fsLR_reg_wf, [
                    ('outputnode.sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR')
                ]),
                (ds_fsLR_reg_wf, fsLR_buffer, [('outputnode.sphere_reg_fsLR', 'sphere_reg_fsLR')]),
            ])  # fmt:skip
        else:
            LOGGER.info('ANAT Stage 9: Found pre-computed fsLR registration sphere')
            fsLR_buffer.inputs.sphere_reg_fsLR = sorted(precomputed['sphere_reg_fsLR'])

    # Stage 10: Cortical surface mask
    if len(precomputed.get('cortex_mask', [])) < 2:
        LOGGER.info('ANAT Stage 11: Creating cortical surface mask')

        cortex_masks_wf = init_cortex_masks_wf()
        ds_cortex_masks_wf = init_ds_surface_masks_wf(
            output_dir=output_dir,
            mask_type='cortex',
            name='ds_cortex_masks_wf',
        )

        workflow.connect([
            (surfaces_buffer, cortex_masks_wf, [
                ('midthickness', 'inputnode.midthickness'),
                ('thickness', 'inputnode.thickness'),
            ]),
            (cortex_masks_wf, ds_cortex_masks_wf, [
                ('outputnode.cortex_masks', 'inputnode.mask_files'),
                ('outputnode.source_files', 'inputnode.source_files'),
            ]),
            (ds_cortex_masks_wf, outputnode, [('outputnode.mask_files', 'cortex_mask')]),
        ])  # fmt:skip
    else:
        LOGGER.info('ANAT Stage 11: Found pre-computed cortical surface mask')
        outputnode.inputs.cortex_mask = sorted(precomputed['cortex_mask'])
    return workflow
