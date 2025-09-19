import typing as ty

import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
from niworkflows.engine import tag
from smriprep.workflows.outputs import (
    init_ds_anat_volumes_wf,
    init_ds_grayord_metrics_wf,
    init_ds_surface_metrics_wf,
    init_ds_surfaces_wf,
)
from smriprep.workflows.surfaces import (
    init_hcp_morphometrics_wf,
    init_morph_grayords_wf,
    init_resample_surfaces_wf,
    init_surface_derivatives_wf,
)

from nibabies import config
from nibabies.workflows.anatomical.outputs import init_ds_seg_wf
from nibabies.workflows.anatomical.surfaces import init_resample_surfaces_dhcp_wf

if ty.TYPE_CHECKING:
    from niworkflows.utils.spaces import SpatialReferences

LOGGER = config.loggers.workflow


@tag('anat.apply')
def init_infant_anat_apply_wf(
    *,
    bids_root: str,
    cifti_output: ty.Literal['91k', '170k', False],
    msm_sulc: bool,
    omp_nthreads: int,
    output_dir: str,
    recon_method: ty.Literal['freesurfer', 'infantfs', 'mcribs', None],
    reference_anat: ty.Literal['T1w', 'T2w'],
    precomputed: dict,
    spaces: 'SpatialReferences',
    name: str = 'infant_anat_apply_wf',
) -> pe.Workflow:
    """Save superfluous outputs."""
    workflow = pe.Workflow(name=name)

    reg_sphere = f'sphere_reg_{"msm" if msm_sulc else "fsLR"}'
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'anat2std_xfm',
                'anat_valid_list',
                'anat_preproc',
                'anat_mask',
                'anat_dseg',
                'anat_tpms',
                'subjects_dir',
                'subject_id',
                'fsnative2anat_xfm',
                'sulc',
                'template',
                'thickness',
                'white',
                'pial',
                'midthickness',
                'cortex_mask',
                reg_sphere,
                # template workflow inputs
                'std_t1w',
                'anat2std_xfm',
                'std_space',
                'std_cohort',
                'std_resolution',
            ]
        ),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['anat_aseg', 'anat_aparc', 'midthickness_fsLR']),
        name='outputnode',
    )

    if spaces.cached.get_spaces(nonstandard=False, dim=(3,)):
        ds_std_volumes_wf = init_ds_anat_volumes_wf(
            bids_root=bids_root,
            output_dir=output_dir,
            name='ds_std_volumes_wf',
        )

        workflow.connect([
            (inputnode, ds_std_volumes_wf, [
                ('anat_valid_list', 'inputnode.source_files'),
                ('anat_preproc', 'inputnode.anat_preproc'),
                ('anat_mask', 'inputnode.anat_mask'),
                ('anat_dseg', 'inputnode.anat_dseg'),
                ('anat_tpms', 'inputnode.anat_tpms'),
            ]),
            (inputnode, ds_std_volumes_wf, [
                ('std_t1w', 'inputnode.ref_file'),
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('std_space', 'inputnode.space'),
                ('std_cohort', 'inputnode.cohort'),
                ('std_resolution', 'inputnode.resolution'),
            ]),
        ])  # fmt:skip

    if recon_method is not None:
        anat_aseg = precomputed.get(f'{reference_anat.lower()}_aseg')
        anat_aparc = precomputed.get(f'{reference_anat.lower()}_aparc')

        surface_derivatives_wf = init_surface_derivatives_wf()
        # Split out segmentations to better save precomputed inputs
        seg_buffer = pe.Node(
            niu.IdentityInterface(fields=['anat_aseg', 'anat_aparc']),
            name='seg_buffer',
        )

        ds_aseg_wf = init_ds_seg_wf(
            output_dir=output_dir,
            seg_type='aseg',
            extra_entities={'space': reference_anat},
        )
        ds_aparc_wf = init_ds_seg_wf(
            output_dir=output_dir,
            seg_type='aparcaseg',
            extra_entities={'space': reference_anat},
        )
        if anat_aseg:
            LOGGER.info('ANAT - using precomputed aseg')
            seg_buffer.inputs.anat_aseg = anat_aseg
        else:
            workflow.connect([
                (surface_derivatives_wf, seg_buffer, [
                    ('outputnode.out_aseg', 'anat_aseg'),
                ]),
            ])  # fmt:skip

        if anat_aparc:
            LOGGER.info('ANAT - using precomputed aparc')
            seg_buffer.inputs.anat_aparc = anat_aparc
        else:
            workflow.connect([
                (surface_derivatives_wf, seg_buffer, [
                    ('outputnode.out_aparc', 'anat_aparc'),
                ]),
            ])  # fmt:skip

        ds_surfaces_wf = init_ds_surfaces_wf(output_dir=output_dir, surfaces=['inflated'])
        ds_curv_wf = init_ds_surface_metrics_wf(
            bids_root=bids_root, output_dir=output_dir, metrics=['curv'], name='ds_curv_wf'
        )

        workflow.connect([
            (inputnode, surface_derivatives_wf, [
                ('anat_preproc', 'inputnode.reference'),
                ('subjects_dir', 'inputnode.subjects_dir'),
                ('subject_id', 'inputnode.subject_id'),
                ('fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (inputnode, ds_surfaces_wf, [
                ('anat_valid_list', 'inputnode.source_files'),
            ]),
            (surface_derivatives_wf, ds_surfaces_wf, [
                ('outputnode.inflated', 'inputnode.inflated'),
            ]),
            (inputnode, ds_curv_wf, [
                ('anat_valid_list', 'inputnode.source_files'),
            ]),
            (surface_derivatives_wf, ds_curv_wf, [
                ('outputnode.curv', 'inputnode.curv'),
            ]),
            (inputnode, ds_aseg_wf, [
                ('anat_valid_list', 'inputnode.source_files'),
            ]),
            (inputnode, ds_aparc_wf, [
                ('anat_valid_list', 'inputnode.source_files'),
            ]),
            (seg_buffer, outputnode, [
                ('anat_aparc', 'anat_aparc'),
                ('anat_aseg', 'anat_aseg'),
            ]),
            (seg_buffer, ds_aparc_wf, [
                ('anat_aparc', 'inputnode.in_seg'),
            ]),
            (seg_buffer, ds_aseg_wf, [
                ('anat_aseg', 'inputnode.in_seg'),
            ])
        ])  # fmt:skip

        if cifti_output:
            hcp_morphometrics_wf = init_hcp_morphometrics_wf(omp_nthreads=omp_nthreads)
            if recon_method == 'mcribs':
                resample_surfaces_wf = init_resample_surfaces_dhcp_wf(
                    surfaces=['white', 'pial', 'midthickness'],
                    grayord_density=cifti_output,
                )
            else:
                resample_surfaces_wf = init_resample_surfaces_wf(
                    surfaces=['white', 'pial', 'midthickness'], grayord_density=cifti_output
                )
            morph_grayords_wf = init_morph_grayords_wf(
                grayord_density=cifti_output, omp_nthreads=omp_nthreads
            )

            ds_fsLR_surfaces_wf = init_ds_surfaces_wf(
                output_dir=output_dir,
                surfaces=['white', 'pial', 'midthickness'],
                entities={
                    'space': 'dhcpAsym' if recon_method == 'mcribs' else 'fsLR',
                    'density': '32k' if cifti_output == '91k' else '59k',
                },
                name='ds_fsLR_surfaces_wf',
            )

            ds_grayord_metrics_wf = init_ds_grayord_metrics_wf(
                bids_root=bids_root,
                output_dir=output_dir,
                metrics=['curv', 'thickness', 'sulc'],
                cifti_output=cifti_output,
            )

            workflow.connect([
                (inputnode, hcp_morphometrics_wf, [
                    ('subject_id', 'inputnode.subject_id'),
                    ('sulc', 'inputnode.sulc'),
                    ('thickness', 'inputnode.thickness'),
                    ('midthickness', 'inputnode.midthickness'),
                ]),
                (surface_derivatives_wf, hcp_morphometrics_wf, [
                    ('outputnode.curv', 'inputnode.curv'),
                ]),
                (inputnode, resample_surfaces_wf, [
                    ('white', 'inputnode.white'),
                    ('pial', 'inputnode.pial'),
                    ('midthickness', 'inputnode.midthickness'),
                    (reg_sphere, 'inputnode.sphere_reg_fsLR'),
                ]),
                (inputnode, morph_grayords_wf, [
                    ('cortex_mask', 'inputnode.roi'),
                    ('midthickness', 'inputnode.midthickness'),
                    (reg_sphere, 'inputnode.sphere_reg_fsLR'),
                ]),
                (hcp_morphometrics_wf, morph_grayords_wf, [
                    ('outputnode.curv', 'inputnode.curv'),
                    ('outputnode.sulc', 'inputnode.sulc'),
                    ('outputnode.thickness', 'inputnode.thickness'),
                ]),
                (resample_surfaces_wf, morph_grayords_wf, [
                    ('outputnode.midthickness_fsLR', 'inputnode.midthickness_fsLR'),
                ]),
                (inputnode, ds_fsLR_surfaces_wf, [
                    ('anat_valid_list', 'inputnode.source_files'),
                ]),
                (resample_surfaces_wf, outputnode, [
                    ('outputnode.midthickness_fsLR', 'midthickness_fsLR'),
                ]),
                (resample_surfaces_wf, ds_fsLR_surfaces_wf, [
                    ('outputnode.white_fsLR', 'inputnode.white'),
                    ('outputnode.pial_fsLR', 'inputnode.pial'),
                    ('outputnode.midthickness_fsLR', 'inputnode.midthickness'),
                ]),
                (inputnode, ds_grayord_metrics_wf, [
                    ('anat_valid_list', 'inputnode.source_files'),
                ]),
                (morph_grayords_wf, ds_grayord_metrics_wf, [
                    ('outputnode.curv_fsLR', 'inputnode.curv'),
                    ('outputnode.curv_metadata', 'inputnode.curv_metadata'),
                    ('outputnode.thickness_fsLR', 'inputnode.thickness'),
                    ('outputnode.thickness_metadata', 'inputnode.thickness_metadata'),
                    ('outputnode.sulc_fsLR', 'inputnode.sulc'),
                    ('outputnode.sulc_metadata', 'inputnode.sulc_metadata'),
                ]),
            ])  # fmt:skip

    return workflow
