# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright The NiPreps Developers <nipreps@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# We support and encourage derived works from this project, please read
# about our expectations at
#
#     https://www.nipreps.org/community/licensing/
#
"""
Orchestrating the BOLD-preprocessing workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_bold_wf
.. autofunction:: init_bold_fit_wf
.. autofunction:: init_bold_native_wf

"""

import typing as ty

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.utils.connections import listify

from nibabies import config
from nibabies._types import Anatomical
from nibabies.interfaces import DerivativesDataSink
from nibabies.utils.misc import estimate_bold_mem_usage

# BOLD workflows
from nibabies.workflows.bold.apply import init_bold_volumetric_resample_wf
from nibabies.workflows.bold.confounds import init_bold_confs_wf, init_carpetplot_wf
from nibabies.workflows.bold.fit import init_bold_fit_wf, init_bold_native_wf
from nibabies.workflows.bold.outputs import (
    init_ds_bold_native_wf,
    init_ds_volumes_wf,
    prepare_timing_parameters,
)
from nibabies.workflows.bold.resampling import init_bold_surf_wf
from nibabies.workflows.bold.t2s import init_t2s_reporting_wf

if ty.TYPE_CHECKING:
    from niworkflows.utils.spaces import SpatialReferences

DEFAULT_DISMISS_ENTITIES = config.DEFAULT_DISMISS_ENTITIES


def init_bold_wf(
    *,
    reference_anat: Anatomical,
    bold_series: list[str],
    precomputed: dict | None = None,
    fieldmap_id: str | None = None,
    spaces: 'SpatialReferences',
) -> pe.Workflow:
    """
    This workflow controls the functional preprocessing stages of *fMRIPrep*.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.tests import mock_config
            from fmriprep import config
            from fmriprep.workflows.bold.base import init_bold_wf
            with mock_config():
                bold_file = config.execution.bids_dir / "sub-01" / "func" \
                    / "sub-01_task-mixedgamblestask_run-01_bold.nii.gz"
                wf = init_bold_wf(
                    bold_series=[str(bold_file)],
                )

    Parameters
    ----------
    bold_series
        List of paths to NIfTI files.
    precomputed
        Dictionary containing precomputed derivatives to reuse, if possible.
    fieldmap_id
        ID of the fieldmap to use to correct this BOLD series. If :obj:`None`,
        no correction will be applied.

    Inputs
    ------
    anat_preproc
        Bias-corrected structural template image
    anat_mask
        Mask of the skull-stripped template image
    anat_dseg
        Segmentation of preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF)
    anat_tpms
        List of tissue probability maps in anatomical space
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    fsnative2anat_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed subject space to anatomical
    white
        FreeSurfer white matter surfaces, in anat space, collated left, then right
    midthickness
        FreeSurfer mid-thickness surfaces, in anat space, collated left, then right
    pial
        FreeSurfer pial surfaces, in anat space, collated left, then right
    sphere_reg_fsLR
        Registration spheres from fsnative to fsLR space, collated left, then right
    anat_ribbon
        Binary cortical ribbon mask in anat space
    fmap_id
        Unique identifiers to select fieldmap files
    fmap
        List of estimated fieldmaps (collated with fmap_id)
    fmap_ref
        List of fieldmap reference files (collated with fmap_id)
    fmap_coeff
        List of lists of spline coefficient files (collated with fmap_id)
    fmap_mask
        List of fieldmap masks (collated with fmap_id)
    sdc_method
        List of fieldmap correction method names (collated with fmap_id)
    anat2std_xfm
        Transform from anatomical space to standard space
    std_t1w
        T1w reference image in standard space
    std_mask
        Brain (binary) mask of the standard reference image
    std_space
        Value of space entity to be used in standard space output filenames
    std_resolution
        Value of resolution entity to be used in standard space output filenames
    std_cohort
        Value of cohort entity to be used in standard space output filenames
    anat2mni6_xfm
        Transform from anatomical space to MNI152NLin6Asym space
    mni6_mask
        Brain (binary) mask of the MNI152NLin6Asym reference image
    mni2009c2anat_xfm
        Transform from MNI152NLin2009cAsym to anatomical space

    Note that ``anat2std_xfm``, ``std_space``, ``std_resolution``,
    ``std_cohort``, ``std_t1w`` and ``std_mask`` are treated as single
    inputs. In order to resample to multiple target spaces, connect
    these fields to an iterable.

    See Also
    --------

    * :func:`~fmriprep.workflows.bold.fit.init_bold_fit_wf`
    * :func:`~fmriprep.workflows.bold.fit.init_bold_native_wf`
    * :func:`~fmriprep.workflows.bold.apply.init_bold_volumetric_resample_wf`
    * :func:`~fmriprep.workflows.bold.outputs.init_ds_bold_native_wf`
    * :func:`~fmriprep.workflows.bold.outputs.init_ds_volumes_wf`
    * :func:`~fmriprep.workflows.bold.t2s.init_t2s_reporting_wf`
    * :func:`~fmriprep.workflows.bold.resampling.init_bold_surf_wf`
    * :func:`~fmriprep.workflows.bold.resampling.init_bold_fsLR_resampling_wf`
    * :func:`~fmriprep.workflows.bold.resampling.init_bold_grayords_wf`
    * :func:`~fmriprep.workflows.bold.confounds.init_bold_confs_wf`
    * :func:`~fmriprep.workflows.bold.confounds.init_carpetplot_wf`

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    if precomputed is None:
        precomputed = {}
    bold_file = bold_series[0]

    output_dir = str(config.execution.nibabies_dir)
    omp_nthreads = config.nipype.omp_nthreads
    all_metadata = [config.execution.layout.get_metadata(file) for file in bold_series]

    nvols, mem_gb = estimate_bold_mem_usage(bold_file)
    if nvols <= 5 - config.execution.sloppy:
        config.loggers.workflow.warning(
            f'Too short BOLD series (<= 5 timepoints). Skipping processing of <{bold_file}>.'
        )
        return

    config.loggers.workflow.debug(
        'Creating bold processing workflow for <%s> (%.2f GB / %d TRs). '
        'Memory resampled/largemem=%.2f/%.2f GB.',
        bold_file,
        mem_gb['filesize'],
        nvols,
        mem_gb['resampled'],
        mem_gb['largemem'],
    )

    workflow = Workflow(name=_get_wf_name(bold_file, 'bold'))
    workflow.__postdesc__ = """\
All resamplings can be performed with *a single interpolation
step* by composing all the pertinent transformations (i.e. head-motion
transform matrices, susceptibility distortion correction when available,
and co-registrations to anatomical and output spaces).
Gridded (volumetric) resamplings were performed using `nitransforms`,
configured with cubic B-spline interpolation.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                # Anatomical coregistration
                'anat_preproc',
                'anat_mask',
                'anat_dseg',
                'anat_tpms',
                'anat_aseg',
                # FreeSurfer outputs
                'subjects_dir',
                'subject_id',
                'fsnative2anat_xfm',
                'white',
                'midthickness',
                'pial',
                'sphere_reg_fsLR',
                'midthickness_fsLR',
                'cortex_mask',
                'anat_ribbon',
                # Fieldmap registration
                'fmap',
                'fmap_ref',
                'fmap_coeff',
                'fmap_mask',
                'fmap_id',
                'sdc_method',
                # Volumetric templates
                'anat2std_xfm',
                'std_t1w',
                'std_mask',
                'std_space',
                'std_resolution',
                'std_cohort',
                # MNI152NLin6Asym warp, for CIFTI use
                # 'anat2mni6_xfm',
                # 'mni6_mask',
                # MNI152NLin2009cAsym inverse warp, for carpetplotting (NOT USED)
                # 'mni2009c2anat_xfm',
                # MNIInfant <cohort> warp, for CIFTI use
                'anat2mniinfant_xfm',
                'mniinfant2anat_xfm',
                'mniinfant_mask',
            ],
        ),
        name='inputnode',
    )

    #
    # Minimal workflow
    #

    bold_fit_wf = init_bold_fit_wf(
        bold_series=bold_series,
        reference_anat=reference_anat,
        precomputed=precomputed,
        fieldmap_id=fieldmap_id,
        omp_nthreads=omp_nthreads,
    )

    workflow.connect([
        (inputnode, bold_fit_wf, [
            ('anat_preproc', 'inputnode.anat_preproc'),
            ('anat_mask', 'inputnode.anat_mask'),
            ('anat_dseg', 'inputnode.anat_dseg'),
            ('subjects_dir', 'inputnode.subjects_dir'),
            ('subject_id', 'inputnode.subject_id'),
            ('fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ('fmap', 'inputnode.fmap'),
            ('fmap_ref', 'inputnode.fmap_ref'),
            ('fmap_coeff', 'inputnode.fmap_coeff'),
            ('fmap_mask', 'inputnode.fmap_mask'),
            ('fmap_id', 'inputnode.fmap_id'),
            ('sdc_method', 'inputnode.sdc_method'),
        ]),
    ])  # fmt:skip

    if config.workflow.level == 'minimal':
        return workflow

    # Now that we're resampling and combining, multiecho matters
    multiecho = len(bold_series) > 2

    nonstd_spaces = set(spaces.get_nonstandard())
    freesurfer_spaces = spaces.get_fs_spaces()

    #
    # Resampling outputs workflow:
    #   - Resample to native
    #   - Save native outputs/echos only if requested
    #

    bold_native_wf = init_bold_native_wf(
        bold_series=bold_series,
        fieldmap_id=fieldmap_id,
        omp_nthreads=omp_nthreads,
    )

    workflow.connect([
        (inputnode, bold_native_wf, [
            ('fmap_ref', 'inputnode.fmap_ref'),
            ('fmap_coeff', 'inputnode.fmap_coeff'),
            ('fmap_id', 'inputnode.fmap_id'),
        ]),
        (bold_fit_wf, bold_native_wf, [
            ('outputnode.coreg_boldref', 'inputnode.boldref'),
            ('outputnode.bold_mask', 'inputnode.bold_mask'),
            ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
            ('outputnode.boldref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
            ('outputnode.dummy_scans', 'inputnode.dummy_scans'),
        ]),
    ])  # fmt:skip

    boldref_out = bool(nonstd_spaces.intersection(('func', 'run', 'bold', 'boldref', 'sbref')))
    boldref_out &= config.workflow.level == 'full'
    echos_out = multiecho and config.execution.me_output_echos

    if boldref_out or echos_out:
        ds_bold_native_wf = init_ds_bold_native_wf(
            bids_root=str(config.execution.bids_dir),
            output_dir=output_dir,
            bold_output=boldref_out,
            echo_output=echos_out,
            multiecho=multiecho,
            all_metadata=all_metadata,
        )
        ds_bold_native_wf.inputs.inputnode.source_files = bold_series

        workflow.connect([
            (bold_fit_wf, ds_bold_native_wf, [
                ('outputnode.bold_mask', 'inputnode.bold_mask'),
                ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
                ('outputnode.boldref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
            ]),
            (bold_native_wf, ds_bold_native_wf, [
                ('outputnode.bold_native', 'inputnode.bold'),
                ('outputnode.bold_echos', 'inputnode.bold_echos'),
                ('outputnode.t2star_map', 'inputnode.t2star'),
            ]),
        ])  # fmt:skip

    if multiecho:
        t2s_reporting_wf = init_t2s_reporting_wf()

        ds_report_t2scomp = pe.Node(
            DerivativesDataSink(
                desc='t2scomp',
                datatype='figures',
                dismiss_entities=DEFAULT_DISMISS_ENTITIES,
            ),
            name='ds_report_t2scomp',
            run_without_submitting=True,
        )

        ds_report_t2star_hist = pe.Node(
            DerivativesDataSink(
                desc='t2starhist',
                datatype='figures',
                dismiss_entities=DEFAULT_DISMISS_ENTITIES,
            ),
            name='ds_report_t2star_hist',
            run_without_submitting=True,
        )

        workflow.connect([
            (inputnode, t2s_reporting_wf, [('anat_dseg', 'inputnode.label_file')]),
            (bold_fit_wf, t2s_reporting_wf, [
                ('outputnode.boldref2anat_xfm', 'inputnode.boldref2anat_xfm'),
                ('outputnode.coreg_boldref', 'inputnode.boldref'),
            ]),
            (bold_native_wf, t2s_reporting_wf, [
                ('outputnode.t2star_map', 'inputnode.t2star_file'),
            ]),
            (t2s_reporting_wf, ds_report_t2scomp, [('outputnode.t2s_comp_report', 'in_file')]),
            (t2s_reporting_wf, ds_report_t2star_hist, [('outputnode.t2star_hist', 'in_file')]),
        ])  # fmt:skip

    if config.workflow.level == 'resampling':
        # Fill-in datasinks of reportlets seen so far
        for node in workflow.list_node_names():
            if node.split('.')[-1].startswith('ds_report'):
                workflow.get_node(node).inputs.base_directory = output_dir
                workflow.get_node(node).inputs.source_file = bold_file
        return workflow

    # Resample to anatomical space
    bold_anat_wf = init_bold_volumetric_resample_wf(
        metadata=all_metadata[0],
        fieldmap_id=fieldmap_id if not multiecho else None,
        omp_nthreads=omp_nthreads,
        mem_gb=mem_gb,
        jacobian='fmap-jacobian' not in config.workflow.ignore,
        name='bold_anat_wf',
    )
    bold_anat_wf.inputs.inputnode.resolution = 'native'

    workflow.connect([
        (inputnode, bold_anat_wf, [
            ('anat_preproc', 'inputnode.target_ref_file'),
            ('anat_mask', 'inputnode.target_mask'),
            ('fmap_ref', 'inputnode.fmap_ref'),
            ('fmap_coeff', 'inputnode.fmap_coeff'),
            ('fmap_id', 'inputnode.fmap_id'),
        ]),
        (bold_fit_wf, bold_anat_wf, [
            ('outputnode.coreg_boldref', 'inputnode.bold_ref_file'),
            ('outputnode.boldref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
            ('outputnode.boldref2anat_xfm', 'inputnode.boldref2anat_xfm'),
        ]),
        (bold_native_wf, bold_anat_wf, [
            ('outputnode.bold_minimal', 'inputnode.bold_file'),
            ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
        ]),
    ])  # fmt:skip

    # Full derivatives, including resampled BOLD series
    if nonstd_spaces.intersection(('anat', 'T1w', 'T2w')):
        ds_bold_anat_wf = init_ds_volumes_wf(
            bids_root=str(config.execution.bids_dir),
            output_dir=output_dir,
            multiecho=multiecho,
            metadata=all_metadata[0],
            name='ds_bold_anat_wf',
        )
        ds_bold_anat_wf.inputs.inputnode.source_files = bold_series
        ds_bold_anat_wf.inputs.inputnode.space = reference_anat

        workflow.connect([
            (bold_fit_wf, ds_bold_anat_wf, [
                ('outputnode.bold_mask', 'inputnode.bold_mask'),
                ('outputnode.coreg_boldref', 'inputnode.bold_ref'),
                ('outputnode.boldref2anat_xfm', 'inputnode.boldref2anat_xfm'),
                ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
                ('outputnode.boldref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
            ]),
            (bold_native_wf, ds_bold_anat_wf, [('outputnode.t2star_map', 'inputnode.t2star')]),
            (bold_anat_wf, ds_bold_anat_wf, [
                ('outputnode.bold_file', 'inputnode.bold'),
                ('outputnode.resampling_reference', 'inputnode.ref_file'),
            ]),
        ])  # fmt:skip

    if spaces.cached.get_spaces(nonstandard=False, dim=(3,)):
        # Missing:
        #  * Clipping BOLD after resampling
        #  * Resampling parcellations
        bold_std_wf = init_bold_volumetric_resample_wf(
            metadata=all_metadata[0],
            fieldmap_id=fieldmap_id if not multiecho else None,
            omp_nthreads=omp_nthreads,
            mem_gb=mem_gb,
            jacobian='fmap-jacobian' not in config.workflow.ignore,
            name='bold_std_wf',
        )
        ds_bold_std_wf = init_ds_volumes_wf(
            bids_root=str(config.execution.bids_dir),
            output_dir=output_dir,
            multiecho=multiecho,
            metadata=all_metadata[0],
            name='ds_bold_std_wf',
        )
        ds_bold_std_wf.inputs.inputnode.source_files = bold_series

        workflow.connect([
            (inputnode, bold_std_wf, [
                ('std_t1w', 'inputnode.target_ref_file'),
                ('std_mask', 'inputnode.target_mask'),
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('std_resolution', 'inputnode.resolution'),
                ('fmap_ref', 'inputnode.fmap_ref'),
                ('fmap_coeff', 'inputnode.fmap_coeff'),
                ('fmap_id', 'inputnode.fmap_id'),
            ]),
            (bold_fit_wf, bold_std_wf, [
                ('outputnode.coreg_boldref', 'inputnode.bold_ref_file'),
                ('outputnode.boldref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
                ('outputnode.boldref2anat_xfm', 'inputnode.boldref2anat_xfm'),
            ]),
            (bold_native_wf, bold_std_wf, [
                ('outputnode.bold_minimal', 'inputnode.bold_file'),
                ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
            ]),
            (inputnode, ds_bold_std_wf, [
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('std_t1w', 'inputnode.template'),
                ('std_space', 'inputnode.space'),
                ('std_resolution', 'inputnode.resolution'),
                ('std_cohort', 'inputnode.cohort'),
            ]),
            (bold_fit_wf, ds_bold_std_wf, [
                ('outputnode.bold_mask', 'inputnode.bold_mask'),
                ('outputnode.coreg_boldref', 'inputnode.bold_ref'),
                ('outputnode.boldref2anat_xfm', 'inputnode.boldref2anat_xfm'),
                ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
                ('outputnode.boldref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
            ]),
            (bold_native_wf, ds_bold_std_wf, [('outputnode.t2star_map', 'inputnode.t2star')]),
            (bold_std_wf, ds_bold_std_wf, [
                ('outputnode.bold_file', 'inputnode.bold'),
                ('outputnode.resampling_reference', 'inputnode.ref_file'),
            ]),
        ])  # fmt:skip

    if config.workflow.surface_recon_method and freesurfer_spaces:
        workflow.__postdesc__ += """\
Non-gridded (surface) resamplings were performed using `mri_vol2surf`
(FreeSurfer).
"""
        config.loggers.workflow.debug('Creating BOLD surface-sampling workflow.')
        bold_surf_wf = init_bold_surf_wf(
            mem_gb=mem_gb['resampled'],
            surface_spaces=freesurfer_spaces,
            medial_surface_nan=config.workflow.medial_surface_nan,
            metadata=all_metadata[0],
            output_dir=output_dir,
            name='bold_surf_wf',
        )
        bold_surf_wf.inputs.inputnode.source_file = bold_file
        workflow.connect([
            (inputnode, bold_surf_wf, [
                ('subjects_dir', 'inputnode.subjects_dir'),
                ('subject_id', 'inputnode.subject_id'),
                ('fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
            ]),
            (bold_anat_wf, bold_surf_wf, [('outputnode.bold_file', 'inputnode.bold_anat')]),
        ])  # fmt:skip

        # sources are bold_file, motion_xfm, boldref2anat_xfm, fsnative2anat_xfm
        merge_surface_sources = pe.Node(
            niu.Merge(4),
            name='merge_surface_sources',
            run_without_submitting=True,
        )
        merge_surface_sources.inputs.in1 = bold_file
        workflow.connect([
            (bold_fit_wf, merge_surface_sources, [
                ('outputnode.motion_xfm', 'in2'),
                ('outputnode.boldref2anat_xfm', 'in3'),
            ]),
            (inputnode, merge_surface_sources, [
                ('fsnative2anat_xfm', 'in4'),
            ]),
            (merge_surface_sources, bold_surf_wf, [
                ('out', 'inputnode.sources'),
            ]),
        ])  # fmt:skip

    cifti_output = config.workflow.cifti_output
    if cifti_output:
        from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms

        from nibabies.workflows.bold.alignment import (
            init_subcortical_mni_alignment_wf,
            init_subcortical_rois_wf,
        )
        from nibabies.workflows.bold.resampling import (
            init_bold_fsLR_resampling_wf,
            init_bold_grayords_wf,
            init_goodvoxels_bold_mask_wf,
        )

        bold_MNIInfant_wf = init_bold_volumetric_resample_wf(
            metadata=all_metadata[0],
            fieldmap_id=fieldmap_id if not multiecho else None,
            omp_nthreads=omp_nthreads,
            mem_gb=mem_gb,
            jacobian='fmap-jacobian' not in config.workflow.ignore,
            name='bold_MNIInfant_wf',
        )

        bold_fsLR_resampling_wf = init_bold_fsLR_resampling_wf(
            grayord_density=cifti_output,
            omp_nthreads=omp_nthreads,
            mem_gb=mem_gb['resampled'],
        )

        if config.workflow.project_goodvoxels:
            goodvoxels_bold_mask_wf = init_goodvoxels_bold_mask_wf(mem_gb['resampled'])

            workflow.connect([
                (inputnode, goodvoxels_bold_mask_wf, [('anat_ribbon', 'inputnode.anat_ribbon')]),
                (bold_anat_wf, goodvoxels_bold_mask_wf, [
                    ('outputnode.bold_file', 'inputnode.bold_file'),
                ]),
                (goodvoxels_bold_mask_wf, bold_fsLR_resampling_wf, [
                    ('outputnode.goodvoxels_mask', 'inputnode.volume_roi'),
                ]),
            ])  # fmt:skip

            bold_fsLR_resampling_wf.__desc__ += """\
A "goodvoxels" mask was applied during volume-to-surface sampling in fsLR space,
excluding voxels whose time-series have a locally high coefficient of variation.
"""

        # MNIInfant -> MNI6 registrations (per ROI)
        MNIInfant_aseg = pe.Node(
            ApplyTransforms(interpolation='MultiLabel'),
            name='MNIInfant_aseg',
            mem_gb=1,
        )

        subcortical_rois_wf = init_subcortical_rois_wf()
        subcortical_mni_alignment_wf = init_subcortical_mni_alignment_wf()

        bold_grayords_wf = init_bold_grayords_wf(
            grayord_density=cifti_output,
            repetition_time=all_metadata[0]['RepetitionTime'],
        )

        ds_bold_cifti = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                dismiss_entities=DEFAULT_DISMISS_ENTITIES,
                space='fsLR',
                density=cifti_output,
                suffix='bold',
                compress=False,
                TaskName=all_metadata[0].get('TaskName'),
                **prepare_timing_parameters(all_metadata[0]),
            ),
            name='ds_bold_cifti',
            run_without_submitting=True,
        )
        ds_bold_cifti.inputs.source_file = bold_file

        mniinfant_res = 2 if cifti_output == '91k' else 1
        inputnode.inputs.mniinfant_mask = get_MNIInfant_mask(spaces, mniinfant_res)

        workflow.connect([
            # Resample BOLD to MNI152NLin6Asym, may duplicate bold_std_wf above
            (inputnode, bold_MNIInfant_wf, [
                ('mniinfant_mask', 'inputnode.target_ref_file'),
                ('mniinfant_mask', 'inputnode.target_mask'),
                ('anat2mniinfant_xfm', 'inputnode.anat2std_xfm'),
                ('fmap_ref', 'inputnode.fmap_ref'),
                ('fmap_coeff', 'inputnode.fmap_coeff'),
                ('fmap_id', 'inputnode.fmap_id'),
            ]),
            (bold_fit_wf, bold_MNIInfant_wf, [
                ('outputnode.coreg_boldref', 'inputnode.bold_ref_file'),
                ('outputnode.boldref2fmap_xfm', 'inputnode.boldref2fmap_xfm'),
                ('outputnode.boldref2anat_xfm', 'inputnode.boldref2anat_xfm'),
            ]),
            (bold_native_wf, bold_MNIInfant_wf, [
                ('outputnode.bold_minimal', 'inputnode.bold_file'),
                ('outputnode.motion_xfm', 'inputnode.motion_xfm'),
            ]),
            (inputnode, MNIInfant_aseg, [
                ('anat2mniinfant_xfm', 'transforms'),
                ('anat_aseg', 'input_image'),
            ]),
            (bold_MNIInfant_wf, MNIInfant_aseg, [
                ('outputnode.resampling_reference', 'reference_image'),
            ]),
            (bold_MNIInfant_wf, subcortical_mni_alignment_wf, [
                ('outputnode.bold_file', 'inputnode.MNIInfant_bold'),
            ]),
            (MNIInfant_aseg, subcortical_rois_wf, [
                ('output_image', 'inputnode.MNIInfant_aseg'),
            ]),
            (subcortical_rois_wf, subcortical_mni_alignment_wf, [
                ('outputnode.MNIInfant_rois', 'inputnode.MNIInfant_rois'),
                ('outputnode.MNI152_rois', 'inputnode.MNI152_rois'),
            ]),
            (subcortical_mni_alignment_wf, bold_grayords_wf, [
                ('outputnode.subcortical_volume', 'inputnode.bold_std'),
                ('outputnode.subcortical_labels', 'inputnode.bold_labels'),
            ]),

            # Resample anat-space BOLD to fsLR surfaces
            (inputnode, bold_fsLR_resampling_wf, [
                ('white', 'inputnode.white'),
                ('pial', 'inputnode.pial'),
                ('midthickness', 'inputnode.midthickness'),
                ('midthickness_fsLR', 'inputnode.midthickness_fsLR'),
                ('sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR'),
                ('cortex_mask', 'inputnode.cortex_mask'),
            ]),
            (bold_anat_wf, bold_fsLR_resampling_wf, [
                ('outputnode.bold_file', 'inputnode.bold_file'),
            ]),
            # (bold_MNI6_wf, bold_grayords_wf, [
            #     ('outputnode.bold_file', 'inputnode.bold_std'),
            # ]),
            (bold_fsLR_resampling_wf, bold_grayords_wf, [
                ('outputnode.bold_fsLR', 'inputnode.bold_fsLR'),
            ]),
            (bold_grayords_wf, ds_bold_cifti, [
                ('outputnode.dtseries', 'in_file'),
                (('outputnode.dtseries_metadata', _read_json), 'meta_dict'),
            ]),
        ])  # fmt:skip

    bold_confounds_wf = init_bold_confs_wf(
        mem_gb=mem_gb['largemem'],
        metadata=all_metadata[0],
        freesurfer=config.workflow.run_reconall,
        regressors_all_comps=config.workflow.regressors_all_comps,
        regressors_fd_th=config.workflow.regressors_fd_th,
        regressors_dvars_th=config.workflow.regressors_dvars_th,
        name='bold_confounds_wf',
    )

    ds_confounds = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='confounds',
            suffix='timeseries',
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
        ),
        name='ds_confounds',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    ds_confounds.inputs.source_file = bold_file

    workflow.connect([
        (inputnode, bold_confounds_wf, [
            ('anat_tpms', 'inputnode.anat_tpms'),
            ('anat_mask', 'inputnode.anat_mask'),
        ]),
        (bold_fit_wf, bold_confounds_wf, [
            ('outputnode.bold_mask', 'inputnode.bold_mask'),
            ('outputnode.movpar_file', 'inputnode.movpar_file'),
            ('outputnode.rmsd_file', 'inputnode.rmsd_file'),
            ('outputnode.boldref2anat_xfm', 'inputnode.boldref2anat_xfm'),
            ('outputnode.dummy_scans', 'inputnode.skip_vols'),
        ]),
        (bold_native_wf, bold_confounds_wf, [
            ('outputnode.bold_native', 'inputnode.bold'),
        ]),
        (bold_confounds_wf, ds_confounds, [
            ('outputnode.confounds_file', 'in_file'),
            ('outputnode.confounds_metadata', 'meta_dict'),
        ]),
    ])  # fmt:skip

    # MG: Carpetplot workflow only work with CIFTI
    if cifti_output:
        carpetplot_wf = init_carpetplot_wf(
            mem_gb=mem_gb['resampled'],
            metadata=all_metadata[0],
            cifti_output=cifti_output,
            name='carpetplot_wf',
        )

        def _last(inlist):
            return inlist[-1]

        workflow.connect([
            (bold_grayords_wf, carpetplot_wf, [
                ('outputnode.dtseries', 'inputnode.cifti_bold'),
            ]),
            (bold_fit_wf, carpetplot_wf, [
                ('outputnode.dummy_scans', 'inputnode.dummy_scans'),
                ('outputnode.bold_mask', 'inputnode.bold_mask'),
                ('outputnode.boldref2anat_xfm', 'inputnode.boldref2anat_xfm'),
            ]),
            (bold_native_wf, carpetplot_wf, [
                ('outputnode.bold_native', 'inputnode.bold'),
            ]),
            (bold_confounds_wf, carpetplot_wf, [
                ('outputnode.confounds_file', 'inputnode.confounds_file'),
                ('outputnode.crown_mask', 'inputnode.crown_mask'),
                (('outputnode.acompcor_masks', _last), 'inputnode.acompcor_mask'),
            ]),
        ])  # fmt:skip

    # Fill-in datasinks of reportlets seen so far
    for node in workflow.list_node_names():
        if node.split('.')[-1].startswith('ds_report'):
            workflow.get_node(node).inputs.base_directory = output_dir
            workflow.get_node(node).inputs.source_file = bold_file

    return workflow


def _get_wf_name(bold_fname, prefix):
    """
    Derive the workflow name for supplied BOLD file.

    >>> _get_wf_name("/completely/made/up/path/sub-01_task-nback_bold.nii.gz", "bold")
    'bold_task_nback_wf'
    >>> _get_wf_name(
    ...     "/completely/made/up/path/sub-01_task-nback_run-01_echo-1_bold.nii.gz",
    ...     "preproc",
    ... )
    'preproc_task_nback_run_01_echo_1_wf'

    """
    from nipype.utils.filemanip import split_filename

    fname = split_filename(bold_fname)[1]
    fname_nosub = '_'.join(fname.split('_')[1:-1])
    return f'{prefix}_{fname_nosub.replace("-", "_")}_wf'


def extract_entities(file_list):
    """
    Return a dictionary of common entities given a list of files.

    Examples
    --------
    >>> extract_entities("sub-01/anat/sub-01_T1w.nii.gz")
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(["sub-01/anat/sub-01_T1w.nii.gz"] * 2)
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(["sub-01/anat/sub-01_run-1_T1w.nii.gz",
    ...                   "sub-01/anat/sub-01_run-2_T1w.nii.gz"])
    {'subject': '01', 'run': [1, 2], 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}

    """
    from collections import defaultdict

    from bids.layout import parse_file_entities

    entities = defaultdict(list)
    for e, v in [
        ev_pair for f in listify(file_list) for ev_pair in parse_file_entities(f).items()
    ]:
        entities[e].append(v)

    def _unique(inlist):
        inlist = sorted(set(inlist))
        if len(inlist) == 1:
            return inlist[0]
        return inlist

    return {k: _unique(v) for k, v in entities.items()}


def _read_json(in_file):
    from json import loads
    from pathlib import Path

    return loads(Path(in_file).read_text())


def get_MNIInfant_mask(spaces: 'SpatialReferences', res: str | int) -> str:
    """Parse spaces and return matching MNIInfant space, including cohort."""
    import templateflow.api as tf

    for ref in spaces.references:
        if ref.space == 'MNIInfant' and f'res-{res}' in str(ref):
            return str(
                tf.get(
                    'MNIInfant',
                    cohort=ref.spec['cohort'],
                    resolution=res,
                    desc='brain',
                    suffix='mask',
                )
            )

    raise FileNotFoundError(f'MNIInfant mask (resolution {res}) not found.')
