# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# STATEMENT OF CHANGES: This file is derived from sources licensed under the Apache-2.0 terms,
# and this file has been changed.
# The original file this work derives from is found at:
# https://github.com/nipreps/fmriprep/blob/a4fd718/fmriprep/workflows/bold/base.py
#
# [January 2022] CHANGES:
#   * `init_nibabies_wf` now takes in a dictionary composed of participant/session key/values.
#   * `init_single_subject_wf` now differentiates between participant sessions.
#     This change is to treat sessions as a "first-class" identifier, to better handle the
#     potential rapid changing of brain morphometry.
#
# Copyright 2023 The NiPreps Developers <nipreps@gmail.com>
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
NiBabies base processing workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_nibabies_wf
.. autofunction:: init_single_subject_wf

"""
from __future__ import annotations

import os
import sys
import typing as ty
from copy import deepcopy

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from packaging.version import Version

from nibabies import config
from nibabies.interfaces import DerivativesDataSink
from nibabies.interfaces.reports import AboutSummary, SubjectSummary
from nibabies.utils.bids import parse_bids_for_age_months
from nibabies.workflows.anatomical.fit import init_infant_anat_fit_wf, init_infant_anat_full_wf

if ty.TYPE_CHECKING:
    from bids.layout import BIDSLayout
    from niworkflows.utils.spaces import SpatialReferences


def init_nibabies_wf(subworkflows_list):
    """
    Build *NiBabies*'s pipeline.

    This workflow organizes the execution of NiBabies, with a sub-workflow for
    each subject.

    If FreeSurfer's ``infant_recon_all`` is to be run, a corresponding folder is created
    and populated with any needed template subjects under the derivatives folder.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from nibabies.workflows.tests import mock_config
            from nibabies.workflows.base import init_nibabies_wf
            with mock_config():
                wf = init_nibabies_wf()

    Parameters
    ----------
    subworkflows_list: :obj:`list` of :obj:`tuple`
        A list of the subworkflows to create.
        Each subject session is run as an individual workflow.
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.bids import BIDSFreeSurferDir

    ver = Version(config.environment.version)
    nibabies_wf = Workflow(name=f'nibabies_{ver.major}_{ver.minor}_wf')
    nibabies_wf.base_dir = config.execution.work_dir

    execution_spaces = init_execution_spaces()

    freesurfer = config.workflow.run_reconall
    if freesurfer:
        fsdir = pe.Node(
            BIDSFreeSurferDir(
                derivatives=config.execution.output_dir,
                freesurfer_home=os.getenv('FREESURFER_HOME'),
                spaces=execution_spaces.get_fs_spaces(),
            ),
            name=f"fsdir_run_{config.execution.run_uuid.replace('-', '_')}",
            run_without_submitting=True,
        )
        if config.execution.fs_subjects_dir is not None:
            fsdir.inputs.subjects_dir = str(config.execution.fs_subjects_dir.absolute())

    for subject_id, session_id in subworkflows_list:
        # Calculate the age and age-specific spaces
        age = parse_bids_for_age_months(config.execution.bids_dir, subject_id, session_id)
        if config.workflow.age_months:
            config.loggers.cli.warning(
                '`--age-months` is deprecated and will be removed in a future release.'
                'Please use a `sessions.tsv` or `participants.tsv` file to track participants age.'
            )
            age = config.workflow.age_months
        if age is None:
            raise RuntimeError(
                'Could not find age for sub-{subject}{session}'.format(
                    subject=subject_id, session=f'_ses-{session_id}' if session_id else ''
                )
            )
        output_spaces = init_workflow_spaces(execution_spaces, age)

        bids_level = [f'sub-{subject_id}']
        if session_id is not None:
            bids_level.append(f'ses-{session_id}')

        log_dir = (
            config.execution.nibabies_dir.joinpath(*bids_level) / 'log' / config.execution.run_uuid
        )
        # skull strip template cohort
        single_subject_wf = init_single_subject_wf(
            subject_id=subject_id,
            age=age,
            session_id=session_id,
            spaces=output_spaces,
        )

        single_subject_wf.config['execution']['crashdump_dir'] = str(log_dir)
        for node in single_subject_wf._get_all_nodes():
            node.config = deepcopy(single_subject_wf.config)
        if freesurfer:
            nibabies_wf.connect(fsdir, 'subjects_dir', single_subject_wf, 'inputnode.subjects_dir')
        else:
            nibabies_wf.add_nodes([single_subject_wf])

        # Dump a copy of the config file into the log directory
        log_dir.mkdir(exist_ok=True, parents=True)
        config.to_filename(log_dir / 'nibabies.toml')

    return nibabies_wf


def init_single_subject_wf(
    *,
    subject_id: str,
    age: int,
    session_id: str | None = None,
    spaces: SpatialReferences | None = None,
):
    """
    Organize the preprocessing pipeline for a single subject, at a single session.

    It collects and reports information about the subject, and prepares
    sub-workflows to perform anatomical and functional preprocessing.
    Anatomical preprocessing is performed in a single workflow, regardless of
    the number of sessions.
    Functional preprocessing is performed using a separate workflow for each
    individual BOLD series.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from nibabies.workflows.tests import mock_config
            from nibabies.workflows.base import init_single_subject_wf
            with mock_config():
                wf = init_single_subject_wf('01')

    Parameters
    ----------
    subject_id : :obj:`str`
        Subject label for this single-subject workflow.
    session_id : :obj:`str` or None
        Session identifier.
    age: :obj:`int` or None
        Age (in months) of subject.

    Inputs
    ------
    subjects_dir : :obj:`str`
        FreeSurfer's ``$SUBJECTS_DIR``.

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.bids import BIDSDataGrabber, BIDSInfo
    from niworkflows.interfaces.nilearn import NILEARN_VERSION
    from niworkflows.utils.bids import collect_data
    from niworkflows.utils.spaces import Reference

    from ..utils.misc import fix_multi_source_name

    subject_session_id = _subject_session_id(subject_id, session_id)
    print(f'{subject_session_id=}')
    workflow = Workflow(name=f'single_subject_{subject_session_id}_wf')
    workflow.__desc__ = f"""
Results included in this manuscript come from preprocessing
performed using *NiBabies* {config.environment.version},
derived from fMRIPrep (@fmriprep1; @fmriprep2; RRID:SCR_016216).
The underlying workflow engine used is *Nipype* {config.environment.nipype_version}
(@nipype1; @nipype2; RRID:SCR_002502).

"""
    workflow.__postdesc__ = f"""

Many internal operations of *NiBabies* use
*Nilearn* {NILEARN_VERSION} [@nilearn, RRID:SCR_001362],
mostly within the functional processing workflow.
For more details of the pipeline, see [the section corresponding
to workflows in *nibabies*'s documentation]\
(https://nibabies.readthedocs.io/en/latest/workflows.html \
"NiBabies's documentation").


### Copyright Waiver

The above boilerplate text was automatically generated by NiBabies
with the express intention that users should copy and paste this
text into their manuscripts *unchanged*.
It is released under the [CC0]\
(https://creativecommons.org/publicdomain/zero/1.0/) license.

### References

"""

    subject_data = collect_data(
        config.execution.layout,
        subject_id,
        session_id=session_id,
        task=config.execution.task_id,
        echo=config.execution.echo_idx,
        bids_filters=config.execution.bids_filters,
    )[0]

    if 'flair' in config.workflow.ignore:
        subject_data['flair'] = []
    if 't1w' in config.workflow.ignore:
        subject_data['t1w'] = []
    if 't2w' in config.workflow.ignore:
        subject_data['t2w'] = []

    anat_only = config.workflow.anat_only
    # Make sure we always go through these two checks
    if not anat_only and not subject_data['bold']:
        task_id = config.execution.task_id
        raise RuntimeError(
            'No BOLD images found for participant {} and task {}. '
            'All workflows require BOLD images.'.format(
                subject_id, task_id if task_id else '<all>'
            )
        )

    # bold_runs = [
    #     sorted(
    #         listify(run),
    #         key=lambda fl: config.execution.layout.get_metadata(fl).get('EchoTime', 0),
    #     )
    #     for run in subject_data['bold']
    # ]

    # if subject_data['roi']:
    #     warnings.warn(
    #         f"Lesion mask {subject_data['roi']} found. "
    #         "Future versions of fMRIPrep will use alternative conventions. "
    #         "Please refer to the documentation before upgrading.",
    #         FutureWarning,
    #         stacklevel=1,
    #     )

    recon_method = config.workflow.surface_recon_method
    msm_sulc = False

    anatomical_cache = {}
    if config.execution.derivatives:
        from nibabies.utils.derivatives import collect_anatomical_derivatives

        std_spaces = spaces.get_spaces(nonstandard=False, dim=(3,))
        std_spaces.append('fsnative')
        for deriv_dir in config.execution.derivatives.values():
            anatomical_cache.update(
                collect_anatomical_derivatives(
                    derivatives_dir=deriv_dir,
                    subject_id=subject_id,
                    session_id=session_id,
                    std_spaces=std_spaces,
                )
            )

    # Determine some session level options here, as we should have
    # all the required information
    if recon_method == 'auto':
        if age <= 8:
            recon_method = 'mcribs'
        elif age <= 24:
            recon_method = 'infantfs'
        else:
            recon_method = 'freesurfer'

    preferred_anat = config.execution.reference_anat
    t1w = subject_data['t1w']
    t2w = subject_data['t2w']
    if not t1w and t2w:
        reference_anat = 'T1w' if t1w else 'T2w'
        if preferred_anat and reference_anat != preferred_anat:
            raise AttributeError(
                f'Requested to use {preferred_anat} as anatomical reference but none available'
            )
    else:
        if not (reference_anat := preferred_anat):
            reference_anat = 'T2w' if age <= 8 else 'T1w'
    anat = reference_anat.lower()  # To be used for workflow connections

    bids_root = str(config.execution.bids_dir)
    output_dir = str(config.execution.nibabies_dir)
    omp_nthreads = config.nipype.omp_nthreads

    inputnode = pe.Node(niu.IdentityInterface(fields=['subjects_dir']), name='inputnode')

    bidssrc = pe.Node(
        BIDSDataGrabber(
            subject_data=subject_data,
            anat_only=config.workflow.anat_only,
            anat_derivatives=anatomical_cache or None,
            subject_id=subject_id,
        ),
        name='bidssrc',
    )

    bids_info = pe.Node(
        BIDSInfo(bids_dir=config.execution.bids_dir, bids_validate=False),
        name='bids_info',
    )

    summary = pe.Node(
        SubjectSummary(
            anatomical_reference=reference_anat,
            recon_method=recon_method,
            std_spaces=spaces.get_spaces(nonstandard=False),
            nstd_spaces=spaces.get_spaces(standard=False),
        ),
        name='summary',
        run_without_submitting=True,
    )

    about = pe.Node(
        AboutSummary(version=config.environment.version, command=' '.join(sys.argv)),
        name='about',
        run_without_submitting=True,
    )

    ds_report_summary = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='summary',
            datatype='figures',
            dismiss_entities=('echo',),
        ),
        name='ds_report_summary',
        run_without_submitting=True,
    )

    ds_report_about = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='about',
            datatype='figures',
            dismiss_entities=('echo',),
        ),
        name='ds_report_about',
        run_without_submitting=True,
    )

    wf_args = {
        'age_months': age,
        't1w': t1w,
        't2w': t2w,
        'flair': subject_data['flair'],
        'bids_root': bids_root,
        'longitudinal': config.workflow.longitudinal,
        'msm_sulc': msm_sulc,
        'omp_nthreads': omp_nthreads,
        'output_dir': config.execution.nibabies_dir,
        'precomputed': anatomical_cache,
        'segmentation_atlases': config.execution.segmentation_atlases_dir,
        'skull_strip_fixed_seed': config.workflow.skull_strip_fixed_seed,
        'skull_strip_mode': config.workflow.skull_strip_anat,
        'skull_strip_template': Reference.from_string(config.workflow.skull_strip_template)[0],
        'recon_method': recon_method,
        'reference_anat': reference_anat,
        'sloppy': config.execution.sloppy,
        'spaces': spaces,
        'cifti_output': config.workflow.cifti_output,
    }

    anat_wf = (
        init_infant_anat_full_wf(**wf_args)
        if config.workflow.level == 'full'
        else init_infant_anat_fit_wf(**wf_args)
    )

    # allow to run with anat-fast-track on fMRI-only dataset
    if (
        't1w_preproc' in anatomical_cache or 't2w_preproc' in anatomical_cache
    ) and not subject_data['t1w']:
        workflow.connect([
            (bidssrc, bids_info, [(('bold', fix_multi_source_name), 'in_file')]),
            (anat_wf, summary, [('outputnode.anat_preproc', anat)]),
            (anat_wf, ds_report_summary, [('outputnode.anat_preproc', 'source_file')]),
            (anat_wf, ds_report_about, [('outputnode.anat_preproc', 'source_file')]),
        ])  # fmt:skip
    else:
        workflow.connect([
            (bidssrc, bids_info, [(('t1w', fix_multi_source_name), 'in_file')]),
            (bidssrc, summary, [('t1w', 't1w')]),
            (bidssrc, ds_report_summary, [(('t1w', fix_multi_source_name), 'source_file')]),
            (bidssrc, ds_report_about, [(('t1w', fix_multi_source_name), 'source_file')]),
        ])  # fmt:skip

    workflow.connect([
        (inputnode, anat_wf, [('subjects_dir', 'inputnode.subjects_dir')]),
        (bidssrc, anat_wf, [
            ('t1w', 'inputnode.t1w'),
            ('t2w', 'inputnode.t2w'),
            ('roi', 'inputnode.roi'),
            ('flair', 'inputnode.flair'),
        ]),
        # Reporting connections
        (inputnode, summary, [('subjects_dir', 'subjects_dir')]),
        (bidssrc, summary, [('t2w', 't2w'), ('bold', 'bold')]),
        (bids_info, summary, [('subject', 'subject_id')]),
        (summary, ds_report_summary, [('out_report', 'in_file')]),
        (summary, anat_wf, [('subject_id', 'inputnode.subject_id')]),
        (about, ds_report_about, [('out_report', 'in_file')]),
    ])  # fmt:skip

    # template_iterator_wf = None
    # select_MNI2009c_xfm = None
    if config.workflow.level == 'full':
        # Much of the logic here is extracted into a separate, fuller anatomical workflow
        # TODO:
        # - Grab template_iterator_wf workflow
        # - Grab select_MNI2009c_xfm node
        pass

        # if 'MNI152NLin2009cAsym' in spaces.get_spaces():
        #     select_MNI2009c_xfm = pe.Node(
        #         KeySelect(fields=['std2anat_xfm'], key='MNI152NLin2009cAsym'),
        #         name='select_MNI2009c_xfm',
        #         run_without_submitting=True,
        #     )
        #     workflow.connect([
        #         (anat_fit_wf, select_MNI2009c_xfm, [
        #             ('outputnode.std2anat_xfm', 'std2anat_xfm'),
        #             ('outputnode.template', 'keys'),
        #         ]),
        #     ])  # fmt:skip

        # Thread MNI152NLin6Asym standard outputs to CIFTI subworkflow, skipping
        # the iterator, which targets only output spaces.
        # This can lead to duplication in the working directory if people actually
        # want MNI152NLin6Asym outputs, but we'll live with it.
        # if config.workflow.cifti_output:
        #     from smriprep.interfaces.templateflow import TemplateFlowSelect

        #     ref = Reference(
        #         'MNI152NLin6Asym',
        #         {'res': 2 if config.workflow.cifti_output == '91k' else 1},
        #     )

        #     select_MNI6_xfm = pe.Node(
        #         KeySelect(fields=['anat2std_xfm'], key=ref.fullname),
        #         name='select_MNI6',
        #         run_without_submitting=True,
        #     )
        #     select_MNI6_tpl = pe.Node(
        #         TemplateFlowSelect(template=ref.fullname, resolution=ref.spec['res']),
        #         name='select_MNI6_tpl',
        #     )
        #     workflow.connect([
        #         (anat_fit_wf, select_MNI6_xfm, [
        #             ('outputnode.anat2std_xfm', 'anat2std_xfm'),
        #             ('outputnode.template', 'keys'),
        #         ]),
        #     ])  # fmt:skip

    if config.workflow.anat_only:
        return clean_datasinks(workflow)

    # TODO: FMAP, BOLD PROCESSING
    return workflow

    # # fmt: off
    # workflow.connect([
    #     (inputnode, anat_preproc_wf, [
    #         ('subjects_dir', 'inputnode.subjects_dir'),
    #     ]),
    #     (inputnode, summary, [
    #         ('subjects_dir', 'subjects_dir'),
    #     ]),
    #     (bidssrc, summary, [
    #         ('bold', 'bold'),
    #     ]),
    #     (bids_info, summary, [
    #         ('subject', 'subject_id'),
    #     ]),
    #     (bids_info, anat_preproc_wf, [
    #         (('subject', _prefix), 'inputnode.subject_id'),
    #     ]),
    #     (bidssrc, anat_preproc_wf, [
    #         ('t1w', 'inputnode.t1w'),
    #         ('t2w', 'inputnode.t2w'),
    #     ]),
    #     (summary, ds_report_summary, [
    #         ('out_report', 'in_file'),
    #     ]),
    #     (about, ds_report_about, [
    #         ('out_report', 'in_file'),
    #     ]),
    # ])

    # workflow.connect([
    #     (bidssrc, bids_info, [
    #         ((contrast.lower(), fix_multi_source_name), 'in_file'),
    #     ]),
    #     (bidssrc, summary, [
    #         ('t1w', 't1w'),
    #         ('t2w', 't2w'),
    #     ]),
    #     (bidssrc, ds_report_summary, [
    #         ((contrast.lower(), fix_multi_source_name), 'source_file'),
    #     ]),
    #     (bidssrc, ds_report_about, [
    #         ((contrast.lower(), fix_multi_source_name), 'source_file'),
    #     ]),
    # ])
    # # fmt: on

    # # Overwrite ``out_path_base`` of smriprep's DataSinks
    # for node in workflow.list_node_names():
    #     if node.split('.')[-1].startswith('ds_'):
    #         workflow.get_node(node).interface.out_path_base = ''

    # if anat_only:
    #     return workflow

    # Susceptibility distortion correction


#     fmap_estimators = None
#     if any((config.workflow.use_syn_sdc, config.workflow.force_syn)):
#         config.loggers.workflow.critical('SyN processing is not yet implemented.')

#     if 'fieldmaps' not in config.workflow.ignore:
#         from sdcflows.utils.wrangler import find_estimators

#         # SDC Step 1: Run basic heuristics to identify available data for fieldmap estimation
#         # For now, no fmapless
#         fmap_estimators = find_estimators(
#             layout=config.execution.layout,
#             subject=subject_id,
#             sessions=[session_id],
#             fmapless=False,  # config.workflow.use_syn,
#             force_fmapless=False,  # config.workflow.force_syn,
#         )

#     # Append the functional section to the existing anatomical exerpt
#     # That way we do not need to stream down the number of bold datasets
#     anat_preproc_wf.__postdesc__ = anat_preproc_wf.__postdesc__ or ''
#     func_pre_desc = f"""

# Functional data preprocessing

# : For each of the {len(subject_data['bold'])} BOLD runs found per subject (across all
# tasks and sessions), the following preprocessing was performed."""

#     func_preproc_wfs = []
#     has_fieldmap = bool(fmap_estimators)
#     for bold_file in subject_data['bold']:
#         func_preproc_wf = init_func_preproc_wf(bold_file, spaces, has_fieldmap=has_fieldmap)
#         if func_preproc_wf is None:
#             continue

#         func_preproc_wf.__desc__ = func_pre_desc + (func_preproc_wf.__desc__ or '')
#         # fmt:off
#         workflow.connect([
#             (anat_preproc_wf, func_preproc_wf, [
#                 ('outputnode.anat_preproc', 'inputnode.anat_preproc'),
#                 ('outputnode.anat_mask', 'inputnode.anat_mask'),
#                 ('outputnode.anat_brain', 'inputnode.anat_brain'),
#                 ('outputnode.anat_dseg', 'inputnode.anat_dseg'),
#                 ('outputnode.anat_aseg', 'inputnode.anat_aseg'),
#                 ('outputnode.anat_aparc', 'inputnode.anat_aparc'),
#                 ('outputnode.anat_tpms', 'inputnode.anat_tpms'),
#                 ('outputnode.template', 'inputnode.template'),
#                 ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
#                 ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
#                 # Undefined if --fs-no-reconall, but this is safe
#                 ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
#                 ('outputnode.subject_id', 'inputnode.subject_id'),
#                 ('outputnode.anat2fsnative_xfm', 'inputnode.t1w2fsnative_xfm'),
#                 ('outputnode.fsnative2anat_xfm', 'inputnode.fsnative2t1w_xfm'),
#                 ('outputnode.surfaces', 'inputnode.surfaces'),
#                 ('outputnode.morphometrics', 'inputnode.morphometrics'),
#                 ('outputnode.anat_ribbon', 'inputnode.anat_ribbon'),
#                 ('outputnode.sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR'),
#                 ('outputnode.midthickness_fsLR', 'inputnode.midthickness_fsLR'),
#             ]),
#         ])
#         # fmt:on
#         func_preproc_wfs.append(func_preproc_wf)

#     if not has_fieldmap:
#         config.loggers.workflow.warning(
#             'Data for fieldmap estimation not present. Please note that these data '
#             'will not be corrected for susceptibility distortions.'
#         )
#         return workflow

#     config.loggers.workflow.info(
#         f'Fieldmap estimators found: {[e.method for e in fmap_estimators]}'
#     )

#     from sdcflows import fieldmaps as fm
#     from sdcflows.workflows.base import init_fmap_preproc_wf

#     fmap_wf = init_fmap_preproc_wf(
#         sloppy=bool(config.execution.sloppy),
#         debug='fieldmaps' in config.execution.debug,
#         estimators=fmap_estimators,
#         omp_nthreads=config.nipype.omp_nthreads,
#         output_dir=nibabies_dir,
#         subject=subject_id,
#     )
#     fmap_wf.__desc__ = f"""

# Preprocessing of B<sub>0</sub> inhomogeneity mappings

# : A total of {len(fmap_estimators)} fieldmaps were found available within the input
# BIDS structure for this particular subject.
# """

#     for func_preproc_wf in func_preproc_wfs:
#         # fmt: off
#         workflow.connect([
#             (fmap_wf, func_preproc_wf, [
#                 ('outputnode.fmap', 'inputnode.fmap'),
#                 ('outputnode.fmap_ref', 'inputnode.fmap_ref'),
#                 ('outputnode.fmap_coeff', 'inputnode.fmap_coeff'),
#                 ('outputnode.fmap_mask', 'inputnode.fmap_mask'),
#                 ('outputnode.fmap_id', 'inputnode.fmap_id'),
#                 ('outputnode.method', 'inputnode.sdc_method'),
#             ]),
#         ])
#         # fmt: on

#     # Overwrite ``out_path_base`` of sdcflows's DataSinks
#     for node in fmap_wf.list_node_names():
#         if node.split('.')[-1].startswith('ds_'):
#             fmap_wf.get_node(node).interface.out_path_base = ''

#     # Step 3: Manually connect PEPOLAR
#     for estimator in fmap_estimators:
#         config.loggers.workflow.info(
#             f"""\
# Setting-up fieldmap "{estimator.bids_id}" ({estimator.method}) with \
# <{', '.join(s.path.name for s in estimator.sources)}>"""
#         )
#         if estimator.method in (fm.EstimatorType.MAPPED, fm.EstimatorType.PHASEDIFF):
#             continue

#         suffices = [s.suffix for s in estimator.sources]

#         if estimator.method == fm.EstimatorType.PEPOLAR:
#             if set(suffices) == {'epi'} or sorted(suffices) == ['bold', 'epi']:
#                 fmap_wf_inputs = getattr(fmap_wf.inputs, f'in_{estimator.bids_id}')
#                 fmap_wf_inputs.in_data = [str(s.path) for s in estimator.sources]
#                 fmap_wf_inputs.metadata = [s.metadata for s in estimator.sources]
#             else:
#                 raise NotImplementedError(
#                     'Sophisticated PEPOLAR schemes (e.g., using DWI+EPI) are unsupported.'
#                 )

#     return workflow


def _subject_session_id(subject_id: str, session_id: str | None) -> str:
    """
    Combine a subject ID with a session ID (if available).

    >>> _subject_session_id('01', None)
    'sub-01'
    >>> _subject_session_id('sub-01', '01')
    'sub-01_ses-01'
    >>> _subject_session_id('01', 'ses-03')
    'sub-01_ses-03'
    """
    entities = []
    entities.append(f'sub-{subject_id}' if not subject_id.startswith('sub-') else subject_id)
    if session_id is not None:
        entities.append(f'ses-{session_id}' if not session_id.startswith('ses-') else session_id)
    return '_'.join(entities)


def clean_datasinks(workflow: pe.Workflow) -> pe.Workflow:
    # Overwrite ``out_path_base`` of smriprep's DataSinks
    for node in workflow.list_node_names():
        if node.split('.')[-1].startswith('ds_'):
            workflow.get_node(node).interface.out_path_base = ''
    return workflow


def init_workflow_spaces(execution_spaces: SpatialReferences, age_months: int):
    """
    Create output spaces at a per-subworkflow level.

    This address the case where a multi-session subject is run,
    and requires separate template cohorts.
    """
    from niworkflows.utils.spaces import Reference

    from nibabies.utils.misc import cohort_by_months

    spaces = deepcopy(execution_spaces)

    if age_months is None:
        raise RuntimeError('Participant age (in months) is required.')

    if not spaces.references:
        # Ensure age specific template is added if nothing is present
        cohort = cohort_by_months('MNIInfant', age_months)
        spaces.add(('MNIInfant', {'res': 'native', 'cohort': cohort}))

    if not spaces.is_cached():
        spaces.checkpoint()

    # Ensure user-defined spatial references for outputs are correctly parsed.
    # Certain options require normalization to a space not explicitly defined by users.
    # These spaces will not be included in the final outputs.
    if config.workflow.use_aroma:
        # Make sure there's a normalization to FSL for AROMA to use.
        spaces.add(Reference('MNI152NLin6Asym', {'res': '2'}))

    if config.workflow.cifti_output:
        # CIFTI grayordinates to corresponding FSL-MNI resolutions.
        vol_res = '2' if config.workflow.cifti_output == '91k' else '1'
        spaces.add(Reference('MNI152NLin6Asym', {'res': vol_res}))
        # Ensure a non-native version of MNIInfant is added as a target
        cohort = cohort_by_months('MNIInfant', age_months)
        spaces.add(Reference('MNIInfant', {'cohort': cohort}))

    return spaces


def init_execution_spaces():
    from niworkflows.utils.spaces import Reference, SpatialReferences

    spaces = config.execution.output_spaces or SpatialReferences()
    if not isinstance(spaces, SpatialReferences):
        spaces = SpatialReferences(
            [ref for s in spaces.split(' ') for ref in Reference.from_string(s)]
        )
    return spaces


def map_fieldmap_estimation(
    layout: BIDSLayout,
    subject_id: str,
    bold_data: list[list[str]],
    ignore_fieldmaps: bool,
    use_syn: bool | str,
    force_syn: bool,
    filters: dict | None,
) -> tuple[list, dict]:
    if not any((not ignore_fieldmaps, use_syn, force_syn)):
        return [], {}

    from sdcflows import fieldmaps as fm
    from sdcflows.utils.wrangler import find_estimators

    # In the case where fieldmaps are ignored and `--use-syn-sdc` is requested,
    # SDCFlows `find_estimators` still receives a full layout (which includes the fmap modality)
    # and will not calculate fmapless schemes.
    # Similarly, if fieldmaps are ignored and `--force-syn` is requested,
    # `fmapless` should be set to True to ensure BOLD targets are found to be corrected.
    fmap_estimators = find_estimators(
        layout=layout,
        subject=subject_id,
        fmapless=bool(use_syn) or ignore_fieldmaps and force_syn,
        force_fmapless=force_syn or ignore_fieldmaps and use_syn,
        bids_filters=filters,
    )

    if not fmap_estimators:
        if use_syn:
            message = (
                'Fieldmap-less (SyN) estimation was requested, but PhaseEncodingDirection '
                'information appears to be absent.'
            )
            config.loggers.workflow.error(message)
            if use_syn == 'error':
                raise ValueError(message)
        return [], {}

    if ignore_fieldmaps and any(f.method == fm.EstimatorType.ANAT for f in fmap_estimators):
        config.loggers.workflow.info(
            'Option "--ignore fieldmaps" was set, but either "--use-syn-sdc" '
            'or "--force-syn" were given, so fieldmap-less estimation will be executed.'
        )
        fmap_estimators = [f for f in fmap_estimators if f.method == fm.EstimatorType.ANAT]

    # Pare down estimators to those that are actually used
    # If fmap_estimators == [], all loops/comprehensions terminate immediately
    all_ids = {fmap.bids_id for fmap in fmap_estimators}
    bold_files = (bold_series[0] for bold_series in bold_data)

    all_estimators = {
        bold_file: [fmap_id for fmap_id in get_estimator(layout, bold_file) if fmap_id in all_ids]
        for bold_file in bold_files
    }

    for bold_file, estimator_key in all_estimators.items():
        if len(estimator_key) > 1:
            config.loggers.workflow.warning(
                f"Several fieldmaps <{', '.join(estimator_key)}> are "
                f"'IntendedFor' <{bold_file}>, using {estimator_key[0]}"
            )
            estimator_key[1:] = []

    # Final, 1-1 map, dropping uncorrected BOLD
    estimator_map = {
        bold_file: estimator_key[0]
        for bold_file, estimator_key in all_estimators.items()
        if estimator_key
    }

    fmap_estimators = [f for f in fmap_estimators if f.bids_id in estimator_map.values()]

    return fmap_estimators, estimator_map


def get_estimator(layout, fname):
    field_source = layout.get_metadata(fname).get('B0FieldSource')
    if isinstance(field_source, str):
        field_source = (field_source,)

    if field_source is None:
        import re
        from pathlib import Path

        from sdcflows.fieldmaps import get_identifier

        # Fallback to IntendedFor
        intended_rel = re.sub(r'^sub-[a-zA-Z0-9]*/', '', str(Path(fname).relative_to(layout.root)))
        field_source = get_identifier(intended_rel)

    return field_source
