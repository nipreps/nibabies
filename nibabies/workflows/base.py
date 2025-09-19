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
NiBabies base processing workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_nibabies_wf
.. autofunction:: init_single_subject_wf

"""

from __future__ import annotations

import os
import pprint
import sys
import typing as ty
import warnings
from copy import deepcopy

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.interfaces.utility import KeySelect
from niworkflows.utils.connections import listify
from packaging.version import Version
from smriprep.workflows.outputs import init_template_iterator_wf

from nibabies import config
from nibabies.interfaces import DerivativesDataSink
from nibabies.interfaces.reports import AboutSummary, SubjectSummary
from nibabies.utils.bids import parse_bids_for_age_months
from nibabies.workflows.anatomical.apply import init_infant_anat_apply_wf
from nibabies.workflows.anatomical.fit import (
    init_infant_anat_fit_wf,
    init_infant_single_anat_fit_wf,
)
from nibabies.workflows.bold.base import init_bold_wf

LOGGER = config.loggers.workflow

if ty.TYPE_CHECKING:
    from bids.layout import BIDSLayout
    from niworkflows.utils.spaces import SpatialReferences

    SubjectSession = tuple[str, str | None]


MCRIBS_RECOMMEND_AGE_CAP = 3


def init_nibabies_wf(subworkflows_list: list[SubjectSession]):
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
                wf = init_nibabies_wf(['01', None])

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
    surface_recon = config.workflow.surface_recon_method is not None
    if surface_recon:
        fsdir = pe.Node(
            BIDSFreeSurferDir(
                derivatives=config.execution.output_dir,
                freesurfer_home=os.getenv('FREESURFER_HOME'),
                spaces=execution_spaces.get_fs_spaces(),
            ),
            name=f'fsdir_run_{config.execution.run_uuid.replace("-", "_")}',
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
        if surface_recon:
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

    bold_runs = [
        sorted(
            listify(run),
            key=lambda fl: config.execution.layout.get_metadata(fl).get('EchoTime', 0),
        )
        for run in subject_data['bold']
    ]

    if subject_data['roi']:
        warnings.warn(
            f'Lesion mask {subject_data["roi"]} found. '
            'Future versions of NiBabies will use alternative conventions. '
            'Please refer to the documentation before upgrading.',
            FutureWarning,
            stacklevel=1,
        )

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

        if config.execution.copy_derivatives:
            from nibabies.utils.derivatives import copy_derivatives

            LOGGER.info('Copying found anat derivatives into output directory')
            copy_derivatives(
                derivs=anatomical_cache,
                outdir=config.execution.nibabies_dir,
                modality='anat',
                subject_id=f'sub-{subject_id}',
                session_id=f'ses-{session_id}' if session_id else None,
                config_hash=config.execution.parameters_hash
                if config.execution.output_layout == 'multiverse'
                else None,
            )

    # Determine some session level options here, as we should have
    # all the required information
    if recon_method == 'auto':
        if age <= MCRIBS_RECOMMEND_AGE_CAP and anatomical_cache.get('t2w_aseg'):
            # do not force mcribs without a vetted segmentation
            recon_method = 'mcribs'
        elif age <= 24:
            recon_method = 'infantfs'
        else:
            recon_method = 'freesurfer'

    requested_anat = config.execution.reference_anat
    t1w = subject_data['t1w']
    t2w = subject_data['t2w']
    single_anat = False

    if not (t1w and t2w):
        single_anat = True
        reference_anat = 'T1w' if t1w else 'T2w'
        if requested_anat and reference_anat != requested_anat:
            raise AttributeError(
                f'Requested to use {requested_anat} as anatomical reference but none available'
            )
    elif (reference_anat := requested_anat) is None:  # Both available with no preference
        reference_anat = (
            'T2w'
            if any(
                (
                    recon_method is None and age <= MCRIBS_RECOMMEND_AGE_CAP,
                    recon_method == 'mcribs',
                )
            )
            else 'T1w'
        )

    anat = reference_anat.lower()  # To be used for workflow connections
    LOGGER.info(
        'Collected the following data for %s:\nRaw:\n%s\n\nDerivatives:\n\n%s\n',
        f'sub-{subject_id}' if not session_id else f'sub-{subject_id}_ses-{session_id}',
        pprint.pformat(subject_data),
        pprint.pformat(anatomical_cache),
    )

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
            require_t1w=False,
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
            age=age,
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

    output_dir = config.execution.nibabies_dir
    sloppy = config.execution.sloppy
    cifti_output = config.workflow.cifti_output

    wf_args = {
        'age_months': age,
        't1w': t1w,
        't2w': t2w,
        'flair': subject_data['flair'],
        'bids_root': bids_root,
        'longitudinal': config.workflow.longitudinal,
        'msm_sulc': msm_sulc,
        'omp_nthreads': omp_nthreads,
        'output_dir': output_dir,
        'precomputed': anatomical_cache,
        'segmentation_atlases': config.execution.segmentation_atlases_dir,
        'skull_strip_fixed_seed': config.workflow.skull_strip_fixed_seed,
        'skull_strip_mode': config.workflow.skull_strip_anat,
        'skull_strip_template': Reference.from_string(config.workflow.skull_strip_template)[0],
        'recon_method': recon_method,
        'reference_anat': reference_anat,
        'sloppy': sloppy,
        'spaces': spaces,
        'cifti_output': cifti_output,
    }

    fit_wf = init_infant_single_anat_fit_wf if single_anat else init_infant_anat_fit_wf
    anat_fit_wf = fit_wf(**wf_args)

    # allow to run with anat-fast-track on fMRI-only dataset
    if f'{anat}_preproc' in anatomical_cache and not subject_data[anat]:
        workflow.connect([
            (bidssrc, bids_info, [(('bold', fix_multi_source_name), 'in_file')]),
            (anat_fit_wf, summary, [('outputnode.anat_preproc', anat)]),
            (anat_fit_wf, ds_report_summary, [('outputnode.anat_preproc', 'source_file')]),
            (anat_fit_wf, ds_report_about, [('outputnode.anat_preproc', 'source_file')]),
        ])  # fmt:skip
    else:
        workflow.connect([
            (bidssrc, bids_info, [((anat, fix_multi_source_name), 'in_file')]),
            (bidssrc, summary, [('t1w', 't1w')]),
            (bidssrc, ds_report_summary, [((anat, fix_multi_source_name), 'source_file')]),
            (bidssrc, ds_report_about, [((anat, fix_multi_source_name), 'source_file')]),
        ])  # fmt:skip

    if single_anat:
        workflow.connect(bidssrc, anat, anat_fit_wf, 'inputnode.anat')

    workflow.connect([
        (inputnode, anat_fit_wf, [('subjects_dir', 'inputnode.subjects_dir')]),
        (bidssrc, anat_fit_wf, [
            ('t1w', 'inputnode.t1w'),
            ('t2w', 'inputnode.t2w'),
            ('roi', 'inputnode.roi'),
            ('flair', 'inputnode.flair'),
        ]),
        # Reporting connections
        (inputnode, summary, [('subjects_dir', 'subjects_dir')]),
        (bidssrc, summary, [('t2w', 't2w'), ('bold', 'bold')]),
        (bids_info, summary, [
            ('subject', 'subject_id'),
            ('session', 'session_id'),
        ]),
        (summary, ds_report_summary, [('out_report', 'in_file')]),
        (summary, anat_fit_wf, [('subject_id', 'inputnode.subject_id')]),
        (about, ds_report_about, [('out_report', 'in_file')]),
    ])  # fmt:skip

    reg_sphere = f'sphere_reg_{"msm" if msm_sulc else "fsLR"}'
    template_iterator_wf = None
    select_MNIInfant_xfm = None
    if config.workflow.level == 'full':
        anat_apply_wf = init_infant_anat_apply_wf(
            bids_root=bids_root,
            cifti_output=cifti_output,
            msm_sulc=msm_sulc,
            omp_nthreads=omp_nthreads,
            output_dir=output_dir,
            precomputed=anatomical_cache,
            recon_method=recon_method,
            reference_anat=reference_anat,
            spaces=spaces,
        )
        workflow.connect([
            (anat_fit_wf, anat_apply_wf, [
                ('outputnode.anat_valid_list', 'inputnode.anat_valid_list'),
                ('outputnode.anat_preproc', 'inputnode.anat_preproc'),
                ('outputnode.anat_mask', 'inputnode.anat_mask'),
                ('outputnode.anat_dseg', 'inputnode.anat_dseg'),
                ('outputnode.anat_tpms', 'inputnode.anat_tpms'),
                ('outputnode.fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
                ('outputnode.white', 'inputnode.white'),
                ('outputnode.pial', 'inputnode.pial'),
                ('outputnode.midthickness', 'inputnode.midthickness'),
                ('outputnode.cortex_mask', 'inputnode.cortex_mask'),
                (f'outputnode.{reg_sphere}', f'inputnode.{reg_sphere}'),
                ('outputnode.sulc', 'inputnode.sulc'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.thickness', 'inputnode.thickness'),
            ]),
        ])  # fmt:skip

        if cifti_output and 'MNIInfant' in [ref.space for ref in spaces.references]:
            mniinfant_res = 2 if config.workflow.cifti_output == '91k' else 1

            select_MNIInfant_xfm = pe.Node(
                KeySelect(
                    fields=['anat2std_xfm', 'std2anat_xfm'],
                    key=get_MNIInfant_key(spaces, mniinfant_res),
                ),
                name='select_MNIInfant_xfm',
                run_without_submitting=True,
            )

            workflow.connect([
                (anat_fit_wf, select_MNIInfant_xfm, [
                    ('outputnode.std2anat_xfm', 'std2anat_xfm'),
                    ('outputnode.anat2std_xfm', 'anat2std_xfm'),
                    ('outputnode.template', 'keys'),
                ]),
            ])  # fmt:skip

        if spaces.cached.get_spaces(nonstandard=False, dim=(3,)):
            template_iterator_wf = init_template_iterator_wf(spaces=spaces, sloppy=sloppy)

            workflow.connect([
                (anat_fit_wf, template_iterator_wf, [
                    ('outputnode.template', 'inputnode.template'),
                    ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
                ]),
                (template_iterator_wf, anat_apply_wf, [
                    ('outputnode.std_t1w', 'inputnode.std_t1w',),
                    ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
                    ('outputnode.space', 'inputnode.std_space'),
                    ('outputnode.cohort', 'inputnode.std_cohort'),
                    ('outputnode.resolution', 'inputnode.std_resolution'),
                ]),
            ])  # fmt:skip

    if config.workflow.anat_only:
        return clean_datasinks(workflow)

    fmap_estimators, estimator_map = map_fieldmap_estimation(
        layout=config.execution.layout,
        subject_id=subject_id,
        bold_data=bold_runs,
        ignore_fieldmaps='fieldmaps' in config.workflow.ignore,
        use_syn=config.workflow.use_syn_sdc,
        force_syn=config.workflow.force_syn,
        filters=config.execution.get().get('bids_filters', {}).get('fmap'),
    )

    if fmap_estimators:
        LOGGER.info(
            'B0 field inhomogeneity map will be estimated with the following '
            f'{len(fmap_estimators)} estimator(s): '
            f'{[e.method for e in fmap_estimators]}.'
        )

        from sdcflows import fieldmaps as fm
        from sdcflows.workflows.base import init_fmap_preproc_wf

        fmap_wf = init_fmap_preproc_wf(
            debug='fieldmaps' in config.execution.debug,
            estimators=fmap_estimators,
            omp_nthreads=omp_nthreads,
            output_dir=output_dir,
            subject=subject_id,
            sd_prior=False,  # No priors for infants yet
        )
        fmap_wf.__desc__ = f"""

Preprocessing of B<sub>0</sub> inhomogeneity mappings

: A total of {len(fmap_estimators)} fieldmaps were found available within the input
BIDS structure for this particular subject.
"""

        # Overwrite ``out_path_base`` of sdcflows's DataSinks
        for node in fmap_wf.list_node_names():
            if node.split('.')[-1].startswith('ds_'):
                fmap_wf.get_node(node).interface.out_path_base = ''

        # MG: No prior is used ATM, so no need for xfm
        # fmap_select_std = pe.Node(
        #     KeySelect(fields=['std2anat_xfm'], key='MNI152NLin2009cAsym'),
        #     name='fmap_select_std',
        #     run_without_submitting=True,
        # )
        # if any(estimator.method == fm.EstimatorType.ANAT for estimator in fmap_estimators):
        #     workflow.connect([
        #         (anat_fit_wf, fmap_select_std, [
        #             ('outputnode.std2anat_xfm', 'std2anat_xfm'),
        #             ('outputnode.template', 'keys')]),
        #     ])  # fmt:skip

        for estimator in fmap_estimators:
            LOGGER.info(
                f"""\
Setting-up fieldmap "{estimator.bids_id}" ({estimator.method}) with \
<{', '.join(s.path.name for s in estimator.sources)}>"""
            )

            # Mapped and phasediff can be connected internally by SDCFlows
            if estimator.method in (fm.EstimatorType.MAPPED, fm.EstimatorType.PHASEDIFF):
                continue

            suffices = [s.suffix for s in estimator.sources]

            if estimator.method == fm.EstimatorType.PEPOLAR:
                if len(suffices) == 2 and all(suf in ('epi', 'bold', 'sbref') for suf in suffices):
                    wf_inputs = getattr(fmap_wf.inputs, f'in_{estimator.bids_id}')
                    wf_inputs.in_data = [str(s.path) for s in estimator.sources]
                    wf_inputs.metadata = [s.metadata for s in estimator.sources]
                else:
                    raise NotImplementedError('Sophisticated PEPOLAR schemes are unsupported.')

            elif estimator.method == fm.EstimatorType.ANAT:
                from sdcflows.workflows.fit.syn import init_syn_preprocessing_wf

                sources = [str(s.path) for s in estimator.sources if s.suffix in ('bold', 'sbref')]
                source_meta = [
                    s.metadata for s in estimator.sources if s.suffix in ('bold', 'sbref')
                ]
                syn_preprocessing_wf = init_syn_preprocessing_wf(
                    omp_nthreads=omp_nthreads,
                    debug=config.execution.sloppy,
                    auto_bold_nss=True,
                    t1w_inversion=False,
                    sd_prior=False,
                    name=f'syn_preprocessing_{estimator.bids_id}',
                )
                syn_preprocessing_wf.inputs.inputnode.in_epis = sources
                syn_preprocessing_wf.inputs.inputnode.in_meta = source_meta

                workflow.connect([
                    (anat_fit_wf, syn_preprocessing_wf, [
                        ('outputnode.anat_preproc', 'inputnode.in_anat'),
                        ('outputnode.anat_mask', 'inputnode.mask_anat'),
                    ]),
                    # MG: No prior is used ATM, so no need for xfm
                    # (fmap_select_std, syn_preprocessing_wf, [
                    #     ('std2anat_xfm', 'inputnode.std2anat_xfm'),
                    # ]),
                    (syn_preprocessing_wf, fmap_wf, [
                        ('outputnode.epi_ref', f'in_{estimator.bids_id}.epi_ref'),
                        ('outputnode.epi_mask', f'in_{estimator.bids_id}.epi_mask'),
                        ('outputnode.anat_ref', f'in_{estimator.bids_id}.anat_ref'),
                        ('outputnode.anat_mask', f'in_{estimator.bids_id}.anat_mask'),
                        ('outputnode.sd_prior', f'in_{estimator.bids_id}.sd_prior'),
                    ]),
                ])  # fmt:skip

    # Append the functional section to the existing anatomical excerpt
    # That way we do not need to stream down the number of bold datasets
    func_pre_desc = f"""
Functional data preprocessing

: For each of the {len(bold_runs)} BOLD runs found per subject (across all
tasks and sessions), the following preprocessing was performed.
"""

    # Before initializing BOLD workflow, select/verify anatomical target for coregistration
    if config.workflow.bold2anat_init in ('auto', 't2w'):
        has_t2w = subject_data['t2w'] or 't2w_preproc' in anatomical_cache
        if config.workflow.bold2anat_init == 't2w' and not has_t2w:
            raise OSError(
                'A T2w image is expected for BOLD-to-anatomical coregistration and was not found'
            )
        config.workflow.bold2anat_init = 't2w' if has_t2w else 't1w'

    for bold_series in bold_runs:
        bold_file = bold_series[0]
        fieldmap_id = estimator_map.get(bold_file)

        functional_cache = {}
        if config.execution.derivatives:
            from nibabies.utils.bids import extract_entities
            from nibabies.utils.derivatives import collect_functional_derivatives

            entities = extract_entities(bold_series)

            for deriv_dir in config.execution.derivatives.values():
                functional_cache.update(
                    collect_functional_derivatives(
                        derivatives_dir=deriv_dir,
                        entities=entities,
                        fieldmap_id=fieldmap_id,
                    )
                )

            if config.execution.copy_derivatives:
                from nibabies.utils.derivatives import copy_derivatives

                LOGGER.info('Copying found func derivatives into output directory')
                copy_derivatives(
                    derivs=functional_cache,
                    outdir=config.execution.nibabies_dir,
                    modality='func',
                    subject_id=f'sub-{subject_id}',
                    session_id=f'ses-{session_id}' if session_id else None,
                    config_hash=config.execution.parameters_hash
                    if config.execution.output_layout == 'multiverse'
                    else None,
                )

        bold_wf = init_bold_wf(
            reference_anat=reference_anat,
            bold_series=bold_series,
            precomputed=functional_cache,
            fieldmap_id=fieldmap_id,
            spaces=spaces,
        )
        if bold_wf is None:
            continue

        bold_wf.__desc__ = func_pre_desc + (bold_wf.__desc__ or '')

        workflow.connect([
            (anat_fit_wf, bold_wf, [
                ('outputnode.anat_preproc', 'inputnode.anat_preproc'),
                ('outputnode.anat_mask', 'inputnode.anat_mask'),
                ('outputnode.anat_dseg', 'inputnode.anat_dseg'),
                ('outputnode.anat_tpms', 'inputnode.anat_tpms'),
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
                ('outputnode.white', 'inputnode.white'),
                ('outputnode.pial', 'inputnode.pial'),
                ('outputnode.midthickness', 'inputnode.midthickness'),
                ('outputnode.anat_ribbon', 'inputnode.anat_ribbon'),
                (f'outputnode.{reg_sphere}', 'inputnode.sphere_reg_fsLR'),
            ]),
        ])  # fmt:skip
        if fieldmap_id:
            workflow.connect([
                (fmap_wf, bold_wf, [
                    ('outputnode.fmap', 'inputnode.fmap'),
                    ('outputnode.fmap_ref', 'inputnode.fmap_ref'),
                    ('outputnode.fmap_coeff', 'inputnode.fmap_coeff'),
                    ('outputnode.fmap_mask', 'inputnode.fmap_mask'),
                    ('outputnode.fmap_id', 'inputnode.fmap_id'),
                    ('outputnode.method', 'inputnode.sdc_method'),
                ]),
            ])  # fmt:skip

        if config.workflow.level == 'full':
            if template_iterator_wf is not None:
                workflow.connect([
                    (template_iterator_wf, bold_wf, [
                        ('outputnode.space', 'inputnode.std_space'),
                        ('outputnode.resolution', 'inputnode.std_resolution'),
                        ('outputnode.cohort', 'inputnode.std_cohort'),
                        ('outputnode.std_t1w', 'inputnode.std_t1w'),
                        ('outputnode.std_mask', 'inputnode.std_mask'),
                        ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
                    ]),
                ])  # fmt:skip

            if select_MNIInfant_xfm is not None:
                workflow.connect([
                    (select_MNIInfant_xfm, bold_wf, [
                        ('std2anat_xfm', 'inputnode.mniinfant2anat_xfm'),
                        ('anat2std_xfm', 'inputnode.anat2mniinfant_xfm')
                    ]),
                ])  # fmt:skip

            if config.workflow.cifti_output:
                workflow.connect([
                    (anat_fit_wf, bold_wf, [
                        ('outputnode.cortex_mask', 'inputnode.cortex_mask'),
                    ]),
                    (anat_apply_wf, bold_wf, [
                        ('outputnode.midthickness_fsLR', 'inputnode.midthickness_fsLR'),
                        ('outputnode.anat_aseg', 'inputnode.anat_aseg'),
                    ]),
                ])  # fmt:skip

    return clean_datasinks(workflow)


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
            workflow.get_node(node).interface.inputs.base_directory = config.execution.nibabies_dir

            if config.execution.output_layout == 'multiverse':
                workflow.get_node(node).interface.inputs.hash = config.execution.parameters_hash
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
        LOGGER.debug('No references specified, MNIInfant:cohort-%s as default', cohort)

    if not spaces.is_cached():
        spaces.checkpoint()

    # Ensure one cohort of MNIInfant is always available as an internal space
    if not any(
        space.startswith('MNIInfant') for space in spaces.get_spaces(nonstandard=False, dim=(3,))
    ):
        cohort = cohort_by_months('MNIInfant', age_months)
        spaces.add(Reference('MNIInfant', {'cohort': cohort}))
        LOGGER.debug('Missing internal space, adding MNIInfant:cohort-%s', cohort)

    if config.workflow.cifti_output:
        # CIFTI grayordinates to corresponding FSL-MNI resolutions.
        vol_res = '2' if config.workflow.cifti_output == '91k' else '1'
        spaces.add(Reference('MNI152NLin6Asym', {'res': vol_res}))
        # Ensure a non-native version of MNIInfant is added as a target
        cohort = cohort_by_months('MNIInfant', age_months)
        spaces.add(Reference('MNIInfant', {'cohort': cohort, 'res': vol_res}))
        LOGGER.debug(
            'Adding MNI152NLin6Asym:res-%s, MNIInfant:cohort-%s:res-%s',
            vol_res,
            cohort,
            vol_res,
        )

    LOGGER.debug('Workflow spaces: %s', spaces.get_spaces())
    if not any(
        space.startswith('MNIInfant') for space in spaces.get_spaces(nonstandard=False, dim=(3,))
    ):
        raise RuntimeError(
            'MNIInfant space is required but not found, likely due to a stale templateflow cache. '
            'Clear the cache (default: $HOME/.cache/templateflow) and retry.'
        )
    return spaces


def init_execution_spaces():
    """Initialize the spaces to be saved.

    Either invoked by ``--output-spaces``,
    or an empty :py:class:`~niworkflows.utils.spaces.SpatialReferences`
    """
    from niworkflows.utils.spaces import Reference, SpatialReferences

    spaces = config.execution.output_spaces or SpatialReferences()
    if not isinstance(spaces, SpatialReferences):
        spaces = SpatialReferences(
            [ref for s in spaces.split(' ') for ref in Reference.from_string(s)]
        )
    LOGGER.debug('Execution spaces: %s', spaces.get_spaces())
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
        anat_suffix=['T1w', 'T2w'],
    )

    if not fmap_estimators:
        if use_syn:
            message = (
                'Fieldmap-less (SyN) estimation was requested, but PhaseEncodingDirection '
                'information appears to be absent.'
            )
            LOGGER.error(message)
            if use_syn == 'error':
                raise ValueError(message)
        return [], {}

    if ignore_fieldmaps:
        if any(f.method == fm.EstimatorType.ANAT for f in fmap_estimators):
            LOGGER.info(
                'Option "--ignore fieldmaps" was set, but either "--use-syn-sdc" '
                'or "--force-syn" were given, so fieldmap-less estimation will be executed.'
            )
            fmap_estimators = [f for f in fmap_estimators if f.method == fm.EstimatorType.ANAT]
        else:
            LOGGER.info('Ignoring fieldmaps - no estimators will be used.')
            return [], {}

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
            LOGGER.warning(
                f'Several fieldmaps <{", ".join(estimator_key)}> are '
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


def get_MNIInfant_key(spaces: SpatialReferences, res: str | int) -> str:
    """Parse spaces and return matching MNIInfant space, including cohort."""
    for ref in spaces.references:
        if ref.space == 'MNIInfant' and f'res-{res}' in str(ref):
            return ref.fullname

    raise KeyError(f'MNIInfant (resolution {res}) not found in SpatialReferences: {spaces}')
