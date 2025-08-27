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
"""Writing out derivative files."""

from __future__ import annotations

import numpy as np
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from niworkflows.utils.images import dseg_label

from nibabies import config
from nibabies._types import Anatomical
from nibabies.config import DEFAULT_DISMISS_ENTITIES, DEFAULT_MEMORY_MIN_GB, dismiss_entities
from nibabies.interfaces import DerivativesDataSink
from nibabies.interfaces.bids import BIDSURI


def prepare_timing_parameters(metadata: dict):
    """Convert initial timing metadata to post-realignment timing metadata

    In particular, SliceTiming metadata is invalid once STC or any realignment is applied,
    as a matrix of voxels no longer corresponds to an acquisition slice.
    Therefore, if SliceTiming is present in the metadata dictionary, and a sparse
    acquisition paradigm is detected, DelayTime or AcquisitionDuration must be derived to
    preserve the timing interpretation.

    Examples
    --------

    .. testsetup::

        >>> from unittest import mock

    If SliceTiming metadata is absent, then the only change is to note that
    STC has not been applied:

    >>> prepare_timing_parameters(dict(RepetitionTime=2))
    {'RepetitionTime': 2, 'SliceTimingCorrected': False}
    >>> prepare_timing_parameters(dict(RepetitionTime=2, DelayTime=0.5))
    {'RepetitionTime': 2, 'DelayTime': 0.5, 'SliceTimingCorrected': False}
    >>> prepare_timing_parameters(dict(VolumeTiming=[0.0, 1.0, 2.0, 5.0, 6.0, 7.0],
    ...                                AcquisitionDuration=1.0))  #doctest: +NORMALIZE_WHITESPACE
    {'VolumeTiming': [0.0, 1.0, 2.0, 5.0, 6.0, 7.0], 'AcquisitionDuration': 1.0,
     'SliceTimingCorrected': False}

    When SliceTiming is available and used, then ``SliceTimingCorrected`` is ``True``
    and the ``StartTime`` indicates a series offset.

    >>> with mock.patch("nibabies.config.workflow.ignore", []):
    ...     prepare_timing_parameters(dict(RepetitionTime=2, SliceTiming=[0.0, 0.2, 0.4, 0.6]))
    {'RepetitionTime': 2, 'SliceTimingCorrected': True, 'DelayTime': 1.2, 'StartTime': 0.3}
    >>> with mock.patch("nibabies.config.workflow.ignore", []):
    ...     prepare_timing_parameters(
    ...         dict(VolumeTiming=[0.0, 1.0, 2.0, 5.0, 6.0, 7.0],
    ...              SliceTiming=[0.0, 0.2, 0.4, 0.6, 0.8]))  #doctest: +NORMALIZE_WHITESPACE
    {'VolumeTiming': [0.0, 1.0, 2.0, 5.0, 6.0, 7.0], 'SliceTimingCorrected': True,
     'AcquisitionDuration': 1.0, 'StartTime': 0.4}

    When SliceTiming is available and not used, then ``SliceTimingCorrected`` is ``False``
    and TA is indicated with ``DelayTime`` or ``AcquisitionDuration``.

    >>> with mock.patch("nibabies.config.workflow.ignore", ["slicetiming"]):
    ...     prepare_timing_parameters(dict(RepetitionTime=2, SliceTiming=[0.0, 0.2, 0.4, 0.6]))
    {'RepetitionTime': 2, 'SliceTimingCorrected': False, 'DelayTime': 1.2}
    >>> with mock.patch("nibabies.config.workflow.ignore", ["slicetiming"]):
    ...     prepare_timing_parameters(
    ...         dict(VolumeTiming=[0.0, 1.0, 2.0, 5.0, 6.0, 7.0],
    ...              SliceTiming=[0.0, 0.2, 0.4, 0.6, 0.8]))  #doctest: +NORMALIZE_WHITESPACE
    {'VolumeTiming': [0.0, 1.0, 2.0, 5.0, 6.0, 7.0], 'SliceTimingCorrected': False,
     'AcquisitionDuration': 1.0}

    If SliceTiming metadata is present but empty, then treat it as missing:

    >>> with mock.patch("nibabies.config.workflow.ignore", []):
    ...     prepare_timing_parameters(dict(RepetitionTime=2, SliceTiming=[]))
    {'RepetitionTime': 2, 'SliceTimingCorrected': False}
    >>> with mock.patch("nibabies.config.workflow.ignore", []):
    ...     prepare_timing_parameters(dict(RepetitionTime=2, SliceTiming=[0.0]))
    {'RepetitionTime': 2, 'SliceTimingCorrected': False}
    """
    timing_parameters = {
        key: metadata[key]
        for key in (
            'RepetitionTime',
            'VolumeTiming',
            'DelayTime',
            'AcquisitionDuration',
            'SliceTiming',
        )
        if key in metadata
    }

    # Treat SliceTiming of [] or length 1 as equivalent to missing and remove it in any case
    slice_timing = timing_parameters.pop('SliceTiming', [])

    run_stc = len(slice_timing) > 1 and 'slicetiming' not in config.workflow.ignore
    timing_parameters['SliceTimingCorrected'] = run_stc

    if len(slice_timing) > 1:
        st = sorted(slice_timing)
        TA = st[-1] + (st[1] - st[0])  # Final slice onset + slice duration
        # For constant TR paradigms, use DelayTime
        if 'RepetitionTime' in timing_parameters:
            TR = timing_parameters['RepetitionTime']
            if not np.isclose(TR, TA) and TA < TR:
                timing_parameters['DelayTime'] = TR - TA
        # For variable TR paradigms, use AcquisitionDuration
        elif 'VolumeTiming' in timing_parameters:
            timing_parameters['AcquisitionDuration'] = TA

        if run_stc:
            first, last = st[0], st[-1]
            frac = config.workflow.slice_time_ref
            tzero = float(np.round(first + frac * (last - first), 3))
            timing_parameters['StartTime'] = tzero

    return timing_parameters


def init_func_fit_reports_wf(
    *,
    reference_anat: Anatomical,
    sdc_correction: bool,
    output_dir: str,
    name='func_fit_reports_wf',
) -> pe.Workflow:
    """
    Set up a battery of datasinks to store reports in the right location.

    Parameters
    ----------
    reference_anat : :obj:`str`
        Reference anatomical image type (T1w, T2w)
    sdc_correction : :obj:`bool`
        SDC was performed
    output_dir : :obj:`str`
        Directory in which to save derivatives
    name : :obj:`str`
        Workflow name (default: func_fit_reports_wf)

    Inputs
    ------
    source_file
        Input BOLD images
    subject_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    anat_preproc
        The anatomical reference map, which is calculated as the average of bias-corrected
        and preprocessed anatomical images, defining the anatomical space.
    anat_dseg
        Segmentation in anatomical reference space
    anat_mask
        Brain (binary) mask estimated by brain extraction.
    """
    from nireports.interfaces.reporting.base import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from sdcflows.interfaces.reportlets import FieldmapReportlet

    workflow = pe.Workflow(name=name)

    inputfields = [
        'source_file',
        'sdc_boldref',
        'coreg_boldref',
        'bold_mask',
        'boldref2anat_xfm',
        'boldref2fmap_xfm',
        'anat_preproc',
        'anat_mask',
        'anat_dseg',
        'fieldmap',
        'fmap_ref',
        # May be missing
        'subject_id',
        'subjects_dir',
        # Report snippets
        'summary_report',
        'validation_report',
    ]
    inputnode = pe.Node(niu.IdentityInterface(fields=inputfields), name='inputnode')

    ds_summary = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='summary',
            datatype='figures',
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
        ),
        name='ds_report_summary',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    ds_validation = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='validation',
            datatype='figures',
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
        ),
        name='ds_report_validation',
        run_without_submitting=True,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )

    # Resample anatomical references into BOLD space for plotting
    anat_boldref = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            float=True,
            invert_transform_flags=[True],
            interpolation='LanczosWindowedSinc',
        ),
        name='anat_boldref',
        mem_gb=1,
    )

    anat_wm = pe.Node(
        niu.Function(function=dseg_label),
        name='anat_wm',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    # TODO: Verify this is consistent with T1/T2
    anat_wm.inputs.label = 2  # BIDS default is WM=2

    boldref_wm = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            invert_transform_flags=[True],
            interpolation='NearestNeighbor',
        ),
        name='boldref_wm',
        mem_gb=1,
    )

    workflow.connect([
        (inputnode, ds_summary, [
            ('source_file', 'source_file'),
            ('summary_report', 'in_file'),
        ]),
        (inputnode, ds_validation, [
            ('source_file', 'source_file'),
            ('validation_report', 'in_file'),
        ]),
        (inputnode, anat_boldref, [
            ('anat_preproc', 'input_image'),
            ('coreg_boldref', 'reference_image'),
            ('boldref2anat_xfm', 'transforms'),
        ]),
        (inputnode, anat_wm, [('anat_dseg', 'in_seg')]),
        (inputnode, boldref_wm, [
            ('coreg_boldref', 'reference_image'),
            ('boldref2anat_xfm', 'transforms'),
        ]),
        (anat_wm, boldref_wm, [('out', 'input_image')]),
    ])  # fmt:skip

    # Reportlets follow the structure of init_bold_fit_wf stages
    # - SDC1:
    #       Before: Pre-SDC boldref
    #       After: Fieldmap reference resampled on boldref
    #       Three-way: Fieldmap resampled on boldref
    # - SDC2:
    #       Before: Pre-SDC boldref with white matter mask
    #       After: Post-SDC boldref with white matter mask
    # - EPI-anat registration:
    #       Before: anat brain with white matter mask
    #       After: Resampled boldref with white matter mask

    if sdc_correction:
        fmapref_boldref = pe.Node(
            ApplyTransforms(
                dimension=3,
                default_value=0,
                float=True,
                invert_transform_flags=[True],
                interpolation='LanczosWindowedSinc',
            ),
            name='fmapref_boldref',
            mem_gb=1,
        )

        # SDC1
        sdcreg_report = pe.Node(
            FieldmapReportlet(
                reference_label='BOLD reference',
                moving_label='Fieldmap reference',
                show='both',
            ),
            name='sdecreg_report',
            mem_gb=0.1,
        )

        ds_sdcreg_report = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc='fmapCoreg',
                suffix='bold',
                datatype='figures',
                dismiss_entities=DEFAULT_DISMISS_ENTITIES,
            ),
            name='ds_sdcreg_report',
        )

        # SDC2
        sdc_report = pe.Node(
            SimpleBeforeAfter(
                before_label='Distorted',
                after_label='Corrected',
                dismiss_affine=True,
            ),
            name='sdc_report',
            mem_gb=0.1,
        )

        ds_sdc_report = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc='sdc',
                suffix='bold',
                datatype='figures',
                dismiss_entities=DEFAULT_DISMISS_ENTITIES,
            ),
            name='ds_sdc_report',
        )

        workflow.connect([
            (inputnode, fmapref_boldref, [
                ('fmap_ref', 'input_image'),
                ('coreg_boldref', 'reference_image'),
                ('boldref2fmap_xfm', 'transforms'),
            ]),
            (inputnode, sdcreg_report, [
                ('sdc_boldref', 'reference'),
                ('fieldmap', 'fieldmap'),
                ('bold_mask', 'mask'),
            ]),
            (fmapref_boldref, sdcreg_report, [('output_image', 'moving')]),
            (inputnode, ds_sdcreg_report, [('source_file', 'source_file')]),
            (sdcreg_report, ds_sdcreg_report, [('out_report', 'in_file')]),
            (inputnode, sdc_report, [
                ('sdc_boldref', 'before'),
                ('coreg_boldref', 'after'),
            ]),
            (boldref_wm, sdc_report, [('output_image', 'wm_seg')]),
            (inputnode, ds_sdc_report, [('source_file', 'source_file')]),
            (sdc_report, ds_sdc_report, [('out_report', 'in_file')]),
        ])  # fmt:skip

    # EPI-anat registration
    # Resample anat image onto EPI-space

    epi_anat_report = pe.Node(
        SimpleBeforeAfter(
            before_label=reference_anat,
            after_label='EPI',
            dismiss_affine=True,
        ),
        name='epi_anat_report',
        mem_gb=0.1,
    )

    ds_epi_anat_report = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='coreg',
            suffix='bold',
            datatype='figures',
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
        ),
        name='ds_epi_anat_report',
    )

    workflow.connect([
        (inputnode, epi_anat_report, [('coreg_boldref', 'after')]),
        (anat_boldref, epi_anat_report, [('output_image', 'before')]),
        (boldref_wm, epi_anat_report, [('output_image', 'wm_seg')]),
        (inputnode, ds_epi_anat_report, [('source_file', 'source_file')]),
        (epi_anat_report, ds_epi_anat_report, [('out_report', 'in_file')]),
    ])  # fmt:skip

    return workflow


def init_ds_boldref_wf(
    *,
    output_dir,
    desc: str,
    name='ds_boldref_wf',
) -> pe.Workflow:
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['source_files', 'boldref']),
        name='inputnode',
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=['boldref']), name='outputnode')

    sources = pe.Node(
        BIDSURI(
            numinputs=1,
            dataset_links=config.execution.dataset_links,
            out_dir=str(config.execution.output_dir.absolute()),
        ),
        name='sources',
    )

    ds_boldref = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc=desc,
            suffix='boldref',
            compress=True,
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
        ),
        name='ds_boldref',
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, sources, [('source_files', 'in1')]),
        (inputnode, ds_boldref, [('boldref', 'in_file'),
                                 ('source_files', 'source_file')]),
        (sources, ds_boldref, [('out', 'Sources')]),
        (ds_boldref, outputnode, [('out_file', 'boldref')]),
    ])
    # fmt:on

    return workflow


def init_ds_registration_wf(
    *,
    output_dir: str,
    source: str,
    dest: str,
    name: str,
) -> pe.Workflow:
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['source_files', 'xform']),
        name='inputnode',
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=['xform']), name='outputnode')

    sources = pe.Node(
        BIDSURI(
            numinputs=1,
            dataset_links=config.execution.dataset_links,
            out_dir=str(config.execution.output_dir.absolute()),
        ),
        name='sources',
    )

    ds_xform = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            mode='image',
            suffix='xfm',
            extension='.txt',
            dismiss_entities=dismiss_entities(['part']),
            **{'from': source, 'to': dest},
        ),
        name='ds_xform',
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    # fmt:off
    workflow.connect([
        (inputnode, sources, [('source_files', 'in1')]),
        (inputnode, ds_xform, [('xform', 'in_file'),
                               ('source_files', 'source_file')]),
        (sources, ds_xform, [('out', 'Sources')]),
        (ds_xform, outputnode, [('out_file', 'xform')]),
    ])
    # fmt:on

    return workflow


def init_ds_hmc_wf(
    *,
    output_dir,
    name='ds_hmc_wf',
) -> pe.Workflow:
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['source_files', 'xforms']),
        name='inputnode',
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=['xforms']), name='outputnode')

    sources = pe.Node(
        BIDSURI(
            numinputs=1,
            dataset_links=config.execution.dataset_links,
            out_dir=str(config.execution.output_dir.absolute()),
        ),
        name='sources',
    )

    ds_xforms = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='hmc',
            suffix='xfm',
            extension='.txt',
            compress=True,
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
            **{'from': 'orig', 'to': 'boldref'},
        ),
        name='ds_xforms',
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, sources, [('source_files', 'in1')]),
        (inputnode, ds_xforms, [('xforms', 'in_file'),
                                ('source_files', 'source_file')]),
        (sources, ds_xforms, [('out', 'Sources')]),
        (ds_xforms, outputnode, [('out_file', 'xforms')]),
    ])
    # fmt:on

    return workflow


def init_ds_bold_native_wf(
    *,
    bids_root: str,
    output_dir: str,
    multiecho: bool,
    bold_output: bool,
    echo_output: bool,
    all_metadata: list[dict],
    name='ds_bold_native_wf',
) -> pe.Workflow:
    metadata = all_metadata[0]
    timing_parameters = prepare_timing_parameters(metadata)

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'source_files',
                'bold',
                'bold_mask',
                'bold_echos',
                't2star',
                # Transforms previously used to generate the outputs
                'motion_xfm',
                'boldref2fmap_xfm',
            ]
        ),
        name='inputnode',
    )

    sources = pe.Node(
        BIDSURI(
            numinputs=3,
            dataset_links=config.execution.dataset_links,
            out_dir=str(config.execution.output_dir.absolute()),
        ),
        name='sources',
    )
    workflow.connect([
        (inputnode, sources, [
            ('source_files', 'in1'),
            ('motion_xfm', 'in2'),
            ('boldref2fmap_xfm', 'in3'),
        ]),
    ])  # fmt:skip

    # Masks should be output if any other derivatives are output
    ds_bold_mask = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='brain',
            suffix='mask',
            compress=True,
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
        ),
        name='ds_bold_mask',
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (inputnode, ds_bold_mask, [
            ('source_files', 'source_file'),
            ('bold_mask', 'in_file'),
        ]),
        (sources, ds_bold_mask, [('out', 'Sources')]),
    ])  # fmt:skip

    if bold_output:
        ds_bold = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc='preproc',
                compress=True,
                SkullStripped=multiecho,
                TaskName=metadata.get('TaskName'),
                dismiss_entities=DEFAULT_DISMISS_ENTITIES,
                **timing_parameters,
            ),
            name='ds_bold',
            mem_gb=DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, ds_bold, [
                ('source_files', 'source_file'),
                ('bold', 'in_file'),
            ]),
            (sources, ds_bold, [('out', 'Sources')]),
        ])  # fmt:skip

    if bold_output and multiecho:
        t2star_meta = {
            'Units': 's',
            'EstimationReference': 'doi:10.1002/mrm.20900',
            'EstimationAlgorithm': 'monoexponential decay model',
        }
        ds_t2star = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                space='boldref',
                suffix='T2starmap',
                compress=True,
                dismiss_entities=DEFAULT_DISMISS_ENTITIES,
                **t2star_meta,
            ),
            name='ds_t2star_bold',
            run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB,
        )
        workflow.connect([
            (inputnode, ds_t2star, [
                ('source_files', 'source_file'),
                ('t2star', 'in_file'),
            ]),
            (sources, ds_t2star, [('out', 'Sources')]),
        ])  # fmt:skip

    if echo_output:
        ds_bold_echos = pe.MapNode(
            DerivativesDataSink(
                base_directory=output_dir,
                desc='preproc',
                compress=True,
                SkullStripped=False,
                TaskName=metadata.get('TaskName'),
                **timing_parameters,
            ),
            iterfield=['source_file', 'in_file', 'meta_dict'],
            name='ds_bold_echos',
            run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB,
        )
        ds_bold_echos.inputs.meta_dict = [{'EchoTime': md['EchoTime']} for md in all_metadata]
        workflow.connect([
            (inputnode, ds_bold_echos, [
                ('source_files', 'source_file'),
                ('bold_echos', 'in_file'),
            ]),
        ])  # fmt:skip

    return workflow


def init_ds_volumes_wf(
    *,
    bids_root: str,
    output_dir: str,
    multiecho: bool,
    metadata: list[dict],
    name='ds_volumes_wf',
) -> pe.Workflow:
    timing_parameters = prepare_timing_parameters(metadata)

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'source_files',
                'ref_file',
                'bold',  # Resampled into target space
                'bold_mask',  # boldref space
                'bold_ref',  # boldref space
                't2star',  # boldref space
                'template',  # target reference image from original transform
                # Anatomical
                'boldref2anat_xfm',
                # Template
                'anat2std_xfm',
                # Entities
                'space',
                'cohort',
                'resolution',
                # Transforms previously used to generate the outputs
                'motion_xfm',
                'boldref2fmap_xfm',
            ]
        ),
        name='inputnode',
    )

    sources = pe.Node(
        BIDSURI(
            numinputs=6,
            dataset_links=config.execution.dataset_links,
            out_dir=str(config.execution.output_dir.absolute()),
        ),
        name='sources',
    )
    boldref2target = pe.Node(niu.Merge(2), name='boldref2target')

    # BOLD is pre-resampled
    ds_bold = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='preproc',
            compress=True,
            SkullStripped=multiecho,
            TaskName=metadata.get('TaskName'),
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
            **timing_parameters,
        ),
        name='ds_bold',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    workflow.connect([
        (inputnode, sources, [
            ('source_files', 'in1'),
            ('motion_xfm', 'in2'),
            ('boldref2fmap_xfm', 'in3'),
            ('boldref2anat_xfm', 'in4'),
            ('anat2std_xfm', 'in5'),
            ('template', 'in6'),
        ]),
        (inputnode, boldref2target, [
            # Note that ANTs expects transforms in target-to-source order
            # Reverse this for nitransforms-based resamplers
            ('anat2std_xfm', 'in1'),
            ('boldref2anat_xfm', 'in2'),
        ]),
        (inputnode, ds_bold, [
            ('source_files', 'source_file'),
            ('bold', 'in_file'),
            ('space', 'space'),
            ('cohort', 'cohort'),
            ('resolution', 'resolution'),
        ]),
        (sources, ds_bold, [('out', 'Sources')]),
    ])  # fmt:skip

    resample_ref = pe.Node(
        ApplyTransforms(
            dimension=3,
            default_value=0,
            float=True,
            interpolation='LanczosWindowedSinc',
        ),
        name='resample_ref',
    )
    resample_mask = pe.Node(ApplyTransforms(interpolation='MultiLabel'), name='resample_mask')
    resamplers = [resample_ref, resample_mask]

    workflow.connect([
        (inputnode, resample_ref, [('bold_ref', 'input_image')]),
        (inputnode, resample_mask, [('bold_mask', 'input_image')]),
    ])  # fmt:skip

    ds_ref = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            suffix='boldref',
            compress=True,
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
        ),
        name='ds_ref',
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    ds_mask = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='brain',
            suffix='mask',
            compress=True,
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
        ),
        name='ds_mask',
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    datasinks = [ds_ref, ds_mask]

    if multiecho:
        t2star_meta = {
            'Units': 's',
            'EstimationReference': 'doi:10.1002/mrm.20900',
            'EstimationAlgorithm': 'monoexponential decay model',
        }
        resample_t2star = pe.Node(
            ApplyTransforms(
                dimension=3,
                default_value=0,
                float=True,
                interpolation='LanczosWindowedSinc',
            ),
            name='resample_t2star',
        )
        ds_t2star = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                suffix='T2starmap',
                compress=True,
                dismiss_entities=DEFAULT_DISMISS_ENTITIES,
                **t2star_meta,
            ),
            name='ds_t2star_std',
            run_without_submitting=True,
            mem_gb=DEFAULT_MEMORY_MIN_GB,
        )
        resamplers.append(resample_t2star)
        datasinks.append(ds_t2star)

        workflow.connect([(inputnode, resample_t2star, [('t2star', 'input_image')])])

    workflow.connect(
        [
            (inputnode, resampler, [('ref_file', 'reference_image')])
            for resampler in resamplers
        ] + [
            (boldref2target, resampler, [('out', 'transforms')])
            for resampler in resamplers
        ] + [
            (inputnode, datasink, [
                ('source_files', 'source_file'),
                ('space', 'space'),
                ('cohort', 'cohort'),
                ('resolution', 'resolution'),
            ])
            for datasink in datasinks
        ] + [
            (sources, datasink, [('out', 'Sources')])
            for datasink in datasinks
        ] + [
            (resampler, datasink, [('output_image', 'in_file')])
            for resampler, datasink in zip(resamplers, datasinks, strict=False)
        ]
    )  # fmt:skip

    return workflow


def init_bold_preproc_report_wf(
    mem_gb: float,
    reportlets_dir: str,
    name: str = 'bold_preproc_report_wf',
):
    """
    Generate a visual report.

    This workflow generates and saves a reportlet showing the effect of resampling
    the BOLD signal using the standard deviation maps.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from nibabies.workflows.bold.resampling import init_bold_preproc_report_wf
            wf = init_bold_preproc_report_wf(mem_gb=1, reportlets_dir='.')

    Parameters
    ----------
    mem_gb : :obj:`float`
        Size of BOLD file in GB
    reportlets_dir : :obj:`str`
        Directory in which to save reportlets
    name : :obj:`str`, optional
        Workflow name (default: bold_preproc_report_wf)

    Inputs
    ------
    in_pre
        BOLD time-series, before resampling
    in_post
        BOLD time-series, after resampling
    name_source
        BOLD series NIfTI file
        Used to recover original information lost during processing

    """
    from nipype.algorithms.confounds import TSNR
    from nireports.interfaces.reporting.base import SimpleBeforeAfterRPT
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    from ...interfaces import DerivativesDataSink

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['in_pre', 'in_post', 'name_source']), name='inputnode'
    )

    pre_tsnr = pe.Node(TSNR(), name='pre_tsnr', mem_gb=mem_gb * 4.5)
    pos_tsnr = pe.Node(TSNR(), name='pos_tsnr', mem_gb=mem_gb * 4.5)

    bold_rpt = pe.Node(SimpleBeforeAfterRPT(), name='bold_rpt', mem_gb=0.1)
    ds_report_bold = pe.Node(
        DerivativesDataSink(
            base_directory=reportlets_dir,
            desc='preproc',
            datatype='figures',
            dismiss_entities=DEFAULT_DISMISS_ENTITIES,
        ),
        name='ds_report_bold',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )
    # fmt:off
    workflow.connect([
        (inputnode, ds_report_bold, [('name_source', 'source_file')]),
        (inputnode, pre_tsnr, [('in_pre', 'in_file')]),
        (inputnode, pos_tsnr, [('in_post', 'in_file')]),
        (pre_tsnr, bold_rpt, [('stddev_file', 'before')]),
        (pos_tsnr, bold_rpt, [('stddev_file', 'after')]),
        (bold_rpt, ds_report_bold, [('out_report', 'in_file')]),
    ])
    # fmt:on

    return workflow
