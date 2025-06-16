from copy import deepcopy
from pathlib import Path
from unittest.mock import patch

import bids
import nibabel as nb
import numpy as np
import pandas as pd
import pytest
from nipype.pipeline.engine.utils import generate_expanded_graph
from niworkflows.utils.testing import generate_bids_skeleton
from sdcflows.fieldmaps import clear_registry
from sdcflows.utils.wrangler import find_estimators

from ... import config
from ..base import get_estimator, init_nibabies_wf
from ..tests import mock_config

BASE_LAYOUT = {
    '01': {
        'anat': [
            {'run': 1, 'suffix': 'T1w'},
            {'run': 2, 'suffix': 'T1w'},
            {'suffix': 'T2w'},
        ],
        'func': [
            *(
                {
                    'task': 'rest',
                    'run': i,
                    'suffix': suffix,
                    'metadata': {
                        'RepetitionTime': 2.0,
                        'PhaseEncodingDirection': 'j',
                        'TotalReadoutTime': 0.6,
                        'EchoTime': 0.03,
                        'SliceTiming': [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8],
                    },
                }
                for suffix in ('bold', 'sbref')
                for i in range(1, 3)
            ),
            *(
                {
                    'task': 'nback',
                    'echo': i,
                    'suffix': 'bold',
                    'metadata': {
                        'RepetitionTime': 2.0,
                        'PhaseEncodingDirection': 'j',
                        'TotalReadoutTime': 0.6,
                        'EchoTime': 0.015 * i,
                        'SliceTiming': [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8],
                    },
                }
                for i in range(1, 4)
            ),
        ],
        'fmap': [
            {'suffix': 'phasediff', 'metadata': {'EchoTime1': 0.005, 'EchoTime2': 0.007}},
            {'suffix': 'magnitude1', 'metadata': {'EchoTime': 0.005}},
            {
                'suffix': 'epi',
                'direction': 'PA',
                'metadata': {'PhaseEncodingDirection': 'j', 'TotalReadoutTime': 0.6},
            },
            {
                'suffix': 'epi',
                'direction': 'AP',
                'metadata': {'PhaseEncodingDirection': 'j-', 'TotalReadoutTime': 0.6},
            },
        ],
    },
}


T1W_ONLY = BASE_LAYOUT.copy()
T1W_ONLY['01']['anat'] = [{'suffix': 'T1w'}]

T2W_ONLY = T1W_ONLY.copy()
T2W_ONLY['01']['anat'] = [{'suffix': 'T2w'}]


@pytest.fixture(scope='module', autouse=True)
def _quiet_logger():
    import logging

    logger = logging.getLogger('nipype.workflow')
    old_level = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield
    logger.setLevel(old_level)


@pytest.fixture(autouse=True)
def _reset_sdcflows_registry():
    yield
    clear_registry()


@pytest.fixture(scope='module', params=[BASE_LAYOUT, T1W_ONLY, T2W_ONLY])
def bids_root(tmp_path_factory, request):
    """
    Create a BIDS skeleton for various input types:
    - Full inputs (T1w + T2w)
    - T1w only
    - T2w only
    """
    base = tmp_path_factory.mktemp('base')
    bids_dir = base / 'bids'
    generate_bids_skeleton(bids_dir, request.param)

    # Ensure age information is available
    df = pd.DataFrame({'participant_id': ['sub-01'], 'age_months': ['1']})
    df.to_csv(str(bids_dir / 'participants.tsv'), sep='\t', index=False)

    img = nb.Nifti1Image(np.zeros((10, 10, 10, 10)), np.eye(4))

    for bold_path in bids_dir.glob('sub-01/*/*.nii.gz'):
        img.to_filename(bold_path)

    return bids_dir


def _make_params(
    bold2anat_init: str = 'auto',
    use_bbr: bool | None = None,
    dummy_scans: int | None = None,
    me_output_echos: bool = False,
    medial_surface_nan: bool = False,
    project_goodvoxels: bool = False,
    cifti_output: bool | str = False,
    run_msmsulc: bool = True,
    skull_strip_anat: str = 'auto',
    use_syn_sdc: str | bool = False,
    force_syn: bool = False,
    surface_recon_method: str | None = 'auto',
    ignore: list[str] = None,
    bids_filters: dict = None,
    norm_csf: bool = False,
    multi_step_reg: bool = True,
):
    if ignore is None:
        ignore = []
    if bids_filters is None:
        bids_filters = {}
    return (
        bold2anat_init,
        use_bbr,
        dummy_scans,
        me_output_echos,
        medial_surface_nan,
        project_goodvoxels,
        cifti_output,
        run_msmsulc,
        skull_strip_anat,
        use_syn_sdc,
        force_syn,
        surface_recon_method,
        ignore,
        bids_filters,
        norm_csf,
        multi_step_reg,
    )


@pytest.mark.parametrize('level', ['minimal', 'resampling', 'full'])
@pytest.mark.parametrize('anat_only', [False, True])
@pytest.mark.parametrize(
    (
        'bold2anat_init',
        'use_bbr',
        'dummy_scans',
        'me_output_echos',
        'medial_surface_nan',
        'project_goodvoxels',
        'cifti_output',
        'run_msmsulc',
        'skull_strip_anat',
        'use_syn_sdc',
        'force_syn',
        'surface_recon_method',
        'ignore',
        'bids_filters',
        'norm_csf',
        'multi_step_reg',
    ),
    [
        _make_params(),
        _make_params(bold2anat_init='t1w'),
        _make_params(bold2anat_init='t2w'),
        _make_params(bold2anat_init='header'),
        _make_params(use_bbr=True),
        _make_params(use_bbr=False),
        _make_params(bold2anat_init='header', use_bbr=True),
        # Currently disabled
        # _make_params(bold2anat_init="header", use_bbr=False),
        _make_params(dummy_scans=2),
        _make_params(me_output_echos=True),
        _make_params(medial_surface_nan=True),
        _make_params(cifti_output='91k'),
        _make_params(cifti_output='91k', project_goodvoxels=True),
        _make_params(cifti_output='91k', run_msmsulc=False),
        _make_params(skull_strip_anat='force'),
        _make_params(skull_strip_anat='skip'),
        # _make_params(use_syn_sdc='warn', force_syn=True, ignore=['fieldmaps']),
        _make_params(surface_recon_method=None),
        _make_params(surface_recon_method=None, use_bbr=True),
        _make_params(surface_recon_method=None, use_bbr=False),
        _make_params(surface_recon_method='auto'),
        # _make_params(surface_recon_method='mcribs'),  # Requires precomputed segmentation
        _make_params(surface_recon_method='infantfs'),
        _make_params(surface_recon_method='freesurfer'),
        # Currently unsupported:
        # _make_params(freesurfer=False, bold2anat_init="header"),
        # _make_params(freesurfer=False, bold2anat_init="header", use_bbr=True),
        # _make_params(freesurfer=False, bold2anat_init="header", use_bbr=False),
        # Regression test for gh-3154:
        _make_params(bids_filters={'sbref': {'suffix': 'sbref'}}),
        _make_params(norm_csf=True),
        _make_params(multi_step_reg=True),
        _make_params(multi_step_reg=False),
    ],
)
def test_init_nibabies_wf(
    monkeypatch,
    bids_root: Path,
    tmp_path: Path,
    level: str,
    anat_only: bool,
    bold2anat_init: str,
    use_bbr: bool | None,
    dummy_scans: int | None,
    me_output_echos: bool,
    medial_surface_nan: bool,
    project_goodvoxels: bool,
    cifti_output: bool | str,
    run_msmsulc: bool,
    skull_strip_anat: str,
    use_syn_sdc: str | bool,
    force_syn: bool,
    surface_recon_method: str | None,
    ignore: list[str],
    bids_filters: dict,
    norm_csf: bool,
    multi_step_reg: bool,
):
    monkeypatch.setenv('SUBJECTS_DIR', '/opt/freesurfer/subjects')
    with mock_config(bids_dir=bids_root):
        config.workflow.level = level
        config.workflow.anat_only = anat_only
        config.workflow.bold2anat_init = bold2anat_init
        config.workflow.use_bbr = use_bbr
        config.workflow.dummy_scans = dummy_scans
        config.execution.me_output_echos = me_output_echos
        config.workflow.medial_surface_nan = medial_surface_nan
        config.workflow.project_goodvoxels = project_goodvoxels
        config.workflow.norm_csf = norm_csf
        config.workflow.multi_step_reg = multi_step_reg
        # config.workflow.run_msmsulc = run_msmsulc
        config.workflow.skull_strip_anat = skull_strip_anat
        config.workflow.cifti_output = cifti_output
        config.workflow.surface_recon_method = surface_recon_method
        config.workflow.ignore = ignore
        with patch.dict('nibabies.config.execution.bids_filters', bids_filters):
            wf = init_nibabies_wf(
                [
                    ('01', None),
                ]
            )

    generate_expanded_graph(wf._create_flat_graph())


def test_get_estimator_none(tmp_path):
    bids_dir = tmp_path / 'bids'

    # No IntendedFors/B0Fields
    generate_bids_skeleton(bids_dir, BASE_LAYOUT)
    layout = bids.BIDSLayout(bids_dir)
    bold_files = sorted(
        layout.get(suffix='bold', task='rest', extension='.nii.gz', return_type='file')
    )

    assert get_estimator(layout, bold_files[0]) == ()
    assert get_estimator(layout, bold_files[1]) == ()


def test_get_estimator_b0field_and_intendedfor(tmp_path):
    bids_dir = tmp_path / 'bids'

    # Set B0FieldSource for run 1
    spec = deepcopy(BASE_LAYOUT)
    spec['01']['func'][0]['metadata']['B0FieldSource'] = 'epi'
    spec['01']['fmap'][2]['metadata']['B0FieldIdentifier'] = 'epi'
    spec['01']['fmap'][3]['metadata']['B0FieldIdentifier'] = 'epi'

    # Set IntendedFor for run 2
    spec['01']['fmap'][0]['metadata']['IntendedFor'] = 'func/sub-01_task-rest_run-2_bold.nii.gz'

    generate_bids_skeleton(bids_dir, spec)
    layout = bids.BIDSLayout(bids_dir)
    _ = find_estimators(layout=layout, subject='01')

    bold_files = sorted(
        layout.get(suffix='bold', task='rest', extension='.nii.gz', return_type='file')
    )

    assert get_estimator(layout, bold_files[0]) == ('epi',)
    assert get_estimator(layout, bold_files[1]) == ('auto_00000',)


def test_get_estimator_overlapping_specs(tmp_path):
    bids_dir = tmp_path / 'bids'

    # Set B0FieldSource for both runs
    spec = deepcopy(BASE_LAYOUT)
    spec['01']['func'][0]['metadata']['B0FieldSource'] = 'epi'
    spec['01']['func'][1]['metadata']['B0FieldSource'] = 'epi'
    spec['01']['fmap'][2]['metadata']['B0FieldIdentifier'] = 'epi'
    spec['01']['fmap'][3]['metadata']['B0FieldIdentifier'] = 'epi'

    # Set IntendedFor for both runs
    spec['01']['fmap'][0]['metadata']['IntendedFor'] = [
        'func/sub-01_task-rest_run-1_bold.nii.gz',
        'func/sub-01_task-rest_run-2_bold.nii.gz',
    ]

    generate_bids_skeleton(bids_dir, spec)
    layout = bids.BIDSLayout(bids_dir)
    _ = find_estimators(layout=layout, subject='01')

    bold_files = sorted(
        layout.get(suffix='bold', task='rest', extension='.nii.gz', return_type='file')
    )

    # B0Fields take precedence
    assert get_estimator(layout, bold_files[0]) == ('epi',)
    assert get_estimator(layout, bold_files[1]) == ('epi',)


def test_get_estimator_multiple_b0fields(tmp_path):
    bids_dir = tmp_path / 'bids'

    # Set B0FieldSource for both runs
    spec = deepcopy(BASE_LAYOUT)
    spec['01']['func'][0]['metadata']['B0FieldSource'] = ('epi', 'phasediff')
    spec['01']['func'][1]['metadata']['B0FieldSource'] = 'epi'
    spec['01']['fmap'][0]['metadata']['B0FieldIdentifier'] = 'phasediff'
    spec['01']['fmap'][1]['metadata']['B0FieldIdentifier'] = 'phasediff'
    spec['01']['fmap'][2]['metadata']['B0FieldIdentifier'] = 'epi'
    spec['01']['fmap'][3]['metadata']['B0FieldIdentifier'] = 'epi'

    generate_bids_skeleton(bids_dir, spec)
    layout = bids.BIDSLayout(bids_dir)
    _ = find_estimators(layout=layout, subject='01')

    bold_files = sorted(
        layout.get(suffix='bold', task='rest', extension='.nii.gz', return_type='file')
    )

    # Always get an iterable; don't care if it's a list or tuple
    assert get_estimator(layout, bold_files[0]) == ['epi', 'phasediff']
    assert get_estimator(layout, bold_files[1]) == ('epi',)
