from pathlib import Path
from unittest.mock import patch

import pytest
from niworkflows.utils.testing import generate_bids_skeleton

from nibabies import config
from nibabies.workflows.anatomical.fit import (
    init_infant_anat_fit_wf,
    init_infant_single_anat_fit_wf,
)
from nibabies.workflows.tests import mock_config

T1_ONLY = {
    '01': {
        'anat': [{'suffix': 'T1w'}],
    }
}

T2_ONLY = {
    '01': {
        'anat': [{'suffix': 'T2w'}],
    }
}

FULL_ANAT = {
    '01': {
        'anat': [{'suffix': 'T1w'}, {'suffix': 'T2w'}],
    },
}


@pytest.fixture(autouse=True)
def _quiet_logger():
    import logging

    logger = logging.getLogger('nipype.workflow')
    old_level = logger.getEffectiveLevel()
    logger.setLevel(logging.ERROR)
    yield
    logger.setLevel(old_level)


def mock_fs_isRunning():
    return False


def _get_default_args(bids_dir, tmp_path):
    from niworkflows.utils.spaces import Reference, SpatialReferences

    return {
        'age_months': 12,
        'bids_root': str(bids_dir),
        'longitudinal': False,
        'omp_nthreads': 1,
        'output_dir': str(tmp_path / 'output'),
        'segmentation_atlases': None,
        'skull_strip_mode': 'skip',
        'skull_strip_template': Reference('MNIInfant', {'cohort': '1'}),
        'skull_strip_fixed_seed': False,
        'sloppy': True,
        'spaces': SpatialReferences([Reference('MNI152NLin6Asym', {})]),
        'cifti_output': False,
    }


@pytest.mark.parametrize('recon_method', ['mcribs', 'freesurfer', 'infantfs'])
@pytest.mark.parametrize('use_aseg', [True, False])
def test_dual_anat_fit_wf(monkeypatch, tmp_path, recon_method, use_aseg):
    monkeypatch.setenv('SUBJECTS_DIR', str(tmp_path / 'subjects'))
    (tmp_path / 'subjects').mkdir()
    bids_dir = tmp_path / 'bids'
    generate_bids_skeleton(bids_dir, FULL_ANAT)

    with mock_config(bids_dir=bids_dir):
        # M-CRIB-S requires an existing output directory
        config.execution.mcribs_dir = tmp_path / 'mcribs'
        config.execution.mcribs_dir.mkdir()

        if use_aseg:
            precomputed = {
                't1w_aseg': str(tmp_path / 'sub-01_space-T1w_desc-aseg_dseg.nii.gz'),
                't2w_aseg': str(tmp_path / 'sub-01_space-T2w_desc-aseg_dseg.nii.gz'),
            }
            Path(precomputed['t1w_aseg']).touch()
            Path(precomputed['t2w_aseg']).touch()
        else:
            precomputed = {}

        kwargs = _get_default_args(bids_dir, tmp_path)
        kwargs.update(
            {
                't1w': [str(bids_dir / 'sub-01' / 'anat' / 'sub-01_T1w.nii.gz')],
                't2w': [str(bids_dir / 'sub-01' / 'anat' / 'sub-01_T2w.nii.gz')],
                'flair': [],
                'precomputed': precomputed,
                'recon_method': recon_method,
                'reference_anat': 'T1w',
            }
        )

        # M-CRIB-S requires aseg to initialize
        if recon_method == 'mcribs' and not use_aseg:
            with patch('smriprep.utils.misc.fs_isRunning', new=mock_fs_isRunning):
                with pytest.raises(NotImplementedError, match='segmentation is required'):
                    init_infant_anat_fit_wf(**kwargs)
            return

        with patch('smriprep.utils.misc.fs_isRunning', new=mock_fs_isRunning):
            init_infant_anat_fit_wf(**kwargs)


@pytest.mark.parametrize('ref_anat', ['T1w', 'T2w'])
@pytest.mark.parametrize('recon_method', ['mcribs', 'freesurfer', 'infantfs'])
@pytest.mark.parametrize('use_aseg', [True, False])
def test_single_anat_fit_wf(monkeypatch, tmp_path, ref_anat, recon_method, use_aseg):
    monkeypatch.setenv('SUBJECTS_DIR', str(tmp_path / 'subjects'))
    (tmp_path / 'subjects').mkdir()
    bids_dir = tmp_path / 'bids'
    layout = T1_ONLY if ref_anat == 'T1w' else T2_ONLY
    generate_bids_skeleton(bids_dir, layout)

    with mock_config(bids_dir=bids_dir):
        # M-CRIB-S requires an existing output directory
        config.execution.mcribs_dir = tmp_path / 'mcribs'
        config.execution.mcribs_dir.mkdir()

        if use_aseg:
            precomputed = {
                f'{ref_anat.lower()}_aseg': str(
                    tmp_path / f'sub-01_space-{ref_anat}_desc-aseg_dseg.nii.gz'
                )
            }
            Path(precomputed[f'{ref_anat.lower()}_aseg']).touch()
        else:
            precomputed = {}

        t1w = (
            [str(bids_dir / 'sub-01' / 'anat' / 'sub-01_T1w.nii.gz')] if ref_anat == 'T1w' else []
        )
        t2w = (
            [str(bids_dir / 'sub-01' / 'anat' / 'sub-01_T2w.nii.gz')] if ref_anat == 'T2w' else []
        )

        kwargs = _get_default_args(bids_dir, tmp_path)
        kwargs.update(
            {
                't1w': t1w,
                't2w': t2w,
                'flair': [],
                'precomputed': precomputed,
                'recon_method': recon_method,
                'reference_anat': ref_anat,
            }
        )

        # M-CRIB-S requires a T2w dataset
        if recon_method == 'mcribs' and ref_anat == 'T1w':
            with patch('smriprep.utils.misc.fs_isRunning', new=mock_fs_isRunning):
                with pytest.raises(ValueError, match='requires use of a T2w image'):
                    init_infant_single_anat_fit_wf(**kwargs)
            return

        # M-CRIB-S requires aseg to initialize
        if recon_method == 'mcribs' and not use_aseg:
            with patch('smriprep.utils.misc.fs_isRunning', new=mock_fs_isRunning):
                with pytest.raises(NotImplementedError, match='segmentation is required'):
                    init_infant_single_anat_fit_wf(**kwargs)
            return

        with patch('smriprep.utils.misc.fs_isRunning', new=mock_fs_isRunning):
            init_infant_single_anat_fit_wf(**kwargs)
