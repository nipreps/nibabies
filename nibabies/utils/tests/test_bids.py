from pathlib import Path

import nibabel as nb
import numpy as np
import pytest

from nibabies.utils.bids import (
    _get_age_from_tsv,
    bold_template_zooms,
    is_valid_bold_template,
)


class _MockLayout:
    """Minimal BIDSLayout stub for is_valid_bold_template tests."""

    def __init__(self, pe_map: dict):
        self._pe_map = pe_map

    def get_metadata(self, f):
        pe = self._pe_map.get(f)
        return {'PhaseEncodingDirection': pe} if pe is not None else {}


TEST_ZOOMS = (2.4, 2.4, 2.4)


# Each run is a (fieldmap, phase-encoding direction, voxel size) tuple
@pytest.mark.parametrize(
    ('runs', 'expected'),
    [
        pytest.param([], False, id='no_runs'),
        pytest.param([('fmap1', 'j', TEST_ZOOMS)], False, id='one_run'),
        pytest.param([('fmap1', 'j', TEST_ZOOMS), ('fmap2', 'j', TEST_ZOOMS)], True, id='all_sdc'),
        pytest.param([('fmap1', 'j', TEST_ZOOMS), (None, 'j', TEST_ZOOMS)], False, id='mixed_sdc'),
        pytest.param(
            [(None, 'j', TEST_ZOOMS), (None, 'j', TEST_ZOOMS)], True, id='no_sdc_single_pe'
        ),
        pytest.param(
            [(None, 'j', TEST_ZOOMS), (None, 'j-', TEST_ZOOMS)], False, id='no_sdc_opposing_pe'
        ),
        pytest.param(
            [(None, None, TEST_ZOOMS), (None, None, TEST_ZOOMS)], True, id='no_sdc_no_pe'
        ),
        pytest.param(
            [(None, 'j', TEST_ZOOMS), (None, None, TEST_ZOOMS)], False, id='no_sdc_missing_pe'
        ),
        pytest.param(
            [('fmap1', 'j', TEST_ZOOMS), ('fmap2', 'j', (0.8, 0.8, 0.8))],
            False,
            id='mixed_resolution',
        ),
        pytest.param(
            [('fmap1', 'j', TEST_ZOOMS), ('fmap2', 'j', (2.4, 2.4, 3.0))],
            False,
            id='mixed_slice_thickness',
        ),
        pytest.param(
            [('fmap1', 'j', TEST_ZOOMS), ('fmap2', 'j', (2.4002, 2.4002, 2.4002))],
            True,
            id='negligible_zoom_difference',
        ),
    ],
)
def test_is_valid_bold_template(tmp_path: Path, runs, expected):
    bold_runs, estimator_map, pe_map = [], {}, {}
    for i, (fmapid, pe_dir, zooms) in enumerate(runs):
        path = str(tmp_path / f'run{i}.nii.gz')
        nb.Nifti1Image(np.zeros((2, 2, 2), dtype='uint8'), np.diag((*zooms, 1.0))).to_filename(
            path
        )
        bold_runs.append([path])
        estimator_map[path] = fmapid
        pe_map[path] = pe_dir

    layout = _MockLayout(pe_map)
    assert is_valid_bold_template(bold_runs, estimator_map, layout) is expected


@pytest.mark.parametrize(
    ('zooms', 'expected'),
    [
        pytest.param((2.4, 2.4, 2.4), (1.2, 1.2, 1.2), id='coarse'),
        pytest.param((4.0, 4.0, 4.0), (1.2, 1.2, 1.2), id='very_coarse'),
        pytest.param((3.0, 3.0, 4.0), (1.2, 1.2, 1.2), id='anisotropic'),
        pytest.param((1.0, 1.0, 2.0), (1.0, 1.0, 1.0), id='finest_axis_is_finer'),
        pytest.param((2.0, 2.0, 0.8), (0.8, 0.8, 0.8), id='finest_axis_is_much_finer'),
        pytest.param((1.2, 1.2, 1.2), None, id='at_target'),
        pytest.param((0.8, 0.8, 0.8), None, id='finer_than_target'),
    ],
)
def test_bold_template_zooms(tmp_path: Path, zooms, expected):
    bold_files = []
    for i in range(2):
        path = str(tmp_path / f'run{i}.nii.gz')
        nb.Nifti1Image(np.zeros((2, 2, 2), dtype='uint8'), np.diag((*zooms, 1.0))).to_filename(
            path
        )
        bold_files.append(path)

    assert bold_template_zooms(bold_files) == expected


def create_tsv(data: dict, out_file: Path) -> None:
    import pandas as pd

    pd.DataFrame(data).to_csv(out_file, index=False, sep='\t')


def create_sidecar(tsv_file: Path, units) -> None:
    import json

    out_file = tsv_file.with_suffix('.json')
    data = {'age': {'Units': units}}
    out_file.write_text(json.dumps(data))


age = {'age': [4, 4, 4]}
age_weeks = {'age_weeks': [4, 8, 12]}
age_months = {'age_months': [3, 6, 9]}
age_years = {'age_years': [1, 1, 2]}

participants = {'participant_id': ['sub-1', 'sub-2', 'sub-11']}
sessions = {'session_id': ['ses-1', 'ses-2', 'ses-3']}
scans = {
    'filename': [
        'dwi/sub-01_dwi.nii.gz',
        'anat/sub-01_T1w.nii.gz',
        'func/sub-01_task-rest_bold.nii.gz',
    ]
}


@pytest.mark.parametrize(
    ('idx_col', 'idx_val', 'data', 'units', 'expected'),
    [
        ('session_id', 'ses-1', age, 'months', 4),
        ('session_id', 'ses-1', age, 'weeks', 1),  # Convert from 4 weeks -> 1 month
        ('session_id', 'ses-1', age, 'week', 1),  # The BIDS example uses the singular
        ('session_id', 'ses-1', age, 'Month', 4),
        ('session_id', 'ses-2', age_weeks, False, 2),
        ('participant_id', 'sub-1', age_months, False, 3),
        ('participant_id', 'sub-11', age_years, False, 24),
        ('session_id', 'ses-3', {**age_months, **age}, False, 9),
        ('filename', r'^anat.*', age_months, False, 6),
    ],
)
def test_get_age_from_tsv(tmp_path, idx_col, idx_val, data, units, expected):
    tsv_file = tmp_path / 'test-age-parsing.tsv'

    if idx_col == 'participant_id':
        base = participants
    elif idx_col == 'session_id':
        base = sessions
    elif idx_col == 'filename':
        base = scans

    create_tsv({**base, **data}, tsv_file)
    if units:
        create_sidecar(tsv_file, units)

    res = _get_age_from_tsv(tsv_file, idx_col, idx_val)
    assert res == expected


def test_get_age_from_tsv_error(tmp_path):
    tsv_file = tmp_path / 'participants.tsv'

    create_tsv({**participants, **age}, tsv_file)
    with pytest.raises(FileNotFoundError):
        _get_age_from_tsv(tsv_file, 'participant_id', 'sub-1')


def test_get_age_from_tsv_warning(tmp_path):
    tsv_file = tmp_path / 'participants.tsv'
    dual_participants = {'participant_id': ['sub-1', 'sub-2', 'sub-2']}
    create_tsv({**dual_participants, **age_months}, tsv_file)

    with pytest.warns(UserWarning, match='Multiple matches for participant_id:sub-2'):
        _get_age_from_tsv(tsv_file, 'participant_id', 'sub-2')
