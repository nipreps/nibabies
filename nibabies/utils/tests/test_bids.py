from pathlib import Path

import pytest

from nibabies.utils.bids import _get_age_from_tsv, is_valid_bold_template


class _MockLayout:
    """Minimal BIDSLayout stub for bold_coreg_compat tests."""

    def __init__(self, pe_map: dict):
        self._pe_map = pe_map

    def get_metadata(self, f):
        pe = self._pe_map.get(f)
        return {'PhaseEncodingDirection': pe} if pe is not None else {}


@pytest.mark.parametrize(
    ('bold_runs', 'estimator_map', 'pe_map', 'expected'),
    [
        ([], {}, {}, False),
        ([['a.nii']], {}, {'a.nii': 'j'}, False),
        ([['a.nii'], ['b.nii']], {'a.nii': 'fmap1', 'b.nii': 'fmap2'}, {}, True),
        ([['a.nii'], ['b.nii']], {'a.nii': 'fmap1'}, {'a.nii': 'j', 'b.nii': 'j'}, False),
        ([['a.nii'], ['b.nii']], {}, {'a.nii': 'j', 'b.nii': 'j'}, True),
        ([['a.nii'], ['b.nii']], {}, {'a.nii': 'j', 'b.nii': 'j-'}, False),
        ([['a.nii'], ['b.nii']], {}, {}, True),
        ([['a.nii'], ['b.nii']], {}, {'a.nii': 'j'}, False),
    ],
    ids=[
        'false-no_runs',
        'false-one_run',
        'true-all_sdc',
        'false-mixed_sdc',
        'true-no_sdc_single_pe',
        'false-no_sdc_opposing_pe',
        'true-no_sdc_no_pe',
        'false-no_sdc_missing_pe',
    ],
)
def test_is_valid_bold_template(bold_runs, estimator_map, pe_map, expected):
    layout = _MockLayout(pe_map)
    assert is_valid_bold_template(bold_runs, estimator_map, layout) is expected


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
