from pathlib import Path

import pytest

from nibabies.utils.bids import _get_age_from_tsv


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


@pytest.mark.parametrize(
    ('idx_col', 'idx_val', 'data', 'sidecar', 'expected'),
    [
        ('session_id', 'x1', age, False, None),
        ('session_id', 'x1', age, 'months', 4),
        ('session_id', 'x1', age, 'weeks', 1),  # Convert from 4 weeks -> 1 month
        ('session_id', 'x1', age, ['months', 'weeks'], None),
        ('session_id', 'x2', age_weeks, False, 2),
        ('participant_id', 'x1', age_months, False, 3),
        ('participant_id', 'x3', age_years, False, 24),
        ('session_id', 'x3', {**age_months, **age}, False, 9),
        (None, None, age_months, False, 3),
    ],
)
def test_get_age_from_tsv(tmp_path, idx_col, idx_val, data, sidecar, expected):
    tsv_file = tmp_path / 'test-age-parsing.tsv'
    base = {}
    if idx_col is not None:
        base[idx_col] = ['x1', 'x2', 'x3']
    create_tsv({**base, **data}, tsv_file)

    if sidecar:
        create_sidecar(tsv_file, sidecar)

    res = _get_age_from_tsv(tsv_file, idx_col, idx_val)
    assert res == expected
