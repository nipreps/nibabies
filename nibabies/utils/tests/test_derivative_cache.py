import logging
from pathlib import Path

import pytest

from nibabies.utils import derivatives
from nibabies.utils.testing import write_derivatives

ENTITIES = {
    'subject': '01',
    'task': 'rest',
    'run': 1,
    'suffix': 'bold',
    'extension': '.nii.gz',
}


def _files(root: Path, datatype: str = 'func') -> list[str]:
    return sorted(str(p) for p in root.rglob(f'{datatype}/*') if p.is_file())


def _entities(session):
    return {**ENTITIES, 'session': session} if session else ENTITIES


BOLDREFS = [
    'hmc_boldref',
    'hmc_boldref_legacy',
    'run_boldref',
    'run_boldref_legacy',
    'session_boldref',
    'subject_boldref',
]
TRANSFORMS = [
    'hmc',
    'hmc_legacy',
    'run2anat',
    'run2anat_legacy',
    'run2fmap',
    'run2fmap_legacy',
    'run2session',
    'run2subject',
    'session2anat',
    'subject2anat',
]
# Groups collected under another name
KEYS = {'run2session': 'run2template', 'run2subject': 'run2template'}
# The run-level reference is collected as ``coreg_boldref``, aliased as ``run_boldref``
BOLDREF_KEYS = {'run_boldref': ('coreg_boldref', 'run_boldref')}


def _key(group: str) -> str:
    return KEYS.get(group, group.removesuffix('_legacy'))


@pytest.mark.parametrize('session', [None, 'A'])
@pytest.mark.parametrize('group', BOLDREFS)
def test_boldref_found_as_str(tmp_path: Path, group: str, session):
    """Generate a single boldref and verify it's found.

    `<x>_legacy` is defined in data/tests/derivatives.yml to generate files,
    but the collector loads `<x>` and `<x>_legacy` as `<x>`.
    """
    root = write_derivatives(tmp_path / 'deriv', [group], run=1, session=session)
    [found] = _files(root)

    key = group.removesuffix('_legacy')
    derivs = derivatives.collect_func_derivatives(
        [root], _entities(session), fieldmap_id='auto_00000'
    )
    assert derivs == {
        **dict.fromkeys(BOLDREF_KEYS.get(key, (key,)), found),
        'transforms': {},
    }


@pytest.mark.parametrize('session', [None, 'A'])
@pytest.mark.parametrize('group', TRANSFORMS)
def test_transform_found_as_str(tmp_path: Path, group: str, session):
    """Generate a single transform and verify it's found.

    See above RE legacy derivatives. ``KEYS`` lists the other groups
    that are collected under another name.
    """
    root = write_derivatives(tmp_path / 'deriv', [group], run=1, session=session)
    [found] = _files(root)

    key = _key(group)
    derivs = derivatives.collect_func_derivatives(
        [root], _entities(session), fieldmap_id='auto_00000'
    )
    assert derivs == {'transforms': {key: found}}

    # Collector function normalizes fieldmap IDs internally
    if key == 'run2fmap':
        assert '_to-auto00000_' in found


@pytest.mark.parametrize('session', [None, 'A'])
def test_defaults_collected_once_each(tmp_path: Path, session):
    """Default (None) includes loads all keys."""
    root = write_derivatives(tmp_path / 'deriv', None, run=1, session=session)

    derivs = derivatives.collect_func_derivatives(
        [root], _entities(session), fieldmap_id='auto_00000'
    )
    transforms = derivs.pop('transforms')
    assert derivs.keys() == {
        'hmc_boldref',
        'coreg_boldref',
        'run_boldref',
        'session_boldref',
        'subject_boldref',
    }
    assert transforms.keys() == {_key(group) for group in TRANSFORMS}
    assert {*derivs.values(), *transforms.values()} == set(_files(root))


@pytest.mark.xfail(reason='existing collector returns both the current and the legacy file')
@pytest.mark.parametrize('group', ['hmc_boldref', 'run_boldref', 'hmc', 'run2anat', 'run2fmap'])
def test_current_name_preferred_over_legacy(tmp_path: Path, group: str):
    """Generate current and legacy files and verify only current is returned."""
    # Each current-name group is collected under its own name
    [current] = _files(write_derivatives(tmp_path / 'current', [group], run=1))
    root = write_derivatives(tmp_path / 'deriv', [group, f'{group}_legacy'], run=1)

    derivs = derivatives.collect_func_derivatives([root], ENTITIES, fieldmap_id='auto_00000')
    found = derivs['transforms'][group] if group in derivs['transforms'] else derivs[group]
    assert found == str(root / 'sub-01' / 'func' / Path(current).name)


@pytest.mark.xfail(reason='existing run2fmap query matches any transform from the run')
def test_no_run2fmap_without_fieldmap(tmp_path: Path):
    root = write_derivatives(tmp_path / 'deriv', ['run2anat'], run=1)

    derivs = derivatives.collect_func_derivatives([root], ENTITIES)
    assert 'run2fmap' not in derivs['transforms']


def _find(root: Path, pattern: str) -> str:
    [path] = root.rglob(pattern)
    return str(path)


def test_collect_func_derivatives_merges_per_key(tmp_path: Path, caplog):
    first = write_derivatives(tmp_path / 'first', ['hmc_boldref', 'hmc', 'run2anat'], run=1)
    second = write_derivatives(tmp_path / 'second', ['run_boldref', 'hmc'], run=1)

    with caplog.at_level(logging.DEBUG, logger='nipype.utils'):
        derivs = derivatives.collect_func_derivatives(
            [first, second], ENTITIES, fieldmap_id='auto_00000'
        )

    # hmc is overridden
    run_boldref = _find(second, '*_space-run_boldref.nii.gz')
    assert derivs == {
        'hmc_boldref': _find(first, '*_desc-hmc_boldref.nii.gz'),
        'coreg_boldref': run_boldref,
        'run_boldref': run_boldref,
        'transforms': {
            'hmc': _find(second, '*_desc-hmc_xfm.txt'),
            'run2anat': _find(first, '*_to-T1w_*xfm.txt'),
        },
    }

    # One DEBUG message about HMC
    [record] = [r for r in caplog.records if 'replacing' in r.getMessage()]
    assert record.levelno == logging.DEBUG
    first_path = _find(first, '*_desc-hmc_xfm.txt')
    message = f'Precomputed transform hmc found in {second}, replacing {first_path}'
    assert record.getMessage() == message


def test_collect_func_derivatives_empty_later_dataset(tmp_path: Path, caplog):
    """Empty datasets do not clobber populated keys."""
    first = write_derivatives(tmp_path / 'first', ['run_boldref', 'hmc'], run=1)
    empty = write_derivatives(tmp_path / 'empty', [])

    with caplog.at_level(logging.DEBUG, logger='nipype.utils'):
        derivs = derivatives.collect_func_derivatives([first, empty], ENTITIES)
    assert derivs == derivatives.collect_func_derivatives([first], ENTITIES)
    assert derivs['transforms']
    assert not [r for r in caplog.records if 'replacing' in r.getMessage()]


def test_collect_func_derivatives_nothing():
    assert derivatives.collect_func_derivatives([], ENTITIES) == {'transforms': {}}


def _run_caches(tmp_path: Path, shared: list[str], per_run: list[str], session='A'):
    """Collect caches for runs 1 and 2, with ``per_run`` derivatives in a dataset per run."""
    shared_dir = write_derivatives(tmp_path / 'shared', shared, session=session)
    run_dirs = [
        write_derivatives(tmp_path / f'run-{run}', per_run, session=session, run=run)
        for run in (1, 2)
    ]
    return [
        derivatives.collect_func_derivatives(
            [shared_dir, run_dir], {**_entities(session), 'run': run}
        )
        for run, run_dir in enumerate(run_dirs, start=1)
    ]


def test_aggregate_coreg_precomputed_run(tmp_path: Path):
    caches = _run_caches(tmp_path, ['session2anat'], ['run2anat', 'run2session'])

    assert derivatives.aggregate_coreg_precomputed(caches, 'run') == {
        'template2anat_xfm': [c['transforms']['run2anat'] for c in caches],
    }


def test_aggregate_coreg_precomputed_group(tmp_path: Path):
    caches = _run_caches(tmp_path, ['session_boldref', 'session2anat'], ['run2session'])
    session2anat = caches[0]['transforms']['session2anat']

    assert derivatives.aggregate_coreg_precomputed(caches, 'session') == {
        'template2anat_xfm': [session2anat, session2anat],
        'run2template_xfms': [c['transforms']['run2template'] for c in caches],
        'boldref_template': caches[0]['session_boldref'],
    }
    assert caches[0]['transforms']['run2template'] != caches[1]['transforms']['run2template']
