from __future__ import annotations

import json
import re
import shutil
from collections import defaultdict
from pathlib import Path

from bids.layout import BIDSLayout
from niworkflows.data import load as nwf_load

from nibabies.data import load


def collect_anatomical_derivatives(
    derivatives_dir: Path | str,
    subject_id: str,
    std_spaces: list,
    session_id: str | None,
    spec: dict | None = None,
    patterns: list | None = None,
):
    """
    Collect outputs from across processing stages.

    Potential files:
    - T1w preproc
    - T2w preproc
    - T1w mask
    - T2w mask


    """

    if spec is None or patterns is None:
        _spec, _patterns = tuple(json.loads(load('io_spec_anat.json').read_text()).values())

        if spec is None:
            spec = _spec
        if patterns is None:
            patterns = _patterns

    deriv_config = nwf_load('nipreps.json')
    layout = BIDSLayout(derivatives_dir, config=deriv_config, validate=False)
    derivs_cache = {}

    base_qry = {
        'subject': subject_id,
    }
    if session_id is not None:
        base_qry['session'] = session_id

    for key, qry in spec['baseline'].items():
        qry.update(base_qry)
        item = layout.get(return_type='filename', **qry)
        if not item:
            continue

        derivs_cache[key] = item[0] if len(item) == 1 else item

    for key, qry in spec['coreg'].items():  # T1w->T2w, T2w->T1w
        qry.update(base_qry)
        item = layout.get(return_type='filename', **qry)
        if not item:
            continue
        derivs_cache[key] = item[0] if len(item) == 1 else item

    transforms = derivs_cache.setdefault('transforms', {})
    for _space in std_spaces:
        space = _space.replace(':cohort-', '+')
        for key, qry in spec['transforms'].items():
            qry = qry.copy()
            qry.update(base_qry)
            qry['from'] = qry['from'] or space
            qry['to'] = qry['to'] or space
            item = layout.get(return_type='filename', **qry)
            if not item:
                continue
            transforms.setdefault(_space, {})[key] = item[0] if len(item) == 1 else item

    for key, qry in spec['surfaces'].items():
        qry.update(base_qry)
        item = layout.get(return_type='filename', **qry)
        if not item or len(item) != 2:
            continue

        derivs_cache[key] = sorted(item)

    return derivs_cache


def collect_functional_derivatives(
    derivatives_dir: Path,
    entities: dict,
    fieldmap_id: str | None,
    spec: dict | None = None,
    patterns: list[str] | None = None,
):
    """Gather existing derivatives and compose a cache."""
    if spec is None or patterns is None:
        _spec, _patterns = tuple(
            json.loads(load.readable('io_spec_func.json').read_text()).values()
        )

        if spec is None:
            spec = _spec
        if patterns is None:
            patterns = _patterns

    derivs_cache = defaultdict(list, {})
    deriv_config = nwf_load('nipreps.json')
    layout = BIDSLayout(derivatives_dir, config=deriv_config, validate=False)
    derivatives_dir = Path(derivatives_dir)

    from nibabies.utils.bids import GROUP_DISMISS_ENTITIES

    # Session- and subject-level templates are written once, dropping the
    # run-varying entities; relax the query by those entities to find them.
    # TODO: Move towards expected filenames for these group-level files
    level_relax = {
        'run': (),
        'session': GROUP_DISMISS_ENTITIES,
        'subject': (*GROUP_DISMISS_ENTITIES, 'session'),
    }

    def _query(q, relax=()):
        query = {k: v for k, v in {**entities, **q}.items() if k not in relax}
        item = layout.get(return_type='filename', **query)
        if not item:
            return None
        return item[0] if len(item) == 1 else item

    # BOLD references, one group per coregistration level. Legacy naming
    # (space-less/desc-coreg references, from-boldref transforms) is captured by
    # the list-valued queries in ``io_spec_func.json``.
    for level, relax in level_relax.items():
        for name, q in spec.get(level, {}).items():
            if q.get('suffix') != 'boldref':
                continue
            item = _query(q, relax)
            if not item:
                continue
            if name == 'hmc':
                key = 'hmc_boldref'
            elif level == 'run':
                key = 'coreg_boldref'
            else:
                key = f'{level}_boldref'
            derivs_cache[key] = item

    # run_boldref is an alias for the run-native coregistration reference
    if 'coreg_boldref' in derivs_cache and 'run_boldref' not in derivs_cache:
        derivs_cache['run_boldref'] = derivs_cache['coreg_boldref']

    # Per-run transforms. Transform extensions/suffixes will not match the provided
    #   entities (e.g., ".txt" vs ".nii.gz", "xfm" vs "bold"); the queries override them.
    transforms_cache = {}
    for xfm, q in spec['transforms'].items():
        query = {**entities, **q}
        if xfm == 'run2fmap' and fieldmap_id:
            # fieldmaps have non-alphanumeric characters removed from their IDs in filenames
            query['to'] = re.sub(r'[^a-zA-Z0-9]', '', fieldmap_id)
        item = layout.get(return_type='filename', **query)
        if not item:
            continue
        transforms_cache[xfm] = item[0] if len(item) == 1 else item

    # Session-/subject-level template->anat transforms are written once
    for level in ('session', 'subject'):
        xfm_q = spec.get(level, {}).get('xfm')
        if not xfm_q:
            continue
        item = _query(xfm_q, level_relax[level])
        if item:
            transforms_cache[f'{level}2anat'] = item

    derivs_cache['transforms'] = transforms_cache
    return derivs_cache


def aggregate_coreg_precomputed(caches: list[dict], level: str) -> dict:
    """Aggregate coregistration precomputed inputs from per-run caches."""

    def get_xfm(cache, key):
        return cache.get('transforms', {}).get(key)

    precomputed = {'template2anat_xfm': [get_xfm(c, f'{level}2anat') for c in caches]}
    if level != 'run':
        precomputed['run2template_xfms'] = [get_xfm(c, 'run2template') for c in caches]
        precomputed['boldref_template'] = next(
            (c[f'{level}_boldref'] for c in caches if c.get(f'{level}_boldref')),
            None,
        )
    return precomputed


def copy_derivatives(
    derivs: dict,
    outdir: Path,
    modality: str,
    subject_id: str,
    session_id: str | None = None,
    config_hash: str | None = None,
) -> None:
    """
    Creates a copy of any found derivatives into output directory.

    Attempts to preserve file metadata to distinguish from generated files.
    """
    out_levels = [subject_id, modality]
    if session_id:
        out_levels.insert(1, session_id)

    outpath = outdir.joinpath(*out_levels)
    outpath.mkdir(parents=True, exist_ok=True)

    # Flatten one level so nested caches (e.g. the ``transforms`` sub-dict) are copied too
    candidates = []
    for value in derivs.values():
        if isinstance(value, dict):
            candidates.extend(value.values())
        else:
            candidates.append(value)

    for deriv in candidates:
        # Skip empty, lists
        if not isinstance(deriv, str):
            continue
        deriv = Path(deriv)
        outname = deriv.name

        if config_hash:
            ents = outname.split('_')
            if any(ent.startswith('hash-') for ent in ents):
                # Avoid adding another hash
                pass
            else:
                idx = 2 if ents[1].startswith('ses-') else 1
                ents.insert(idx, f'hash-{config_hash}')
                outname = '_'.join(ents)

        shutil.copy2(deriv, outpath / outname)
        json = deriv.parent / (outname.split('.')[0] + '.json')
        if json.exists():
            shutil.copy2(json, outpath / json.name)
