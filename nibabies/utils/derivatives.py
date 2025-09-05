from __future__ import annotations

import json
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

    # search for both boldrefs
    for key, qry in spec['baseline'].items():
        query = {**qry, **entities}
        item = layout.get(return_type='filename', **query)
        if not item:
            continue
        derivs_cache[f'{key}_boldref'] = item[0] if len(item) == 1 else item

    for xfm, qry in spec['transforms'].items():
        query = {**qry, **entities}
        if xfm == 'boldref2fmap':
            query['to'] = fieldmap_id
        item = layout.get(return_type='filename', **query)
        if not item:
            continue
        derivs_cache[xfm] = item[0] if len(item) == 1 else item
    return derivs_cache


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

    for deriv in derivs.values():
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
