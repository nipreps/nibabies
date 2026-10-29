from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path

from bids.layout import BIDSLayout
from nipost.bids import collect_derivatives, load_spec, sanitize_fieldmap_id
from niworkflows.data import load as nwf_load

from nibabies.data import load

TYPE_CHECKING = False
if TYPE_CHECKING:
    from typing import NotRequired, TypedDict

    class AggregateCache(TypedDict):
        template2anat_xfm: list[str]
        run2template_xfms: NotRequired[list[str]]
        boldref_template: NotRequired[str | None]


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


def _merge(merged: dict, found: dict, kind: str, deriv_dir) -> None:
    """Update ``merged`` with ``found``, logging each entry that is replaced."""
    for name, value in found.items():
        if name in merged:
            logging.getLogger('nipype.utils').debug(
                f'Precomputed {kind} {name} found in {deriv_dir}, replacing {merged[name]}'
            )
        merged[name] = value


def collect_anat_derivatives(
    derivatives: list[Path],
    subject_id: str,
    std_spaces: list,
    session_id: str | None,
) -> dict:
    """Gather precomputed anatomical derivatives from one or more datasets.

    Parameters
    ----------
    derivatives
        Derivatives datasets to search. Later datasets take precedence for each
        derivative and each space's transforms they provide.
        Each derivative a later dataset replaces is logged at DEBUG level.
    subject_id
        Subject label, without ``sub-``.
    std_spaces
        Standard spaces to find transforms for.
    session_id
        Session label, without ``ses-``.

    Returns
    -------
    dict
        Paths to precomputed derivatives, as returned by
        :func:`collect_anatomical_derivatives`.
    """
    merged = {}
    transforms = {}
    for deriv_dir in derivatives:
        found = collect_anatomical_derivatives(
            derivatives_dir=deriv_dir,
            subject_id=subject_id,
            session_id=session_id,
            std_spaces=std_spaces,
        )
        _merge(transforms, found.pop('transforms', {}), 'transform', deriv_dir)
        _merge(merged, found, 'derivative', deriv_dir)
    return {**merged, 'transforms': transforms}


def collect_func_derivatives(
    derivatives: list[Path],
    entities: dict,
    fieldmap_id: str | None = None,
) -> dict:
    """Gather precomputed functional derivatives from one or more datasets.

    Parameters
    ----------
    derivatives
        Derivatives datasets to search. Later datasets take precedence for each
        derivative and each transform they provide.
        Each derivative a later dataset replaces is logged at DEBUG level.
    entities
        BIDS entities of the BOLD series.
    fieldmap_id
        Fieldmap estimator ID, used to find the BOLD-to-fieldmap transform.
        Without one, no BOLD-to-fieldmap transform is collected.

    Returns
    -------
    dict
        ``boldrefs`` and ``transforms`` map query names in ``func_spec.yml`` to paths.
        ``run2fmap``, from the spec's ``fmap_transforms`` group for the requested
        fieldmap ID, is merged into ``transforms``.

    Raises
    ------
    ValueError
        If a query matches more than one file in a dataset.
    """
    spec = load_spec(load('func_spec.yml'))
    fmapid = sanitize_fieldmap_id(fieldmap_id) if fieldmap_id else None
    merged = {'boldrefs': {}, 'transforms': {}}
    for deriv_dir in derivatives:
        found = collect_derivatives(
            deriv_dir, spec=spec, entities=entities, params={'fmapid': fmapid}
        )
        # A run has at most one estimator, so its transform joins the others
        transforms = {**found['transforms'], **found['fmap_transforms'].get(fmapid, {})}
        _merge(merged['boldrefs'], found['boldrefs'], 'boldref', deriv_dir)
        _merge(merged['transforms'], transforms, 'transform', deriv_dir)
    return merged


def aggregate_coreg_precomputed(caches: list[dict], level: str) -> AggregateCache:
    """Aggregate coregistration precomputed inputs from per-run caches."""

    def get_xfm(cache, key):
        return cache.get('transforms', {}).get(key)

    precomputed: AggregateCache = {
        'template2anat_xfm': [get_xfm(c, f'{level}2anat') for c in caches]
    }
    if level != 'run':
        precomputed['run2template_xfms'] = [get_xfm(c, f'run2{level}') for c in caches]
        precomputed['boldref_template'] = next(
            (ref for c in caches if (ref := c.get('boldrefs', {}).get(level))),
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
