from __future__ import annotations

import logging
import shutil
from pathlib import Path

from nipost.bids import (
    collect_derivatives,
    collect_fieldmaps,
    load_spec,
    sanitize_fieldmap_id,
    sanitize_space,
)

from nibabies.data import load

TYPE_CHECKING = False
if TYPE_CHECKING:
    from typing import NotRequired, TypedDict

    class AggregateCache(TypedDict):
        template2anat_xfm: list[str]
        run2template_xfms: NotRequired[list[str]]
        boldref_template: NotRequired[str | None]


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
        image, coregistration transform and surface they provide, and for each space
        they have transforms for. Forward and reverse transforms for a space are never
        combined across datasets.
        Each derivative a later dataset replaces is logged at DEBUG level.
    subject_id
        Subject label, without ``sub-``.
    std_spaces
        Standard spaces to find transforms for, in TemplateFlow form
        (e.g., ``MNIInfant:cohort-1``). Transforms to ``fsnative`` are always searched.
    session_id
        Session label, without ``ses-``.

    Returns
    -------
    dict
        ``images``, ``coreg`` and ``surfaces`` map query names in ``anat_spec.yml``
        to paths. ``transforms`` maps each space with precomputed transforms to its
        ``forward`` and/or ``reverse`` transform.

    Raises
    ------
    ValueError
        If a query matches more than one file in a dataset.
    """
    entities = {'subject': subject_id}
    if session_id:
        entities['session'] = session_id
    # Filenames use the BIDS form of cohort spaces; workflows use the TemplateFlow form
    spaces = {sanitize_space(space): space for space in [*std_spaces, 'fsnative']}
    params = {'space': list(spaces)}
    spec = load_spec(load('anat_spec.yml'))

    deriv_cache = {'images': {}, 'coreg': {}, 'surfaces': {}, 'transforms': {}}
    for deriv_dir in derivatives:
        collected = collect_derivatives(deriv_dir, spec=spec, entities=entities, params=params)
        _merge(deriv_cache['images'], collected['images'], 'image', deriv_dir)
        _merge(deriv_cache['coreg'], collected['coreg'], 'transform', deriv_dir)
        _merge(deriv_cache['surfaces'], collected['surfaces'], 'surface', deriv_dir)
        # Replace whole spaces, so a transform pair always comes from one dataset.
        # Spaces this dataset has no transforms for are empty and must not replace others.
        transforms = {
            spaces[space]: xfms for space, xfms in collected['transforms'].items() if xfms
        }
        _merge(deriv_cache['transforms'], transforms, 'transforms for', deriv_dir)

    return deriv_cache


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


def collect_fmap_derivatives(
    derivatives: list[Path],
    subject_id: str,
    session_id: str | list[str] | None = None,
) -> dict:
    """Gather precomputed fieldmaps from one or more datasets.

    Parameters
    ----------
    derivatives
        Derivatives datasets to search. A later dataset's files for a fieldmap ID replace
        any earlier dataset's files for that ID; files for one ID are never combined
        across datasets. Each replacement is logged at DEBUG level.
    subject_id
        Subject label, without ``sub-``.
    session_id
        Session label(s), without ``ses-``. Fieldmap IDs may repeat across sessions,
        so only the given sessions are searched. By default, only fieldmaps without
        a session are.

    Returns
    -------
    dict
        Precomputed fieldmaps by fieldmap ID, each mapping the queries in ``fmap_spec.yml``
        (``fieldmap``, ``coeffs`` and ``magnitude``) that found files to those files.
        ``coeffs`` is a list.
    """
    spec = load_spec(load('fmap_spec.yml'))
    entities = {'subject': subject_id, 'session': session_id}
    merged = {}
    for deriv_dir in derivatives:
        found = collect_fieldmaps(deriv_dir, entities=entities, spec=spec)
        # Keep only the IDs with files
        fmaps = {fmapid: files for fmapid, files in found['fieldmaps'].items() if files}
        _merge(merged, fmaps, 'fieldmap', deriv_dir)
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
