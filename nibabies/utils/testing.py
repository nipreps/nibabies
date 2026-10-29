"""Helpers for writing test datasets."""

import re
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import yaml
from niworkflows.utils.testing import generate_bids_skeleton

from nibabies.data import load

_PLACEHOLDER = re.compile(r'\{(\w+)\}')


def _fill(value, params):
    if isinstance(value, str) and (match := _PLACEHOLDER.fullmatch(value)):
        return params[match.group(1)]
    return value


def deriv_skeleton(
    include=None, *, subject='01', session=None, task='rest', run=None, fmapid='auto00000'
):
    """Describe a skeleton of NiBabies functional derivatives.

    Parameters
    ----------
    include : :obj:`list` of :obj:`str` or :obj:`None`
        Groups from ``nibabies/data/tests/derivatives.yml`` to write.
        By default, the groups listed under ``defaults``.
    subject : :obj:`str`
        Subject label, without ``sub-``.
    session, task, run, fmapid
        Values for the ``{name}`` placeholders in the groups.
        Entities whose value is ``None`` are left out.

    Returns
    -------
    :obj:`dict`
        A dataset description for :func:`niworkflows.utils.testing.generate_bids_skeleton`.

    """
    spec = yaml.safe_load(load.readable('tests/derivatives.yml').read_text())
    params = {'session': session, 'task': task, 'run': run, 'fmapid': fmapid}

    sessions = {}
    for name in spec['defaults'] if include is None else include:
        for file in spec['groups'][name]:
            entities = {
                key: interp
                for key, value in file.items()
                if (interp := _fill(value, params)) is not None
            }
            datatype = entities.pop('datatype')
            ses = entities.pop('session', None)
            sessions.setdefault(ses, {}).setdefault(datatype, []).append(entities)

    return {
        'dataset_description': spec['dataset_description'],
        subject: [{'session': ses, **datatypes} for ses, datatypes in sessions.items()],
    }


def write_derivatives(path, include=None, **kwargs) -> Path:
    """Write a skeleton of NiBabies functional derivatives.

    Parameters
    ----------
    path : :obj:`os.PathLike`
        Root of the dataset. If it exists, the files are added to it,
        so that a dataset can hold several runs.
    include, **kwargs
        Passed to :func:`deriv_skeleton`.

    Returns
    -------
    :obj:`~pathlib.Path`
        ``path``.

    """
    path = Path(path)
    skeleton = deriv_skeleton(include, **kwargs)
    if not path.exists():
        generate_bids_skeleton(path, skeleton)
        return path

    # generate_bids_skeleton() only writes new datasets
    with TemporaryDirectory() as tmpdir:
        generate_bids_skeleton(Path(tmpdir) / 'deriv', skeleton)
        shutil.copytree(Path(tmpdir) / 'deriv', path, dirs_exist_ok=True)
    return path
