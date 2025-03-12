"""
This script restructures workflow outputs to be ingested by the HBCD database.

The following changes are made to the outputs:

- FreeSurfer output is changed to follow the BIDS hierarchy:

    freesurfer/
        sub-<subject>
            ses-<session>/
                mri/
                surf/
                ...

- MCRIBS output is changed to follow the BIDS hierarchy:

    mcribs/
        sub-<subject>
            ses-<session>/
                SurfReconDeformable/
                TissueSegDrawEM/
                ...

- Symbolic links are replaced with the files they point to.

WARNING: This alters the directories in place into a structure that the
underlying software used to create them will not recognize. Use with caution.
"""

import argparse
import shutil
from pathlib import Path


def _parser():
    from functools import partial

    from .parser import _path_exists

    parser = argparse.ArgumentParser(
        prog='nibabies-hbcd',
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    PathExists = partial(_path_exists, parser=parser)

    parser.add_argument(
        '--fs',
        type=PathExists,
        help='Path to the FreeSurfer output directory',
    )
    parser.add_argument(
        '--mcribs',
        type=PathExists,
        help='Path to the MCRIBS output directory.',
    )
    return parser


def copy_symlinks(directory: Path):
    for fl in directory.rglob('*'):
        if fl.is_symlink():
            target = fl.resolve()
            print(f'Found symlink {fl} pointing to {target}')
            fl.unlink()
            shutil.copy2(target, fl)


def restructure(directory: Path):
    """Change the structure of a directory in place to resemble BIDS hierarchy."""
    for sid in directory.glob('sub-*'):
        try:
            subject, session = sid.name.split('_', 1)
            print(sid)
        except ValueError:
            continue

        if not subject.startswith('sub-'):
            raise AttributeError(f'Incorrect subject ID {subject}')
        if not session.startswith('ses-'):
            raise AttributeError(f'Incorrect session ID {session}')

        # First traverse and ensure no symbolic links are present
        copy_symlinks(sid)

        target_directory = directory / subject / session
        print(f'Making target directory {target_directory}')
        target_directory.mkdir(parents=True, exist_ok=True)

        print(f'Copying {sid} to {target_directory}')
        shutil.copytree(sid, target_directory, dirs_exist_ok=True)
        shutil.rmtree(sid)
    print(f'Completed restructuring {directory}')


def main(argv=None):
    """Entry point `nibabies-hbcd`."""
    parser = _parser()
    pargs = parser.parse_args(argv)

    fs = pargs.fs
    if fs is None:
        print('FreeSurfer directory not provided. Skipping')
    else:
        restructure(fs)

    mcribs = pargs.mcribs
    if mcribs is None:
        print('MCRIBS directory not provided. Skipping')
    else:
        restructure(mcribs)
