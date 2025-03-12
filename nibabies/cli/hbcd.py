"""
This script restructures workflow outputs to be ingested by the HBCD database.

WARNING: This alters the directories in place into a structure that the underlying software used
to create them will not recognize. Use with caution.

The following changes are made to the outputs:

1. FreeSurfer output is changed to follow the BIDS hierarchy:
    freesurfer/
        sub-<subject>
            ses-<session>/
                mri/
                surf/
                ...

2. MCRIBS output is changed to follow the BIDS hierarchy:
    mcribs/
        sub-<subject>
            ses-<session>/
                SurfReconDeformable/
                TissueSegDrawEM/
                ...

3. Symbolic links are followed and copied.
"""

import argparse
import shutil
from pathlib import Path


def _parser():
    from functools import partial

    from .parser import _path_exists

    parser = argparse.ArgumentParser(description='Prepare outputs for HBCD database ingestion.')

    PathExists = partial(_path_exists, parser=parser)

    parser.add_argument(
        'deriv_dir', type=PathExists, help='Path to the BIDS derivatives directory'
    )
    parser.add_argument(
        '--fs-dir',
        type=PathExists,
        help='Path to the FreeSurfer directory. If not provided, will look for in derivatives.',
    )
    parser.add_argument(
        '--mcribs-dir',
        type=PathExists,
        help='Path to the MCRIBS directory. If not provided, will look for in the derivatives.',
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
        print(sid)
        try:
            subject, session = sid.name.split('_', 1)
        except ValueError:
            print(f'Could not split {sid} into subject and session')
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


def main(argv=None):
    """Entry point `nibabies-hbcd`."""
    parser = _parser()
    args = parser.parse_args(argv)

    derivatives = args.deriv_dir
    fs_dir = args.fs_dir
    mcribs_dir = args.mcribs_dir

    if fs_dir is None:
        fs_dir = derivatives / 'nibabies' / 'sourcedata' / 'freesurfer'
        if not fs_dir.exists():
            raise FileNotFoundError(
                f'Could not find FreeSurfer directory at {fs_dir} - use `--fs-dir`.'
            )

    if mcribs_dir is None:
        mcribs_dir = derivatives / 'nibabies' / 'sourcedata' / 'mcribs'
        if not fs_dir.exists():
            raise FileNotFoundError(
                f'Could not find FreeSurfer directory at {mcribs_dir} - use `--mcribs-dir`.'
            )

    restructure(fs_dir)
    restructure(mcribs_dir)
