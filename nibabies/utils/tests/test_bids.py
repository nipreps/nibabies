from __future__ import annotations

import json
import typing as ty
from pathlib import Path

import pytest

from nibabies.utils import bids


def _create_nifti(filename: str) -> str:
    import nibabel as nb
    import numpy as np

    data = np.zeros((4, 4, 4), dtype='int8')
    nb.Nifti1Image(data, np.eye(4)).to_filename(filename)
    return filename


def _create_bids_dir(root_path: Path):
    if not root_path.exists():
        root_path.mkdir()
    anat_dir = root_path / 'sub-01' / 'anat'
    anat_dir.mkdir(parents=True)
    _create_nifti(str(anat_dir / 'sub-01_T1w.nii.gz'))
    _create_nifti(str(anat_dir / 'sub-01_T2w.nii.gz'))


def _create_bids_derivs(
    root_path: Path,
    *,
    t1w_mask: bool = False,
    t1w_aseg: bool = False,
    t2w_mask: bool = False,
    t2w_aseg: bool = False,
):
    if not root_path.exists():
        root_path.mkdir()
    (root_path / 'dataset_description.json').write_text(
        json.dumps(
            {'Name': 'Derivatives Test', 'BIDSVersion': '1.8.0', 'DatasetType': 'derivative'}
        )
    )
    anat_dir = root_path / 'sub-01' / 'anat'
    anat_dir.mkdir(parents=True)

    def _create_deriv(name: str, modality: ty.Literal['t1w', 't2w']):
        if modality == 't1w':
            reference = 'sub-01/anat/sub-01_T1w.nii.gz'
        elif modality == 't2w':
            reference = 'sub-01/anat/sub-01_T2w.nii.gz'

        _create_nifti(str((anat_dir / name).with_suffix('.nii.gz')))
        (anat_dir / name).with_suffix('.json').write_text(
            json.dumps({'SpatialReference': reference})
        )

    if t1w_mask:
        _create_deriv('sub-01_space-T1w_desc-brain_mask', 't1w')
    if t1w_aseg:
        _create_deriv('sub-01_space-T1w_desc-aseg_dseg', 't1w')
    if t2w_mask:
        _create_deriv('sub-01_space-T2w_desc-brain_mask', 't2w')
    if t2w_aseg:
        _create_deriv('sub-01_space-T2w_desc-aseg_dseg', 't2w')


@pytest.mark.parametrize(
    't1w_mask,t1w_aseg,t2w_mask,t2w_aseg,mask,aseg',
    [
        (True, True, False, False, 't1w_mask', 't1w_aseg'),
        (True, True, True, True, 't1w_mask', 't1w_aseg'),
        (False, False, True, True, 't2w_mask', 't2w_aseg'),
        (True, False, False, True, 't1w_mask', 't2w_aseg'),
        (False, False, False, False, None, None),
    ],
)
def test_derivatives(
    tmp_path: Path,
    t1w_mask: bool,
    t1w_aseg: bool,
    t2w_mask: bool,
    t2w_aseg: bool,
    mask: str | None,
    aseg: str | None,
):
    bids_dir = tmp_path / 'bids'
    _create_bids_dir(bids_dir)
    deriv_dir = tmp_path / 'derivatives'
    _create_bids_derivs(
        deriv_dir, t1w_mask=t1w_mask, t1w_aseg=t1w_aseg, t2w_mask=t2w_mask, t2w_aseg=t2w_aseg
    )

    derivatives = bids.Derivatives(bids_dir)
    assert derivatives.mask is None
    assert derivatives.t1w_mask is None
    assert derivatives.t2w_mask is None
    assert derivatives.aseg is None
    assert derivatives.t1w_aseg is None
    assert derivatives.t2w_aseg is None

    derivatives.populate(deriv_dir, subject_id='01')
    if mask:
        assert derivatives.mask == getattr(derivatives, mask)
        assert derivatives.references[mask]
    else:
        assert derivatives.mask is None
    if aseg:
        assert derivatives.aseg == getattr(derivatives, aseg)
        assert derivatives.references[aseg]
    else:
        assert derivatives.aseg == None
