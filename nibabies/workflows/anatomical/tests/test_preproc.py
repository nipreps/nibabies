import typing as ty
from pathlib import Path

import nibabel as nb
import numpy as np
import pytest

from nibabies.workflows.anatomical.preproc import _normalize_roi, init_csf_norm_wf

EXPECTED_CSF_NORM = np.array([[[10, 73], [73, 29]], [[77, 80], [6, 16]]], dtype='uint8')


@pytest.fixture
def csf_norm_data(tmp_path) -> ty.Generator[tuple[Path, list[Path]], None, None]:
    np.random.seed(10)

    in_file = tmp_path / 'input.nii.gz'
    data = np.random.randint(1, 101, size=(2, 2, 2), dtype='uint8')
    img = nb.Nifti1Image(data, np.eye(4))
    img.to_filename(in_file)

    masks = []
    for tpm in ('gm', 'wm', 'csf'):
        name = tmp_path / f'{tpm}.nii.gz'
        binmask = data > np.random.randint(10, 90)
        masked = (binmask * 1).astype('uint8')
        mask = nb.Nifti1Image(masked, img.affine)
        mask.to_filename(name)
        masks.append(name)

    yield in_file, masks

    in_file.unlink()
    for m in masks:
        m.unlink()


def test_csf_norm_wf(tmp_path, csf_norm_data):
    anat, tpms = csf_norm_data
    wf = init_csf_norm_wf()
    wf.base_dir = tmp_path

    wf.inputs.inputnode.anat_preproc = anat
    wf.inputs.inputnode.anat_tpms = tpms

    # verify workflow runs
    wf.run()

    # verify function works as expected
    outfile = _normalize_roi(anat, tpms[2])
    assert np.array_equal(
        np.asanyarray(nb.load(outfile).dataobj),
        EXPECTED_CSF_NORM,
    )
    Path(outfile).unlink()
