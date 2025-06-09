from pathlib import Path

import nibabel as nb
import numpy as np
import pytest
from nibabel.orientations import aff2axcodes

from nibabies.workflows.anatomical.preproc import (
    _normalize_roi,
    init_conform_derivative_wf,
    init_csf_norm_wf,
)

EXPECTED_CSF_NORM = np.array([[[49, 75], [23, 75]], [[77, 80], [33, 3]]], dtype='uint8')


@pytest.fixture
def anat_file(tmp_path):
    data = np.array([[[49, 73], [23, 73]], [[77, 80], [33, 3]]], dtype='uint8')
    img = nb.Nifti1Image(data, np.eye(4))
    out = tmp_path / 'input.nii.gz'
    img.to_filename(out)
    return out


def test_csf_norm_wf(tmp_path, anat_file):
    tpms = []
    for tpm, thresh in (('gm', 25), ('wm', 75), ('csf', 50)):
        name = tmp_path / f'{tpm}.nii.gz'
        anat_img = nb.load(anat_file)
        anat_data = np.asanyarray(nb.load(anat_file).dataobj)

        binmask = anat_data > thresh
        masked = (binmask * 1).astype('uint8')
        mask = nb.Nifti1Image(masked, anat_img.affine)
        mask.to_filename(name)
        tpms.append(name)

    wf = init_csf_norm_wf()
    wf.base_dir = tmp_path
    wf.inputs.inputnode.anat_preproc = anat_file
    wf.inputs.inputnode.anat_tpms = tpms

    # verify workflow runs
    wf.run()

    # verify function works as expected
    outfile = _normalize_roi(anat_file, tpms[2])
    assert np.array_equal(
        np.asanyarray(nb.load(outfile).dataobj),
        EXPECTED_CSF_NORM,
    )
    Path(outfile).unlink()


@pytest.mark.parametrize(
    ('affine_mismatch', 'ornt_mismatch'),
    [
        (False, False),
        (True, False),
        (False, True),
        (True, True),
    ],
)
def test_conform_derivative_wf(tmp_path, anat_file, affine_mismatch, ornt_mismatch):
    deriv = tmp_path / 'mask.nii.gz'
    ref_img = nb.load(anat_file)
    aff = ref_img.affine.copy()
    if affine_mismatch:
        # Alter affine slightly
        aff[:3, :3] += 0.01
        assert not np.array_equal(aff, ref_img.affine)

    img = ref_img.__class__(ref_img.dataobj, affine=aff)
    if ornt_mismatch:
        from niworkflows.interfaces.nibabel import reorient_image

        img = reorient_image(img, target_ornt='LPI')
        assert aff2axcodes(img.affine) != aff2axcodes(ref_img.affine)

    img.to_filename(deriv)
    wf = init_conform_derivative_wf(in_file=deriv)
    wf.base_dir = tmp_path
    wf.inputs.inputnode.ref_file = anat_file

    wf.run()

    output = list((tmp_path / 'conform_derivative_wf' / 'match_header').glob('*.nii.gz'))
    assert output
    out_file = output[0]
    out_img = nb.load(out_file)
    assert np.array_equal(out_img.affine, ref_img.affine)
    assert aff2axcodes(out_img.affine) == aff2axcodes(ref_img.affine)
