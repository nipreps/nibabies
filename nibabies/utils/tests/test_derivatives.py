from nibabies.utils.derivatives import copy_derivatives


def test_copy_derivatives(tmp_path):
    precomp = tmp_path / 'precomputed'
    precomp.mkdir()
    out = tmp_path / 'out'
    out.mkdir()

    mask = precomp / 'mask.nii.gz'
    mask.touch()
    aseg = precomp / 'aseg.nii.gz'
    aseg.touch()
    aseg_meta = precomp / 'aseg.json'
    aseg_meta.touch()

    derivs = {
        't2w_mask': str(mask),
        't2w_aseg': str(aseg),
        'transforms': {},
    }

    copy_derivatives(derivs, out, 'anat', 'sub-01')
    outpath = out / 'sub-01' / 'anat'

    for fl in (mask, aseg, aseg_meta):
        assert (outpath / fl.name).exists()

    copy_derivatives(derivs, out, 'anat', 'sub-01', 'ses-a')
    outpath = out / 'sub-01' / 'ses-a' / 'anat'

    for fl in (mask, aseg, aseg_meta):
        assert (outpath / fl.name).exists()
