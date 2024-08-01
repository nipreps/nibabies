from niworkflows.utils.testing import generate_bids_skeleton

from nibabies.utils.derivatives import (
    collect_anatomical_derivatives,
    collect_functional_derivatives,
    copy_derivatives,
)

from . import DERIV_SKELETON


def test_collect_derivatives(tmp_path):
    deriv_dir = tmp_path / 'derivatives'
    generate_bids_skeleton(deriv_dir, str(DERIV_SKELETON))
    output_spaces = ['MNIInfant:cohort-1']

    anat_cache = collect_anatomical_derivatives(
        derivatives_dir=deriv_dir,
        subject_id='01',
        session_id=None,
        std_spaces=output_spaces,
    )
    for suffix in ('preproc', 'mask', 'dseg'):
        assert anat_cache[f't2w_{suffix}']
    assert len(anat_cache['t2w_tpms']) == 3
    xfms = anat_cache['transforms']
    for space in output_spaces:
        assert xfms[space]['reverse']
        assert xfms[space]['forward']
    for surface in (
        'white',
        'pial',
        'midthickness',
        'sphere',
        'thickness',
        'sulc',
        'sphere_reg',
        'sphere_reg_fsLR',
    ):
        assert len(anat_cache[surface]) == 2

    func_cache = collect_functional_derivatives(deriv_dir, {'subject': '01'}, None)
    for val in ('hmc_boldref', 'coreg_boldref', 'hmc'):
        assert func_cache[val]


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
