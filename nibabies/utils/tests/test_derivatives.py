from niworkflows.utils.testing import generate_bids_skeleton

from nibabies.utils.derivatives import (
    collect_anat_derivatives,
    collect_func_derivatives,
    copy_derivatives,
)

from . import DERIV_SKELETON


def test_collect_derivatives(tmp_path):
    deriv_dir = tmp_path / 'derivatives'
    generate_bids_skeleton(deriv_dir, str(DERIV_SKELETON))
    output_spaces = ['MNIInfant:cohort-1']

    anat_cache = collect_anat_derivatives(
        [deriv_dir],
        subject_id='01',
        session_id=None,
        std_spaces=output_spaces,
    )
    images = anat_cache['images']
    for suffix in ('preproc', 'mask', 'dseg'):
        assert images[f't2w_{suffix}']
    assert [p.split('label-')[1].split('_')[0] for p in images['t2w_tpms']] == ['GM', 'WM', 'CSF']
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
        assert len(anat_cache['surfaces'][surface]) == 2

    func_cache = collect_func_derivatives([deriv_dir], {'subject': '01'})
    for val in ('hmc', 'run'):
        assert func_cache['boldrefs'][val]
    for xfm in ('hmc', 'run2anat'):
        assert func_cache['transforms'][xfm]


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


def test_collect_anat_derivatives_empty_later_dataset(tmp_path):
    """An empty later dataset does not clobber an earlier dataset's transforms."""
    deriv_dir = tmp_path / 'derivatives'
    generate_bids_skeleton(deriv_dir, str(DERIV_SKELETON))
    empty = tmp_path / 'empty'
    generate_bids_skeleton(
        empty,
        {'dataset_description': {'Name': 'empty', 'DatasetType': 'derivative'}, '01': []},
    )
    kwargs = {'subject_id': '01', 'std_spaces': ['MNIInfant:cohort-1'], 'session_id': None}

    anat_cache = collect_anat_derivatives([deriv_dir, empty], **kwargs)
    assert anat_cache == collect_anat_derivatives([deriv_dir], **kwargs)
    assert anat_cache['transforms']['MNIInfant:cohort-1']['forward']


def _anat_skeleton(root, files):
    generate_bids_skeleton(
        root,
        {
            'dataset_description': {'Name': 'NiBabies', 'DatasetType': 'derivative'},
            '01': [{'session': None, 'anat': files}],
        },
    )
    return root


def test_collect_anat_derivatives_sphere_reg_current_name(tmp_path):
    """``sphere_reg`` is found under its current ``space-fsaverage`` name."""
    sphere_reg = [
        {
            'hemi': hemi,
            'space': 'fsaverage',
            'desc': 'reg',
            'suffix': 'sphere',
            'extension': '.surf.gii',
        }
        for hemi in ('L', 'R')
    ]
    root = _anat_skeleton(tmp_path / 'deriv', sphere_reg)

    surfaces = collect_anat_derivatives([root], '01', [], None)['surfaces']
    assert [p.split('/')[-1] for p in surfaces['sphere_reg']] == [
        'sub-01_hemi-L_space-fsaverage_desc-reg_sphere.surf.gii',
        'sub-01_hemi-R_space-fsaverage_desc-reg_sphere.surf.gii',
    ]


def test_collect_anat_derivatives_coreg(tmp_path):
    """T1w<->T2w transforms are collected under ``coreg``."""
    xfms = [
        {'from': src, 'to': dst, 'mode': 'image', 'suffix': 'xfm', 'extension': '.txt'}
        for src, dst in (('T1w', 'T2w'), ('T2w', 'T1w'))
    ]
    root = _anat_skeleton(tmp_path / 'deriv', xfms)

    cache = collect_anat_derivatives([root], '01', [], None)
    assert cache['coreg']['t1w2t2w_xfm'].endswith('_from-T1w_to-T2w_mode-image_xfm.txt')
    assert cache['coreg']['t2w2t1w_xfm'].endswith('_from-T2w_to-T1w_mode-image_xfm.txt')
    assert cache['transforms'] == {}
