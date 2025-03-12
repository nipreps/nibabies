import pytest

from ..hbcd import main


# create a FreeSurfer / MCRIBS directory structure
@pytest.fixture
def freesurfer(tmp_path):
    root = tmp_path / 'freesurfer'
    for s in ('fsaverage', 'sub-01_ses-1', 'sub-01_ses-2'):
        (root / s / 'mri').mkdir(parents=True, exist_ok=True)
        (root / s / 'mri' / 'T1.mgz').touch()
        (root / s / 'surf').mkdir(parents=True, exist_ok=True)
        (root / s / 'surf' / 'lh.pial').touch()
    return root


@pytest.fixture
def mcribs(tmp_path):
    root = tmp_path / 'mcribs'
    for s in ('sub-01_ses-1', 'sub-01_ses-2'):
        (root / s / 'TissueSeg').mkdir(parents=True, exist_ok=True)
        orig = root / s / 'TissueSeg' / f'{s}_all_labels.nii.gz'
        orig.touch()
        (root / s / 'TissueSeg' / f'{s}_all_labels_manedit.nii.gz').symlink_to(orig)
    return root


def test_hbcd_restructure(freesurfer, mcribs):
    # running without options is fine
    main([])

    main(['--fs', str(freesurfer), '--mcribs', str(mcribs)])
    assert sorted(x.name for x in freesurfer.iterdir()) == ['fsaverage', 'sub-01']
    assert sorted(x.name for x in (freesurfer / 'sub-01').iterdir()) == ['ses-1', 'ses-2']

    assert [x.name for x in mcribs.iterdir()] == ['sub-01']
    linkd = mcribs / 'sub-01' / 'ses-1' / 'TissueSeg' / 'sub-01_ses-1_all_labels_manedit.nii.gz'
    assert linkd.exists()
    assert not linkd.is_symlink()

    # and run again should not fail
    main(['--fs', str(freesurfer), '--mcribs', str(mcribs)])
