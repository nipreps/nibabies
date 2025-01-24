import shutil
from pathlib import Path

import pytest

from nibabies.interfaces.mcribs import MCRIBReconAll

SUBJECT_ID = 'X'


@pytest.fixture
def mcribs_directory(tmp_path):
    def make_tree(path, tree):
        for d, fls in tree.items():
            (path / d).mkdir(exist_ok=True)
            for f in fls:
                (path / d / f).touch()

    root = tmp_path / 'mcribs'
    surfrecon = root / SUBJECT_ID / 'SurfReconDeformable' / SUBJECT_ID
    surfrecon.mkdir(parents=True, exist_ok=True)
    make_tree(surfrecon, MCRIBReconAll._expected_files['surfrecon'])
    autorecon = root / SUBJECT_ID / 'freesurfer' / SUBJECT_ID
    autorecon.mkdir(parents=True, exist_ok=True)
    make_tree(autorecon, MCRIBReconAll._expected_files['autorecon'])

    yield root

    shutil.rmtree(root)


def test_MCRIBReconAll(mcribs_directory):
    t2w = Path('T2w.nii.gz')
    t2w.touch()

    surfrecon = MCRIBReconAll(
        subject_id=SUBJECT_ID,
        surfrecon=True,
        surfrecon_method='Deformable',
        join_thresh=1.0,
        fast_collision=True,
    )

    # Requires T2w input
    with pytest.raises(AttributeError):
        surfrecon.cmdline  # noqa

    surfrecon.inputs.t2w_file = t2w
    # Since no existing directory is found, will run fresh
    assert 'MCRIBReconAll --deformablefastcollision --deformablejointhresh' in surfrecon.cmdline

    # But should not need to run again
    surfrecon.inputs.outdir = mcribs_directory
    assert surfrecon.cmdline == 'echo MCRIBReconAll: nothing to do'

    t2w.unlink()

    autorecon = MCRIBReconAll(
        subject_id=SUBJECT_ID,
        autorecon_after_surf=True,
    )
    # No need for T2w here
    assert autorecon.cmdline == 'MCRIBReconAll --autoreconaftersurf X'
    autorecon.inputs.outdir = mcribs_directory
    assert autorecon.cmdline == 'echo MCRIBReconAll: nothing to do'
