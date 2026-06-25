from pathlib import Path
from uuid import uuid4

import nibabel as nb
import numpy as np
import pytest

from ..nibabel import ReorientImage


def create_save_img(ornt: str):
    data = np.random.rand(2, 2, 2)
    img = nb.Nifti1Image(data, affine=np.eye(4))
    # img will always be in RAS at the start
    ras = nb.orientations.axcodes2ornt('RAS')
    if ornt != 'RAS':
        new = nb.orientations.axcodes2ornt(ornt)
        xfm = nb.orientations.ornt_transform(ras, new)
        img = img.as_reoriented(xfm)
    out_file = f'{uuid4()}.nii.gz'
    img.to_filename(out_file)
    return out_file


@pytest.mark.parametrize(
    ('in_ornt', 'out_ornt'),
    [
        ('RAS', 'RAS'),
        ('RAS', 'LAS'),
        ('LAS', 'RAS'),
        ('RAS', 'RPI'),
        ('LPI', 'RAS'),
    ],
)
def test_reorient_image(tmpdir, in_ornt, out_ornt):
    tmpdir.chdir()

    in_file = create_save_img(ornt=in_ornt)
    in_img = nb.load(in_file)
    assert ''.join(nb.aff2axcodes(in_img.affine)) == in_ornt

    # test string representation
    res = ReorientImage(in_file=in_file, target_orientation=out_ornt).run()
    out_file = res.outputs.out_file
    out_img = nb.load(out_file)
    assert ''.join(nb.aff2axcodes(out_img.affine)) == out_ornt
    Path(out_file).unlink()

    # test with target file
    target_file = create_save_img(ornt=out_ornt)
    target_img = nb.load(target_file)
    assert ''.join(nb.aff2axcodes(target_img.affine)) == out_ornt
    res = ReorientImage(in_file=in_file, target_file=target_file).run()
    out_file = res.outputs.out_file
    out_img = nb.load(out_file)
    assert ''.join(nb.aff2axcodes(out_img.affine)) == out_ornt

    # cleanup
    for f in (in_file, target_file, out_file):
        Path(f).unlink()
