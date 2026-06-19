from pathlib import Path
from uuid import uuid4

import nibabel as nb
import numpy as np
import pytest

from ..nibabel import CropToROI, MergeLabelROIs, ReorientImage


def _make_vol(shape, fill, affine, tmp_path, name):
    data = np.zeros(shape, dtype='int16')
    for idx, value in fill:
        data[idx] = value
    out = str(tmp_path / name)
    nb.Nifti1Image(data, affine).to_filename(out)
    return out


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


def test_crop_to_roi(tmp_path):
    affine = np.diag([2.0, 2.0, 2.0, 1.0])
    affine[:3, 3] = [-10.0, -12.0, -14.0]
    fill = [((i, j, k), 1) for i in (9, 10) for j in (9, 10) for k in (9, 10)]
    roi = _make_vol((20, 20, 20), fill, affine, tmp_path, 'roi.nii.gz')

    res = CropToROI(in_file=roi, margin_mm=2.0).run(cwd=str(tmp_path))

    img = nb.load(res.outputs.out_file)
    data = np.asanyarray(img.dataobj)
    # 2-voxel bbox + 1-voxel (2mm) margin each side; all labeled voxels retained
    assert img.shape == (4, 4, 4)
    assert int(data.sum()) == 8
    # world coordinates of the labeled voxels are unchanged by cropping
    orig = np.asanyarray(nb.load(roi).dataobj)
    cropped_world = sorted(
        nb.affines.apply_affine(img.affine, v).tolist() for v in np.argwhere(data > 0)
    )
    orig_world = sorted(nb.affines.apply_affine(affine, v).tolist() for v in np.argwhere(orig > 0))
    assert np.allclose(cropped_world, orig_world)


def test_merge_label_rois(tmp_path):
    affine = np.diag([2.0, 2.0, 2.0, 1.0])
    template = _make_vol((20, 20, 20), [], affine, tmp_path, 'template.nii.gz')

    aff_a = affine.copy()
    aff_a[:3, 3] = nb.affines.apply_affine(affine, [4, 4, 4])
    a = _make_vol((3, 3, 3), [((0, 0, 0), 26)], aff_a, tmp_path, 'a.nii.gz')
    aff_b = affine.copy()
    aff_b[:3, 3] = nb.affines.apply_affine(affine, [10, 10, 10])
    b = _make_vol((3, 3, 3), [((1, 1, 1), 58)], aff_b, tmp_path, 'b.nii.gz')

    res = MergeLabelROIs(in_files=[a, b], template=template).run(cwd=str(tmp_path))

    img = nb.load(res.outputs.out_file)
    assert img.shape == (20, 20, 20)
    assert np.allclose(img.affine, affine)
    data = np.asanyarray(img.dataobj)
    assert data[4, 4, 4] == 26  # a's voxel
    assert data[11, 11, 11] == 58  # b's voxel


@pytest.mark.parametrize(
    ('overlap_method', 'label'),
    [('first', 11), ('last', 77), ('fail', None)],
)
def test_merge_label_rois_overlap(tmp_path, overlap_method, label):
    affine = np.diag([2.0, 2.0, 2.0, 1.0])
    template = _make_vol((20, 20, 20), [], affine, tmp_path, 'template.nii.gz')
    aff = affine.copy()
    aff[:3, 3] = nb.affines.apply_affine(affine, [5, 5, 5])
    # three ROIs all claiming voxel [5, 5, 5], in order
    first = _make_vol((2, 2, 2), [((0, 0, 0), 11)], aff, tmp_path, 'first.nii.gz')
    second = _make_vol((2, 2, 2), [((0, 0, 0), 99)], aff, tmp_path, 'second.nii.gz')
    third = _make_vol((2, 2, 2), [((0, 0, 0), 77)], aff, tmp_path, 'third.nii.gz')
    merge = MergeLabelROIs(
        in_files=[first, second, third], template=template, overlap_method=overlap_method
    )

    if label is None:
        with pytest.raises(ValueError, match='[Oo]verlap'):
            merge.run(cwd=str(tmp_path))
        return

    res = merge.run(cwd=str(tmp_path))
    seg = np.asanyarray(nb.load(res.outputs.out_file).dataobj)
    assert seg[5, 5, 5] == label
    overlap = np.asanyarray(nb.load(res.outputs.overlap_file).dataobj)
    assert overlap[5, 5, 5] == 3
