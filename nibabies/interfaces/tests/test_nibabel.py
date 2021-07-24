import uuid

import nibabel as nb
import numpy as np
import pytest

from ..nibabel import MergeROIs

@pytest.fixture
def create_roi(tmp_path):
    files = []
    def _create_roi(affine, img_data, roi_index):
        img_data[tuple(roi_index)] = 1
        nii = nb.Nifti1Image(img_data, affine)
        filename = tmp_path / f"{str(uuid.uuid4())}.nii.gz"
        files.append(filename)
        nii.to_filename(filename)
        return filename
    yield _create_roi

    for f in files:
        f.unlink()

# create a slightly off affine
bad_affine = np.eye(4)
bad_affine[0, -1] = -1

@pytest.mark.parametrize("affine,data,roi_index,error,err_message", [
    (np.eye(4), np.zeros((2,2,2,2), dtype=int), [1, 0], None, None),
    (np.eye(4), np.zeros((2,2,3,2), dtype=int), [1, 0], True, "Mismatch in image shape"),
    (bad_affine, np.zeros((2,2,2,2), dtype=int), [1, 0], True, "Mismatch in affine"),
    (np.eye(4), np.zeros((2,2,2,2), dtype=int), [0, 0, 0], True, "Overlapping ROIs"),
])
def test_merge_rois(tmpdir, create_roi, affine, data, roi_index, error, err_message):
    tmpdir.chdir()
    roi0 = create_roi(np.eye(4), np.zeros((2,2,2,2), dtype=int), [0, 0])
    roi1 = create_roi(np.eye(4), np.zeros((2,2,2,2), dtype=int), [0, 1])
    test_roi = create_roi(affine, data, roi_index)

    merge = MergeROIs(in_files=[roi0, roi1, test_roi])
    if error is None:
        merge.run()
        return
    # otherwise check expected exceptions
    with pytest.raises(AssertionError) as err:
        merge.run()
    assert err_message in str(err.value)
