import json
from pathlib import Path
import uuid

import nibabel as nb
import numpy as np
import pytest

from ..nibabel import MergeROIs, MapLabels


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


def create_image(data, filename):
    nb.Nifti1Image(data, affine=np.eye(4)).to_filename(str(filename))
    return filename


# create a slightly off affine
bad_affine = np.eye(4)
bad_affine[0, -1] = -1


@pytest.mark.parametrize(
    "affine, data, roi_index, error, err_message",
    [
        (np.eye(4), np.zeros((2, 2, 2, 2), dtype=int), [1, 0], None, None),
        (
            np.eye(4),
            np.zeros((2, 2, 3, 2), dtype=int),
            [1, 0],
            True,
            "Mismatch in image shape",
        ),
        (
            bad_affine,
            np.zeros((2, 2, 2, 2), dtype=int),
            [1, 0],
            True,
            "Mismatch in affine",
        ),
        (
            np.eye(4),
            np.zeros((2, 2, 2, 2), dtype=int),
            [0, 0, 0],
            True,
            "Overlapping ROIs",
        ),
    ],
)
def test_merge_rois(tmpdir, create_roi, affine, data, roi_index, error, err_message):
    tmpdir.chdir()
    roi0 = create_roi(np.eye(4), np.zeros((2, 2, 2, 2), dtype=int), [0, 0])
    roi1 = create_roi(np.eye(4), np.zeros((2, 2, 2, 2), dtype=int), [0, 1])
    test_roi = create_roi(affine, data, roi_index)

    merge = MergeROIs(in_files=[roi0, roi1, test_roi])
    if error is None:
        merge.run()
        return
    # otherwise check expected exceptions
    with pytest.raises(AssertionError) as err:
        merge.run()
    assert err_message in str(err.value)


DEFAULT_MAPPING = {5: 1, 6: 1, 7: 2}
DEFAULT_INPUT = np.arange(8).reshape(2, 2, 2)
DEFAULT_OUTPUT = np.asarray([0, 1, 2, 3, 4, 1, 1, 2]).reshape(2, 2, 2)


@pytest.mark.parametrize(
    "data,mapping,tojson,expected",
    [
        (DEFAULT_INPUT, DEFAULT_MAPPING, False, DEFAULT_OUTPUT),
        (DEFAULT_INPUT, DEFAULT_MAPPING, True, DEFAULT_OUTPUT),
    ],
)
def test_map_labels(tmpdir, data, mapping, tojson, expected):
    tmpdir.chdir()
    in_file = create_image(data, Path("test.nii.gz"))
    maplbl = MapLabels(in_file=in_file)
    if tojson:
        map_file = Path('mapping.json')
        map_file.write_text(json.dumps(mapping))
        maplbl.inputs.mappings_file = map_file
    else:
        maplbl.inputs.mappings = mapping
    out_file = maplbl.run().outputs.out_file

    orig = nb.load(in_file).get_fdata()
    labels = nb.load(out_file).get_fdata()
    assert orig.shape == labels.shape
    assert np.all(labels == expected)

    Path(in_file).unlink()
    if tojson:
        Path(map_file).unlink()
