"""py.test configuration"""

import json
from pathlib import Path
from tempfile import TemporaryDirectory

import nibabel as nb
import numpy as np
import pytest

from nibabies.data import load as load_data

FILES = (
    'functional.nii',
    'anatomical.nii',
    'func.dlabel.nii',
    'func.dtseries.nii',
    'epi.nii',
    'T1w.nii',
    'func_to_struct.mat',
    'atlas.nii',
    'label_list.txt',
    'sub-01_run-01_echo-1_bold.nii.gz',
    'sub-01_run-01_echo-2_bold.nii.gz',
    'sub-01_run-01_echo-3_bold.nii.gz',
)


@pytest.fixture(scope='package')
def data_dir():
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        for fname in FILES:
            Path.touch(tmp_path / fname)
        yield tmp_path


@pytest.fixture(autouse=True)
def _populate_namespace(doctest_namespace, data_dir):
    doctest_namespace['data_dir'] = data_dir
    doctest_namespace['test_data'] = load_data.cached('../tests/data')
    doctest_namespace['Path'] = Path


@pytest.fixture
def minimal_bids(tmp_path):
    bids = tmp_path / 'bids'
    bids.mkdir()
    Path.write_text(
        bids / 'dataset_description.json', json.dumps({'Name': 'Test DS', 'BIDSVersion': '1.8.0'})
    )
    T1w = bids / 'sub-01' / 'anat' / 'sub-01_T1w.nii.gz'
    T1w.parent.mkdir(parents=True)
    nb.Nifti1Image(np.zeros((5, 5, 5)), np.eye(4)).to_filename(T1w)
    return bids
