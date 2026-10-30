"""py.test configuration"""

import json
from pathlib import Path

import nibabel as nb
import numpy as np
import pytest

from nibabies.data import load as load_data


@pytest.fixture(autouse=True)
def _populate_namespace(doctest_namespace, tmp_path):
    doctest_namespace['testdir'] = tmp_path
    doctest_namespace['datadir'] = load_data()


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
