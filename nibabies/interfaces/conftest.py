"""py.test configuration"""
from pathlib import Path
import pytest
from tempfile import TemporaryDirectory


FILES = (
    'functional.nii',
    'anatomical.nii',
    'epi.nii',
    'T1w.nii',
    'rois.nii',
    'func_to_struct.mat',
    'labels.nii',
    'sub-01_run-01_echo-1_bold.nii.gz',
    'sub-01_run-01_echo-2_bold.nii.gz',
    'sub-01_run-01_echo-3_bold.nii.gz',
)


@pytest.fixture(scope="package")
def data_dir():
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        for fname in FILES:
            Path.touch(tmp_path / fname)
        yield tmp_path


@pytest.fixture(autouse=True)
def set_namespace(doctest_namespace, data_dir):
    doctest_namespace["data_dir"] = data_dir
