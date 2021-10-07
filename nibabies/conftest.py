"""py.test configuration"""
from pathlib import Path
from tempfile import TemporaryDirectory
from pkg_resources import resource_filename

import pytest

FILES = (
    "functional.nii",
    "anatomical.nii",
    "func.dlabel.nii",
    "func.dtseries.nii",
    "epi.nii",
    "T1w.nii",
    "func_to_struct.mat",
    "atlas.nii",
    "label_list.txt",
    "sub-01_run-01_echo-1_bold.nii.gz",
    "sub-01_run-01_echo-2_bold.nii.gz",
    "sub-01_run-01_echo-3_bold.nii.gz",
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
    doctest_namespace["test_data"] = Path(resource_filename("nibabies", "tests/data"))
