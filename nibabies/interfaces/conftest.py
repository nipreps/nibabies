from pathlib import Path
from shutil import copytree
from tempfile import TemporaryDirectory

import pytest

try:
    from contextlib import chdir as _chdir
except ImportError:  # PY310
    import os
    from contextlib import contextmanager

    @contextmanager  # type: ignore
    def _chdir(path):
        cwd = os.getcwd()
        os.chdir(path)
        try:
            yield
        finally:
            os.chdir(cwd)


DATA_FILES = (
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
    'xfm0.h5',
    'xfm1.h5',
)


@pytest.fixture(scope='package')
def data_dir():
    with TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        for fname in DATA_FILES:
            Path.touch(tmp_path / fname)
        yield tmp_path


@pytest.fixture(autouse=True)
def _docdir(data_dir, request, tmp_path):
    # Trigger ONLY for the doctests.
    doctest_plugin = request.config.pluginmanager.getplugin('doctest')
    if isinstance(request.node, doctest_plugin.DoctestItem):
        copytree(data_dir, tmp_path, dirs_exist_ok=True)

        # Chdir only for the duration of the test.
        with _chdir(tmp_path):
            yield

    else:
        # For normal tests, we have to yield, since this is a yield-fixture.
        yield
