# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright 2023 The NiPreps Developers <nipreps@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# We support and encourage derived works from this project, please read
# about our expectations at
#
#     https://www.nipreps.org/community/licensing/
#
"""Check the configuration module and file."""

import os
from unittest.mock import patch

import pytest
from niworkflows.utils.spaces import format_reference
from toml import loads

from nibabies import config
from nibabies.data import load as load_data


def _reset_config():
    """
    Forcibly reload the configuration module to restore defaults.
    .. caution::
      `importlib.reload` creates new sets of objects, but will not remove
      previous references to those objects."""
    import importlib

    importlib.reload(config)


def test_reset_config():
    execution = config.execution
    execution.bids_dir = 'TESTING'
    assert config.execution.bids_dir == 'TESTING'
    _reset_config()
    assert config.execution.bids_dir is None
    # Even though the config module was reset,
    # previous references to config classes
    # have not been touched.
    assert execution.bids_dir == 'TESTING'


def test_config_spaces():
    """Check that all necessary spaces are recorded in the config."""
    settings = loads(load_data.readable('tests/config.toml').read_text())
    for sectionname, configs in settings.items():
        if sectionname != 'environment':
            section = getattr(config, sectionname)
            section.load(configs, init=False)
    config.nipype.init()
    config.loggers.init()
    age = 8
    spaces = _load_spaces(age)
    assert 'MNI152NLin6Asym:res-2' not in [str(s) for s in spaces.get_standard(full_spec=True)]

    assert 'MNI152NLin6Asym_res-2' not in [
        format_reference((s.fullname, s.spec))
        for s in spaces.references
        if s.standard and s.dim == 3
    ]
    # Only enabled if CIFTI is requested
    assert 'MNI152NLin6Asym:res-2' not in [str(s) for s in spaces.get_standard(full_spec=True)]
    assert 'MNI152NLin6Asym_res-2' not in [
        format_reference((s.fullname, s.spec))
        for s in spaces.references
        if s.standard and s.dim == 3
    ]

    config.execution.output_spaces = None

    with pytest.raises(RuntimeError):
        spaces = _load_spaces(None)

    config.execution.output_spaces = None
    config.workflow.cifti_output = '91k'
    spaces = _load_spaces(1)

    assert [str(s) for s in spaces.get_standard(full_spec=True)] == [
        'MNIInfant:cohort-1:res-native',  # Default output space
        'MNI152NLin6Asym:res-2',
        'MNIInfant:cohort-1:res-2',  # CIFTI: MNIInfant (2x2x2) -> MNI152NLin6Asym (2x2x2)
    ]

    assert [
        format_reference((s.fullname, s.spec))
        for s in spaces.references
        if s.standard and s.dim == 3
    ] == ['MNIInfant_cohort-1_res-native', 'MNI152NLin6Asym_res-2', 'MNIInfant_cohort-1_res-2']
    _reset_config()

    config.execution.output_spaces = None
    config.workflow.cifti_output = '170k'
    spaces = _load_spaces(1)

    assert [str(s) for s in spaces.get_standard(full_spec=True)] == [
        'MNIInfant:cohort-1:res-native',  # Default output space
        'MNI152NLin6Asym:res-1',
        'MNIInfant:cohort-1:res-1',
    ]

    assert [
        format_reference((s.fullname, s.spec))
        for s in spaces.references
        if s.standard and s.dim == 3
    ] == ['MNIInfant_cohort-1_res-native', 'MNI152NLin6Asym_res-1', 'MNIInfant_cohort-1_res-1']
    _reset_config()


@pytest.mark.parametrize(
    ('master_seed', 'ants_seed', 'numpy_seed'), [(1, 17612, 8272), (100, 19094, 60232)]
)
def test_prng_seed(master_seed, ants_seed, numpy_seed):
    """Ensure seeds are properly tracked"""
    seeds = config.seeds
    with patch.dict(os.environ, {}):
        seeds.load({'_random_seed': master_seed}, init=True)
        assert seeds.master == master_seed
        assert seeds.ants == ants_seed
        assert seeds.numpy == numpy_seed
        assert os.getenv('ANTS_RANDOM_SEED') == str(ants_seed)

    _reset_config()
    for seed in ('_random_seed', 'master', 'ants', 'numpy'):
        assert getattr(config.seeds, seed) is None


def _load_spaces(age):
    from nibabies.workflows.base import init_execution_spaces, init_workflow_spaces

    # Conditional based on workflow necessities
    spaces = init_workflow_spaces(init_execution_spaces(), age)
    return spaces


def test_hash_config():
    # This may change with changes to config defaults / new attributes!
    expected = 'cfee5aaf'
    assert config.hash_config(config.get()) == expected
    _reset_config()

    config.execution.log_level = 5  # non-vital attributes do not matter
    assert config.hash_config(config.get()) == expected
    _reset_config()

    # but altering a vital attribute will create a new hash
    config.workflow.surface_recon_method = 'mcribs'
    assert config.hash_config(config.get()) != expected
    _reset_config()
