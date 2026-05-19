"""Tests for the NiBabies fmriprep extension descriptor."""

import argparse
from unittest.mock import MagicMock, patch

import pytest


def _make_parser():
    """Return a bare ArgumentParser to extend."""
    return argparse.ArgumentParser()


def test_extend_parser_adds_nibabies_group():
    """extend_parser adds an argument group named 'NiBabies options'."""
    from nibabies.cli import extend_parser

    parser = _make_parser()
    extend_parser(parser)
    group_titles = [g.title for g in parser._action_groups]
    assert 'NiBabies options' in group_titles


def test_extend_parser_age_months_flag():
    """--age-months is present and parses to int."""
    from nibabies.cli import extend_parser

    parser = _make_parser()
    extend_parser(parser)
    opts = parser.parse_args(['--age-months', '6'])
    assert opts.age_months == 6


def test_extend_parser_surface_recon_method_none():
    """--surface-recon-method none/None converts to Python None."""
    from nibabies.cli import extend_parser

    parser = _make_parser()
    extend_parser(parser)
    for val in ('none', 'None'):
        opts = parser.parse_args(['--surface-recon-method', val])
        assert opts.surface_recon_method is None


def test_extend_parser_hmc_bold_frame_auto():
    """--hmc-bold-frame accepts 'auto'."""
    from nibabies.cli import extend_parser

    parser = _make_parser()
    extend_parser(parser)
    opts = parser.parse_args(['--hmc-bold-frame', 'auto'])
    assert opts.hmc_bold_frame == 'auto'


def test_extend_parser_no_multi_step_reg():
    """--no-multi-step-reg sets multi_step_reg to False."""
    from nibabies.cli import extend_parser

    parser = _make_parser()
    extend_parser(parser)
    opts = parser.parse_args(['--no-multi-step-reg'])
    assert opts.multi_step_reg is False


def test_populate_config_renames_recon_method():
    """surface_recon_method on opts becomes recon_method in the namespace."""
    from nibabies.cli import populate_config

    ext = MagicMock()
    opts = argparse.Namespace(
        age_months=6,
        segmentation_atlases_dir=None,
        fd_radius=45,
        surface_recon_method='infantfs',
        reference_anat='T2w',
        hmc_bold_frame=16,
        norm_csf=False,
        multi_step_reg=True,
    )
    populate_config(opts, ext)
    ext.set.assert_any_call('recon_method', 'infantfs')
    ext.set.assert_any_call('reference_anat', 'T2w')


def test_populate_config_writes_all_fields():
    """populate_config writes all 8 nibabies fields to the extension namespace."""
    from nibabies.cli import populate_config

    ext = MagicMock()
    opts = argparse.Namespace(
        age_months=3,
        segmentation_atlases_dir=None,
        fd_radius=45,
        surface_recon_method='auto',
        reference_anat=None,
        hmc_bold_frame=16,
        norm_csf=False,
        multi_step_reg=True,
    )
    populate_config(opts, ext)
    expected_keys = {
        'age_months',
        'segmentation_atlases',
        'fd_radius',
        'recon_method',
        'reference_anat',
        'hmc_bold_frame',
        'norm_csf',
        'multi_step_reg',
    }
    actual_keys = {call.args[0] for call in ext.set.call_args_list}
    assert actual_keys == expected_keys


def test_descriptor_required_attrs():
    """NibabiesExtension declares name, version, fmriprep_compat, contracts."""
    from nibabies.descriptor import NibabiesExtension

    desc = NibabiesExtension()
    assert desc.name == 'nibabies'
    assert desc.version
    assert desc.fmriprep_compat
    assert 'anat_fit' in desc.contracts


def test_descriptor_init_config_raises_without_age():
    """init_config raises ExtensionConfigError when age_months is not set."""
    from fmriprep.extensions.exceptions import ExtensionConfigError

    from nibabies.descriptor import NibabiesExtension

    desc = NibabiesExtension()
    with patch.object(desc, 'get', return_value=None), patch.object(desc, 'set'):
        with pytest.raises(ExtensionConfigError):
            desc.init_config()


def test_descriptor_init_anat_fit_wf_routes_single_anat():
    """init_anat_fit_wf calls single-anat builder when only t1w OR t2w present."""
    from nibabies.descriptor import NibabiesExtension

    desc = NibabiesExtension()

    def _fake_get(key, default=None):
        return {
            'age_months': 6,
            'recon_method': 'auto',
            'segmentation_atlases': None,
            'reference_anat': None,
        }.get(key, default)

    mock_config = MagicMock()
    mock_config.workflow.cifti_output = False
    with (
        patch.object(desc, 'get', side_effect=_fake_get),
        patch('nibabies.descriptor.init_infant_single_anat_fit_wf') as mock_single,
        patch('nibabies.descriptor.init_infant_anat_fit_wf') as mock_dual,
        patch('nibabies.descriptor.config', mock_config),
    ):
        desc.init_anat_fit_wf(
            t1w=['t1.nii.gz'],
            t2w=[],
            precomputed={},
            flair=[],
            bids_root='/data',
            longitudinal=False,
            msm_sulc=False,
            omp_nthreads=1,
            output_dir='/out',
            skull_strip_fixed_seed=False,
            skull_strip_mode='force',
            skull_strip_template=MagicMock(),
            sloppy=False,
            spaces=MagicMock(),
        )
        mock_single.assert_called_once()
        mock_dual.assert_not_called()


def test_descriptor_init_anat_fit_wf_routes_dual_anat():
    """init_anat_fit_wf calls dual-anat builder when both t1w and t2w present."""
    from nibabies.descriptor import NibabiesExtension

    desc = NibabiesExtension()

    def _fake_get(key, default=None):
        return {
            'age_months': 6,
            'recon_method': 'auto',
            'segmentation_atlases': None,
            'reference_anat': None,
        }.get(key, default)

    mock_config = MagicMock()
    mock_config.workflow.cifti_output = False
    with (
        patch.object(desc, 'get', side_effect=_fake_get),
        patch('nibabies.descriptor.init_infant_single_anat_fit_wf') as mock_single,
        patch('nibabies.descriptor.init_infant_anat_fit_wf') as mock_dual,
        patch('nibabies.descriptor.config', mock_config),
    ):
        desc.init_anat_fit_wf(
            t1w=['t1.nii.gz'],
            t2w=['t2.nii.gz'],
            precomputed={},
            flair=[],
            bids_root='/data',
            longitudinal=False,
            msm_sulc=False,
            omp_nthreads=1,
            output_dir='/out',
            skull_strip_fixed_seed=False,
            skull_strip_mode='force',
            skull_strip_template=MagicMock(),
            sloppy=False,
            spaces=MagicMock(),
        )
        mock_dual.assert_called_once()
        mock_single.assert_not_called()
