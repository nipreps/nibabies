
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utilities to handle BIDS inputs."""
import os
import sys
import json
from pathlib import Path


def write_bidsignore(deriv_dir):
    # TODO: Port to niworkflows
    bids_ignore = (
        "*.html", "logs/", "figures/",  # Reports
        "*_xfm.*",  # Unspecified transform files
        "*.surf.gii",  # Unspecified structural outputs
        # Unspecified functional outputs
        "*_boldref.nii.gz", "*_bold.func.gii",
        "*_mixing.tsv", "*_AROMAnoiseICs.csv", "*_timeseries.tsv",
    )
    ignore_file = Path(deriv_dir) / ".bidsignore"

    ignore_file.write_text("\n".join(bids_ignore) + "\n")


def write_derivative_description(bids_dir, deriv_dir):
    from ..__about__ import (
        __version__,
        __packagename__,
        DOWNLOAD_URL,
    )

    bids_dir = Path(bids_dir)
    deriv_dir = Path(deriv_dir)
    desc = {
        'Name': 'NiBabies: Neuroimaging preprocessing workflows for babies',
        'BIDSVersion': '1.4.0',
        'DatasetType': 'derivative',
        'GeneratedBy': [{
            'Name': __packagename__,
            'Version': __version__,
            'CodeURL': DOWNLOAD_URL,
        }],
        'HowToAcknowledge': 'TODO',
    }

    # Keys that can only be set by environment
    if 'NIBABIES_DOCKER_TAG' in os.environ:
        desc['GeneratedBy'][0]['Container'] = {
            "Type": "docker",
            "Tag": f"nipreps/fmriprep:{os.environ['NIBABIES_DOCKER_TAG']}"
        }
    if 'NIBABIES_SINGULARITY_URL' in os.environ:
        desc['GeneratedBy'][0]['Container'] = {
            "Type": "singularity",
            "URI": os.getenv('NIBABIES_SINGULARITY_URL')
        }

    # Keys deriving from source dataset
    orig_desc = {}
    fname = bids_dir / 'dataset_description.json'
    if fname.exists():
        orig_desc = json.loads(fname.read_text())

    if 'DatasetDOI' in orig_desc:
        desc['SourceDatasets'] = [{
            'URL': f'https://doi.org/{orig_desc["DatasetDOI"]}',
            'DOI': orig_desc['DatasetDOI']
        }]
    if 'License' in orig_desc:
        desc['License'] = orig_desc['License']

    Path.write_text(deriv_dir / 'dataset_description.json', json.dumps(desc, indent=4))
