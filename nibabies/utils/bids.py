# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utilities to handle BIDS inputs."""
from __future__ import annotations

import json
import os
import sys
import typing as ty
import warnings
from pathlib import Path

import nibabel as nb
import numpy as np
from bids.layout import BIDSLayout, Query

_spec: dict = {
    't1w_mask': {
        'datatype': 'anat',
        'desc': 'brain',
        'space': 'T1w',
        'suffix': 'mask',
    },
    't1w_aseg': {'datatype': 'anat', 'desc': 'aseg', 'space': 'T1w', 'suffix': 'dseg'},
    't2w_mask': {
        'datatype': 'anat',
        'desc': 'brain',
        'space': 'T2w',
        'suffix': 'mask',
    },
    't2w_aseg': {
        'datatype': 'anat',
        'desc': 'aseg',
        'space': 'T2w',
        'suffix': 'dseg',
    },
}


class Derivatives:
    """A container class for storing precomputed derivatives."""

    def __getattribute__(self, attr):
        """In cases where the spec may change, avoid errors."""
        try:
            return object.__getattribute__(self, attr)
        except AttributeError:
            return None

    def __init__(self, spec: dict | Path | str | None = None, **args):
        self.spec = _spec
        if spec is not None:
            if not isinstance(spec, dict):
                spec: dict = json.loads(Path(spec).read_text())
            self.spec = spec

        self.names = set(self.spec.keys())
        self.references = {name: None for name in self.names}
        for name in self.names:
            setattr(self, name, None)

    def __repr__(self):
        return '\n'.join([name for name in self.names if getattr(self, name)])

    def __contains__(self, val: str):
        return val in self.names

    def __bool__(self):
        return any(getattr(self, name) for name in self.names)

    def populate(
        self, deriv_path, subject_id: str, session_id: str | Query | None = Query.OPTIONAL
    ) -> None:
        """Query a derivatives directory and populate values and references based on the spec."""
        layout = BIDSLayout(deriv_path, validate=False)
        for name, query in self.spec.items():
            items = layout.get(
                subject=subject_id,
                session=session_id,
                extension=['.nii', '.nii.gz'],
                **query,
            )
            if not items or len(items) > 1:
                continue
            item = items[0]

            # Skip if derivative does not have valid metadata
            metadata = item.get_metadata()
            if not metadata or not (reference := metadata.get('SpatialReference')):
                # raise warning
                continue
            if isinstance(reference, list):
                if len(reference) > 1:
                    # raise warning
                    continue
                reference = reference[0]

            reference = (Path(deriv_path) / reference).absolute()
            if not self.validate(item.path, str(reference)):
                # raise warning
                continue

            setattr(self, name, Path(item.path))
            self.references[name] = reference

    @property
    def mask(self) -> str | None:
        return self.t1w_mask or self.t2w_mask

    @property
    def aseg(self) -> str | None:
        return self.t1w_aseg or self.t2w_aseg

    @staticmethod
    def validate(derivative: str, reference: str, atol: float = 1e-5) -> bool:
        anat = nb.load(reference)
        expected_ort = nb.aff2axcodes(anat.affine)
        img = nb.load(derivative)
        if nb.aff2axcodes(img.affine) != expected_ort:
            return False
        if img.shape != anat.shape or not np.allclose(anat.affine, img.affine, atol=atol):
            return False
        return True


def write_bidsignore(deriv_dir):
    # TODO: Port to niworkflows
    bids_ignore = (
        "*.html",
        "logs/",
        "figures/",  # Reports
        "*_xfm.*",  # Unspecified transform files
        "*.surf.gii",  # Unspecified structural outputs
        # Unspecified functional outputs
        "*_boldref.nii.gz",
        "*_bold.func.gii",
        "*_mixing.tsv",
        "*_AROMAnoiseICs.csv",
        "*_timeseries.tsv",
    )
    ignore_file = Path(deriv_dir) / ".bidsignore"

    ignore_file.write_text("\n".join(bids_ignore) + "\n")


def write_derivative_description(bids_dir, deriv_dir):
    from ..__about__ import DOWNLOAD_URL, __packagename__, __version__

    bids_dir = Path(bids_dir)
    deriv_dir = Path(deriv_dir)
    desc = {
        "Name": "NiBabies: Neuroimaging preprocessing workflows for babies",
        "BIDSVersion": "1.4.0",
        "DatasetType": "derivative",
        "GeneratedBy": [
            {
                "Name": __packagename__,
                "Version": __version__,
                "CodeURL": DOWNLOAD_URL,
            }
        ],
        "HowToAcknowledge": "TODO",
    }

    # Keys that can only be set by environment
    if "NIBABIES_DOCKER_TAG" in os.environ:
        desc["GeneratedBy"][0]["Container"] = {
            "Type": "docker",
            "Tag": f"nipreps/nibabies:{os.environ['NIBABIES_DOCKER_TAG']}",
        }
    if "NIBABIES_SINGULARITY_URL" in os.environ:
        desc["GeneratedBy"][0]["Container"] = {
            "Type": "singularity",
            "URI": os.getenv("NIBABIES_SINGULARITY_URL"),
        }

    # Keys deriving from source dataset
    orig_desc = {}
    fname = bids_dir / "dataset_description.json"
    if fname.exists():
        orig_desc = json.loads(fname.read_text())

    if "DatasetDOI" in orig_desc:
        desc["SourceDatasets"] = [
            {"URL": f'https://doi.org/{orig_desc["DatasetDOI"]}', "DOI": orig_desc["DatasetDOI"]}
        ]
    if "License" in orig_desc:
        desc["License"] = orig_desc["License"]

    Path.write_text(deriv_dir / "dataset_description.json", json.dumps(desc, indent=4))


def extract_entities(file_list):
    """
    Return a dictionary of common entities given a list of files.
    Examples
    --------
    >>> extract_entities('sub-01/anat/sub-01_T1w.nii.gz')
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(['sub-01/anat/sub-01_T1w.nii.gz'] * 2)
    {'subject': '01', 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    >>> extract_entities(['sub-01/anat/sub-01_run-1_T1w.nii.gz',
    ...                   'sub-01/anat/sub-01_run-2_T1w.nii.gz'])
    {'subject': '01', 'run': [1, 2], 'suffix': 'T1w', 'datatype': 'anat', 'extension': '.nii.gz'}
    """
    from collections import defaultdict

    from bids.layout import parse_file_entities
    from niworkflows.utils.connections import listify

    entities = defaultdict(list)
    for e, v in [
        ev_pair for f in listify(file_list) for ev_pair in parse_file_entities(f).items()
    ]:
        entities[e].append(v)

    def _unique(inlist):
        inlist = sorted(set(inlist))
        if len(inlist) == 1:
            return inlist[0]
        return inlist

    return {k: _unique(v) for k, v in entities.items()}


def validate_input_dir(exec_env, bids_dir, participant_label):
    # Ignore issues and warnings that should not influence NiBabies
    import subprocess
    import tempfile

    validator_config_dict = {
        "ignore": [
            "EVENTS_COLUMN_ONSET",
            "EVENTS_COLUMN_DURATION",
            "TSV_EQUAL_ROWS",
            "TSV_EMPTY_CELL",
            "TSV_IMPROPER_NA",
            "VOLUME_COUNT_MISMATCH",
            "BVAL_MULTIPLE_ROWS",
            "BVEC_NUMBER_ROWS",
            "DWI_MISSING_BVAL",
            "INCONSISTENT_SUBJECTS",
            "INCONSISTENT_PARAMETERS",
            "BVEC_ROW_LENGTH",
            "B_FILE",
            "PARTICIPANT_ID_COLUMN",
            "PARTICIPANT_ID_MISMATCH",
            "TASK_NAME_MUST_DEFINE",
            "PHENOTYPE_SUBJECTS_MISSING",
            "STIMULUS_FILE_MISSING",
            "DWI_MISSING_BVEC",
            "EVENTS_TSV_MISSING",
            "TSV_IMPROPER_NA",
            "ACQTIME_FMT",
            "Participants age 89 or higher",
            "DATASET_DESCRIPTION_JSON_MISSING",
            "FILENAME_COLUMN",
            "WRONG_NEW_LINE",
            "MISSING_TSV_COLUMN_CHANNELS",
            "MISSING_TSV_COLUMN_IEEG_CHANNELS",
            "MISSING_TSV_COLUMN_IEEG_ELECTRODES",
            "UNUSED_STIMULUS",
            "CHANNELS_COLUMN_SFREQ",
            "CHANNELS_COLUMN_LOWCUT",
            "CHANNELS_COLUMN_HIGHCUT",
            "CHANNELS_COLUMN_NOTCH",
            "CUSTOM_COLUMN_WITHOUT_DESCRIPTION",
            "ACQTIME_FMT",
            "SUSPICIOUSLY_LONG_EVENT_DESIGN",
            "SUSPICIOUSLY_SHORT_EVENT_DESIGN",
            "MALFORMED_BVEC",
            "MALFORMED_BVAL",
            "MISSING_TSV_COLUMN_EEG_ELECTRODES",
            "MISSING_SESSION",
        ],
        "error": ["NO_T1W"],
        "ignoredFiles": ["/dataset_description.json", "/participants.tsv"],
    }
    # Limit validation only to data from requested participants
    if participant_label:
        all_subs = set([s.name[4:] for s in bids_dir.glob("sub-*")])
        selected_subs = set([s[4:] if s.startswith("sub-") else s for s in participant_label])
        bad_labels = selected_subs.difference(all_subs)
        if bad_labels:
            error_msg = (
                "Data for requested participant(s) label(s) not found. Could "
                "not find data for participant(s): %s. Please verify the requested "
                "participant labels."
            )
            if exec_env == "docker":
                error_msg += (
                    " This error can be caused by the input data not being "
                    "accessible inside the docker container. Please make sure all "
                    "volumes are mounted properly (see https://docs.docker.com/"
                    "engine/reference/commandline/run/#mount-volume--v---read-only)"
                )
            if exec_env == "singularity":
                error_msg += (
                    " This error can be caused by the input data not being "
                    "accessible inside the singularity container. Please make sure "
                    "all paths are mapped properly (see https://www.sylabs.io/"
                    "guides/3.0/user-guide/bind_paths_and_mounts.html)"
                )
            raise RuntimeError(error_msg % ",".join(bad_labels))

        ignored_subs = all_subs.difference(selected_subs)
        if ignored_subs:
            for sub in ignored_subs:
                validator_config_dict["ignoredFiles"].append("/sub-%s/**" % sub)
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".json") as temp:
        temp.write(json.dumps(validator_config_dict))
        temp.flush()
        try:
            subprocess.check_call(["bids-validator", str(bids_dir), "-c", temp.name])
        except FileNotFoundError:
            print("bids-validator does not appear to be installed", file=sys.stderr)


def parse_bids_for_age_months(
    bids_root: str | Path,
    subject_id: str,
    session_id: str | None = None,
) -> int | None:
    """
    Given a BIDS root, query the BIDS metadata files for participant age, in months.

    The heuristic followed is:
    1) Check `sub-<subject_id>/sub-<subject_id>_sessions.tsv`
    2) Check `<root>/participants.tsv`
    """
    age = None
    if subject_id.startswith('sub-'):
        subject_id = subject_id[4:]
    if session_id and session_id.startswith('ses-'):
        session_id = session_id[4:]

    sessions_tsv = Path(bids_root) / f'sub-{subject_id}' / f'sub-{subject_id}_sessions.tsv'
    if sessions_tsv.exists() and session_id is not None:
        age = _get_age_from_tsv(sessions_tsv, level='session', key=f'ses-{session_id}')

    participants_tsv = Path(bids_root) / 'participants.tsv'
    if participants_tsv.exists() and age is None:
        age = _get_age_from_tsv(participants_tsv, level='participant', key=f'sub-{subject_id}')

    return age


def _get_age_from_tsv(
    bids_tsv: Path, level: ty.Literal['session', 'participant'], key: str
) -> int | None:
    import pandas as pd

    df = pd.read_csv(str(bids_tsv), sep='\t')
    age_col = None
    # prefer explicit "age_months" over "age"
    for c in ('age_months', 'age'):
        if c in df.columns:
            age_col = c
            break

    if age_col == 'age':
        # verify age is in months
        bids_json = bids_tsv.with_suffix('.json')
        if not _verify_age_json(bids_json):
            warnings.warn(f'Could not verify age column is in months for file: {bids_tsv}')

    # find the relevant row
    if level == 'session':
        mask = df.session_id == key
    elif level == 'participant':
        mask = df.participant_id == key

    try:
        # extract age value from row
        age = int(df.loc[mask, age_col].values[0])
    except Exception:
        age = None
    return age


def _verify_age_json(bids_json: Path) -> bool:
    try:
        data = json.loads(bids_json.read_text())
        return data['age']['Units'].lower() == 'months'
    except Exception:
        return False
