# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Utilities to handle BIDS inputs."""
import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import IO, List, Union


@dataclass
class BOLDGrouping:
    """This class is used to facilitate the grouping of BOLD series."""

    session: Union[str, None]
    pe_dir: str
    readout: float
    multiecho_id: str = None
    files: List[IO] = field(default_factory=list)

    @property
    def name(self) -> str:
        return f"{self.session}-{self.pe_dir}-{self.readout}-{self.multiecho_id}"

    def add_file(self, fl) -> None:
        self.files.append(fl)


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


def group_bolds_ref(*, layout, subject, sessions=None):
    """
    Extracts BOLD files from a BIDS dataset and combines them into buckets.
    Files in a bucket share:
    1) Session
    2) Phase-encoding direction (PEdir)
    3) Total readout time (TRT)

    Any files with missing data for (2) or (3) are put in their own bucket.

    Parameters
    ----------
    layout : pybids.layout.BIDSLayout
        Initialized BIDSLayout
    subject : str
        The subject ID
    sessions : None

    Outputs
    -------
    combinations : list of tuples
        Each tuple is composed of (session, PEdir, TRT)
    files : list of lists
        Files matching each combination.

    Limitations
    -----------
    Single-band reference (sbref) are excluded.
    """
    import re
    from contextlib import suppress
    from itertools import product

    from sdcflows.utils.epimanip import get_trt

    base_entities = {
        "subject": subject,
        "extension": (".nii", ".nii.gz"),
        "scope": "raw",  # Ensure derivatives are not captured
    }
    # dictionary containing unique Groupings and files
    groupings = {}
    # list of all BOLDS encountered
    all_bolds = []

    sessions = sessions if sessions else layout.get_sessions(subject=subject, scope="raw")

    for ses, suffix in sorted(product(sessions or (None,), {"bold"})):
        # bold files same session
        bolds = layout.get(suffix=suffix, session=ses, **base_entities)
        # some sessions may not have BOLD scans
        if bolds is None:
            continue

        for i, bold in enumerate(bolds):
            multiecho_id = None
            # multi-echo should be grouped together
            if 'echo' in bold.entities:
                # create unique id by dropping "_echo-{i}"
                multiecho_id = re.sub(r"_echo-\d+", "", bold.filename)

            # session, pe, ro
            meta = bold.get_metadata()
            pe_dir = meta.get("PhaseEncodingDirection")

            ro = None
            with suppress(ValueError):
                ro = get_trt(meta, bold.path)
            if ro is not None:
                meta.update({"TotalReadoutTime": ro})

            grouping = BOLDGrouping(
                session=ses,
                pe_dir=pe_dir,
                readout=ro,
                multiecho_id=multiecho_id,
            )

            if any(v is None for v in (pe_dir, ro)):
                # cannot be certain so treat as unique
                grouping.add_file(bold.path)
                groupings[f'unknown{i}'] = grouping
            else:
                try:
                    grouping = groupings[grouping.name]
                except KeyError:
                    groupings[grouping.name] = grouping

                grouping.add_file(bold.path)

        all_bolds += bolds

        if len(all_bolds) != sum([len(g.files) for _, g in groupings.items()]):
            msg = f"""Error encountered when grouping BOLD runs.
Combinations: {groupings}
BOLD files: {bolds}

Please file a bug-report with the nibabies developers at:
https://github.com/nipreps/nibabies/issues/new/choose
"""
            raise RuntimeError(msg)

    return groupings


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


def collect_precomputed_derivatives(layout, subject_id, derivatives_filters=None):
    """
    Query and collect precomputed derivatives.

    This function is used to determine which workflow steps can be skipped,
    based on the files found.
    """

    deriv_queries = {
        'anat_mask': {
            'datatype': 'anat',
            'desc': 'brain',
            'space': 'orig',
            'suffix': 'mask',
        },
        'anat_aseg': {
            'datatype': 'anat',
            'desc': 'aseg',
            'space': 'orig',
            'suffix': 'dseg',
        },
    }
    if derivatives_filters is not None:
        deriv_queries.update(derivatives_filters)

    derivatives = {}
    for deriv, query in deriv_queries.items():
        res = layout.get(
            scope='derivatives',
            subject=subject_id,
            extension=['.nii', '.nii.gz'],
            return_type="filename",
            **query,
        )
        if not res:
            continue
        if len(res) > 1:  # Some queries may want multiple results
            raise Exception(
                f"When searching for <{deriv}>, found multiple results: {[f.path for f in res]}"
            )
        derivatives[deriv] = res[0]
    return derivatives
