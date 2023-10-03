# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# STATEMENT OF CHANGES: This file is derived from sources licensed under the Apache-2.0 terms,
# and this file has been changed.
# The original file this work derives from is found at:
# https://github.com/nipreps/fmriprep/blob/a4fd718/fmriprep/workflows/bold/base.py
#
# [January 2022] CHANGES:
#   * `init_nibabies_wf` now takes in a dictionary composed of participant/session key/values.
#   * `init_single_subject_wf` now differentiates between participant sessions.
#     This change is to treat sessions as a "first-class" identifier, to better handle the
#     potential rapid changing of brain morphometry.
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
"""
NiBabies base processing workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_nibabies_wf
.. autofunction:: init_single_subject_wf

"""
from __future__ import annotations

import os
import sys
import typing as ty
from copy import deepcopy

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from packaging.version import Version

from nibabies import config
from nibabies.interfaces import DerivativesDataSink
from nibabies.interfaces.reports import AboutSummary, SubjectSummary
from nibabies.utils.bids import parse_bids_for_age_months
from nibabies.workflows.bold import init_func_preproc_wf

if ty.TYPE_CHECKING:
    from niworkflows.utils.spaces import SpatialReferences


def init_nibabies_wf(subworkflows_list):
    """
    Build *NiBabies*'s pipeline.

    This workflow organizes the execution of NiBabies, with a sub-workflow for
    each subject.

    If FreeSurfer's ``infant_recon_all`` is to be run, a corresponding folder is created
    and populated with any needed template subjects under the derivatives folder.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from nibabies.workflows.tests import mock_config
            from nibabies.workflows.base import init_nibabies_wf
            with mock_config():
                wf = init_nibabies_wf()

    Parameters
    ----------
    subworkflows_list: :obj:`list` of :obj:`tuple`
        A list of the subworkflows to create.
        Each subject session is run as an individual workflow.
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.bids import BIDSFreeSurferDir

    ver = Version(config.environment.version)
    nibabies_wf = Workflow(name=f"nibabies_{ver.major}_{ver.minor}_wf")
    nibabies_wf.base_dir = config.execution.work_dir

    execution_spaces = init_execution_spaces()

    freesurfer = config.workflow.run_reconall
    if freesurfer:
        fsdir = pe.Node(
            BIDSFreeSurferDir(
                derivatives=config.execution.output_dir,
                freesurfer_home=os.getenv("FREESURFER_HOME"),
                spaces=execution_spaces.get_fs_spaces(),
            ),
            name=f"fsdir_run_{config.execution.run_uuid.replace('-', '_')}",
            run_without_submitting=True,
        )
        if config.execution.fs_subjects_dir is not None:
            fsdir.inputs.subjects_dir = str(config.execution.fs_subjects_dir.absolute())

    for subject_id, session_id in subworkflows_list:
        # Calculate the age and age-specific spaces
        age = parse_bids_for_age_months(config.execution.bids_dir, subject_id, session_id)
        if config.workflow.age_months:
            config.loggers.cli.warning(
                "`--age-months` is deprecated and will be removed in a future release."
                "Please use a `sessions.tsv` or `participants.tsv` file to track participants age."
            )
            age = config.workflow.age_months
        if age is None:
            raise RuntimeError(
                "Could not find age for sub-{subject}{session}".format(
                    subject=subject_id, session=f'_ses-{session_id}' if session_id else ''
                )
            )
        output_spaces = init_workflow_spaces(execution_spaces, age)

        # skull strip template cohort
        single_subject_wf = init_single_subject_wf(
            subject_id,
            session_id=session_id,
            age=age,
            spaces=output_spaces,
        )

        bids_level = [f"sub-{subject_id}"]
        if session_id:
            bids_level.append(f"ses-{session_id}")

        log_dir = (
            config.execution.nibabies_dir.joinpath(*bids_level) / "log" / config.execution.run_uuid
        )

        single_subject_wf.config["execution"]["crashdump_dir"] = str(log_dir)
        for node in single_subject_wf._get_all_nodes():
            node.config = deepcopy(single_subject_wf.config)
        if freesurfer:
            nibabies_wf.connect(fsdir, "subjects_dir", single_subject_wf, "inputnode.subjects_dir")
        else:
            nibabies_wf.add_nodes([single_subject_wf])

        # Dump a copy of the config file into the log directory
        log_dir.mkdir(exist_ok=True, parents=True)
        config.to_filename(log_dir / "nibabies.toml")

    return nibabies_wf


def init_single_subject_wf(
    subject_id: str,
    session_id: str | None = None,
    age: int | None = None,
    spaces: SpatialReferences | None = None,
):
    """
    Organize the preprocessing pipeline for a single subject, at a single session.

    It collects and reports information about the subject, and prepares
    sub-workflows to perform anatomical and functional preprocessing.
    Anatomical preprocessing is performed in a single workflow, regardless of
    the number of sessions.
    Functional preprocessing is performed using a separate workflow for each
    individual BOLD series.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from nibabies.workflows.tests import mock_config
            from nibabies.workflows.base import init_single_subject_wf
            with mock_config():
                wf = init_single_subject_wf('01')

    Parameters
    ----------
    subject_id : :obj:`str`
        Subject label for this single-subject workflow.
    session_id : :obj:`str` or None
        Session identifier.
    age: :obj:`int` or None
        Age (in months) of subject.

    Inputs
    ------
    subjects_dir : :obj:`str`
        FreeSurfer's ``$SUBJECTS_DIR``.

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.bids import BIDSDataGrabber, BIDSInfo
    from niworkflows.interfaces.nilearn import NILEARN_VERSION
    from niworkflows.utils.bids import collect_data
    from niworkflows.utils.spaces import Reference

    from ..utils.bids import Derivatives
    from ..utils.misc import fix_multi_source_name
    from .anatomical import init_infant_anat_wf, init_infant_single_anat_wf

    name = (
        f"single_subject_{subject_id}_{session_id}_wf"
        if session_id
        else f"single_subject_{subject_id}_wf"
    )
    subject_data = collect_data(
        config.execution.layout,
        subject_id,
        session_id=session_id,
        task=config.execution.task_id,
        echo=config.execution.echo_idx,
        bids_filters=config.execution.bids_filters,
    )[0]

    if "flair" in config.workflow.ignore:
        subject_data["flair"] = []
    if "t2w" in config.workflow.ignore:
        subject_data["t2w"] = []

    anat_only = config.workflow.anat_only
    derivatives = Derivatives(bids_root=config.execution.layout.root)
    contrast = "T1w" if subject_data["t1w"] else "T2w"
    single_modality = not (subject_data['t1w'] and subject_data['t2w'])
    # Make sure we always go through these two checks
    if not anat_only and not subject_data["bold"]:
        task_id = config.execution.task_id
        raise RuntimeError(
            "No BOLD images found for participant {} and task {}. "
            "All workflows require BOLD images.".format(
                subject_id, task_id if task_id else "<all>"
            )
        )

    if config.execution.derivatives:
        for deriv_path in config.execution.derivatives:
            config.loggers.workflow.info("Searching for derivatives in %s", deriv_path)
            derivatives.populate(
                deriv_path,
                subject_id,
                session_id=session_id,
            )
        config.loggers.workflow.info("Found precomputed derivatives %s", derivatives)

    workflow = Workflow(name=name)
    workflow.__desc__ = """
Results included in this manuscript come from preprocessing
performed using *NiBabies* {nibabies_ver},
derived from fMRIPrep (@fmriprep1; @fmriprep2; RRID:SCR_016216).
The underlying workflow engine used is *Nipype* {nipype_ver}
(@nipype1; @nipype2; RRID:SCR_002502).

""".format(
        nibabies_ver=config.environment.version,
        nipype_ver=config.environment.nipype_version,
    )
    workflow.__postdesc__ = """

Many internal operations of *NiBabies* use
*Nilearn* {nilearn_ver} [@nilearn, RRID:SCR_001362],
mostly within the functional processing workflow.
For more details of the pipeline, see [the section corresponding
to workflows in *nibabies*'s documentation]\
(https://nibabies.readthedocs.io/en/latest/workflows.html \
"NiBabies's documentation").


### Copyright Waiver

The above boilerplate text was automatically generated by NiBabies
with the express intention that users should copy and paste this
text into their manuscripts *unchanged*.
It is released under the [CC0]\
(https://creativecommons.org/publicdomain/zero/1.0/) license.

### References

""".format(
        nilearn_ver=NILEARN_VERSION
    )

    nibabies_dir = str(config.execution.nibabies_dir)

    inputnode = pe.Node(niu.IdentityInterface(fields=["subjects_dir"]), name="inputnode")

    # TODO: Revisit T1w/T2w restrictions for BIDSDataGrabber
    bidssrc = pe.Node(
        BIDSDataGrabber(
            subject_data=subject_data,
            anat_only=anat_only,
            anat_derivatives=False,
            subject_id=subject_id,
        ),
        name="bidssrc",
    )

    bids_info = pe.Node(
        BIDSInfo(bids_dir=config.execution.bids_dir, bids_validate=False),
        name="bids_info",
    )

    summary = pe.Node(
        SubjectSummary(
            std_spaces=spaces.get_spaces(nonstandard=False),
            nstd_spaces=spaces.get_spaces(standard=False),
        ),
        name="summary",
        run_without_submitting=True,
    )

    about = pe.Node(
        AboutSummary(version=config.environment.version, command=" ".join(sys.argv)),
        name="about",
        run_without_submitting=True,
    )

    ds_report_summary = pe.Node(
        DerivativesDataSink(
            base_directory=nibabies_dir,
            desc="summary",
            datatype="figures",
            dismiss_entities=("echo",),
        ),
        name="ds_report_summary",
        run_without_submitting=True,
    )

    ds_report_about = pe.Node(
        DerivativesDataSink(
            base_directory=nibabies_dir,
            desc="about",
            datatype="figures",
            dismiss_entities=("echo",),
        ),
        name="ds_report_about",
        run_without_submitting=True,
    )

    wf_args = dict(
        ants_affine_init=True,
        age_months=age,
        contrast=contrast,
        t1w=subject_data["t1w"],
        t2w=subject_data["t2w"],
        bids_root=config.execution.bids_dir,
        derivatives=derivatives,
        freesurfer=config.workflow.run_reconall,
        hires=config.workflow.hires,
        longitudinal=config.workflow.longitudinal,
        omp_nthreads=config.nipype.omp_nthreads,
        output_dir=nibabies_dir,
        segmentation_atlases=config.execution.segmentation_atlases_dir,
        skull_strip_mode=config.workflow.skull_strip_t1w,
        skull_strip_template=Reference.from_string(config.workflow.skull_strip_template)[0],
        sloppy=config.execution.sloppy,
        spaces=spaces,
        cifti_output=config.workflow.cifti_output,
    )
    anat_preproc_wf = (
        init_infant_anat_wf(**wf_args)
        if not single_modality
        else init_infant_single_anat_wf(**wf_args)
    )

    # fmt: off
    workflow.connect([
        (inputnode, anat_preproc_wf, [
            ('subjects_dir', 'inputnode.subjects_dir'),
        ]),
        (inputnode, summary, [
            ('subjects_dir', 'subjects_dir'),
        ]),
        (bidssrc, summary, [
            ('bold', 'bold'),
        ]),
        (bids_info, summary, [
            ('subject', 'subject_id'),
        ]),
        (bids_info, anat_preproc_wf, [
            (('subject', _prefix), 'inputnode.subject_id'),
        ]),
        (bidssrc, anat_preproc_wf, [
            ('t1w', 'inputnode.t1w'),
            ('t2w', 'inputnode.t2w'),
        ]),
        (summary, ds_report_summary, [
            ('out_report', 'in_file'),
        ]),
        (about, ds_report_about, [
            ('out_report', 'in_file'),
        ]),
    ])

    workflow.connect([
        (bidssrc, bids_info, [
            ((contrast.lower(), fix_multi_source_name), 'in_file'),
        ]),
        (bidssrc, summary, [
            ('t1w', 't1w'),
            ('t2w', 't2w'),
        ]),
        (bidssrc, ds_report_summary, [
            ((contrast.lower(), fix_multi_source_name), 'source_file'),
        ]),
        (bidssrc, ds_report_about, [
            ((contrast.lower(), fix_multi_source_name), 'source_file'),
        ]),
    ])
    # fmt: on

    # Overwrite ``out_path_base`` of smriprep's DataSinks
    for node in workflow.list_node_names():
        if node.split(".")[-1].startswith("ds_"):
            workflow.get_node(node).interface.out_path_base = ""

    if anat_only:
        return workflow

    # Susceptibility distortion correction
    fmap_estimators = None
    if any((config.workflow.use_syn_sdc, config.workflow.force_syn)):
        config.loggers.workflow.critical("SyN processing is not yet implemented.")

    if "fieldmaps" not in config.workflow.ignore:
        from sdcflows.utils.wrangler import find_estimators

        # SDC Step 1: Run basic heuristics to identify available data for fieldmap estimation
        # For now, no fmapless
        fmap_estimators = find_estimators(
            layout=config.execution.layout,
            subject=subject_id,
            sessions=[session_id],
            fmapless=False,  # config.workflow.use_syn,
            force_fmapless=False,  # config.workflow.force_syn,
        )

    # Append the functional section to the existing anatomical exerpt
    # That way we do not need to stream down the number of bold datasets
    anat_preproc_wf.__postdesc__ = getattr(anat_preproc_wf, '__postdesc__') or ''
    func_pre_desc = f"""

Functional data preprocessing

: For each of the {len(subject_data['bold'])} BOLD runs found per subject (across all
tasks and sessions), the following preprocessing was performed."""

    func_preproc_wfs = []
    has_fieldmap = bool(fmap_estimators)
    for bold_file in subject_data['bold']:
        func_preproc_wf = init_func_preproc_wf(bold_file, spaces, has_fieldmap=has_fieldmap)
        if func_preproc_wf is None:
            continue

        func_preproc_wf.__desc__ = func_pre_desc + (getattr(func_preproc_wf, '__desc__') or '')
        # fmt:off
        workflow.connect([
            (anat_preproc_wf, func_preproc_wf, [
                ('outputnode.anat_preproc', 'inputnode.anat_preproc'),
                ('outputnode.anat_mask', 'inputnode.anat_mask'),
                ('outputnode.anat_brain', 'inputnode.anat_brain'),
                ('outputnode.anat_dseg', 'inputnode.anat_dseg'),
                ('outputnode.anat_aseg', 'inputnode.anat_aseg'),
                ('outputnode.anat_aparc', 'inputnode.anat_aparc'),
                ('outputnode.anat_tpms', 'inputnode.anat_tpms'),
                ('outputnode.template', 'inputnode.template'),
                ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
                ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
                # Undefined if --fs-no-reconall, but this is safe
                ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
                ('outputnode.subject_id', 'inputnode.subject_id'),
                ('outputnode.t1w2fsnative_xfm', 'inputnode.t1w2fsnative_xfm'),
                ('outputnode.fsnative2t1w_xfm', 'inputnode.fsnative2t1w_xfm'),
                ('outputnode.surfaces', 'inputnode.surfaces'),
                ('outputnode.morphometrics', 'inputnode.morphometrics'),
                ('outputnode.anat_ribbon', 'inputnode.anat_ribbon'),
                ('outputnode.sphere_reg_fsLR', 'inputnode.sphere_reg_fsLR'),
                ('outputnode.midthickness_fsLR', 'inputnode.midthickness_fsLR'),
            ]),
        ])
        # fmt:on
        func_preproc_wfs.append(func_preproc_wf)

    if not has_fieldmap:
        config.loggers.workflow.warning(
            "Data for fieldmap estimation not present. Please note that these data "
            "will not be corrected for susceptibility distortions."
        )
        return workflow

    config.loggers.workflow.info(
        f"Fieldmap estimators found: {[e.method for e in fmap_estimators]}"
    )

    from sdcflows import fieldmaps as fm
    from sdcflows.workflows.base import init_fmap_preproc_wf

    fmap_wf = init_fmap_preproc_wf(
        sloppy=bool(config.execution.sloppy),
        debug="fieldmaps" in config.execution.debug,
        estimators=fmap_estimators,
        omp_nthreads=config.nipype.omp_nthreads,
        output_dir=nibabies_dir,
        subject=subject_id,
    )
    fmap_wf.__desc__ = f"""

Preprocessing of B<sub>0</sub> inhomogeneity mappings

: A total of {len(fmap_estimators)} fieldmaps were found available within the input
BIDS structure for this particular subject.
"""

    for func_preproc_wf in func_preproc_wfs:
        # fmt: off
        workflow.connect([
            (fmap_wf, func_preproc_wf, [
                ("outputnode.fmap", "inputnode.fmap"),
                ("outputnode.fmap_ref", "inputnode.fmap_ref"),
                ("outputnode.fmap_coeff", "inputnode.fmap_coeff"),
                ("outputnode.fmap_mask", "inputnode.fmap_mask"),
                ("outputnode.fmap_id", "inputnode.fmap_id"),
                ("outputnode.method", "inputnode.sdc_method"),
            ]),
        ])
        # fmt: on

    # Overwrite ``out_path_base`` of sdcflows's DataSinks
    for node in fmap_wf.list_node_names():
        if node.split(".")[-1].startswith("ds_"):
            fmap_wf.get_node(node).interface.out_path_base = ""

    # Step 3: Manually connect PEPOLAR
    for estimator in fmap_estimators:
        config.loggers.workflow.info(
            f"""\
Setting-up fieldmap "{estimator.bids_id}" ({estimator.method}) with \
<{', '.join(s.path.name for s in estimator.sources)}>"""
        )
        if estimator.method in (fm.EstimatorType.MAPPED, fm.EstimatorType.PHASEDIFF):
            continue

        suffices = [s.suffix for s in estimator.sources]

        if estimator.method == fm.EstimatorType.PEPOLAR:
            if set(suffices) == {"epi"} or sorted(suffices) == ["bold", "epi"]:
                fmap_wf_inputs = getattr(fmap_wf.inputs, f"in_{estimator.bids_id}")
                fmap_wf_inputs.in_data = [str(s.path) for s in estimator.sources]
                fmap_wf_inputs.metadata = [s.metadata for s in estimator.sources]
            else:
                raise NotImplementedError(
                    "Sophisticated PEPOLAR schemes (e.g., using DWI+EPI) are unsupported."
                )

    return workflow


def _prefix(subid):
    return subid if subid.startswith("sub-") else f"sub-{subid}"


def init_workflow_spaces(execution_spaces: SpatialReferences, age_months: int):
    """
    Create output spaces at a per-subworkflow level.

    This address the case where a multi-session subject is run, and requires separate template cohorts.
    """
    from niworkflows.utils.spaces import Reference

    from nibabies.utils.misc import cohort_by_months

    spaces = deepcopy(execution_spaces)

    if age_months is None:
        raise RuntimeError("Participant age (in months) is required.")

    if not spaces.references:
        # Ensure age specific template is added if nothing is present
        cohort = cohort_by_months("MNIInfant", age_months)
        spaces.add(("MNIInfant", {"res": "native", "cohort": cohort}))

    if not spaces.is_cached():
        spaces.checkpoint()

    # Ensure user-defined spatial references for outputs are correctly parsed.
    # Certain options require normalization to a space not explicitly defined by users.
    # These spaces will not be included in the final outputs.
    if config.workflow.use_aroma:
        # Make sure there's a normalization to FSL for AROMA to use.
        spaces.add(Reference("MNI152NLin6Asym", {"res": "2"}))

    if config.workflow.cifti_output:
        # CIFTI grayordinates to corresponding FSL-MNI resolutions.
        vol_res = "2" if config.workflow.cifti_output == "91k" else "1"
        spaces.add(Reference("MNI152NLin6Asym", {"res": vol_res}))
        # Ensure a non-native version of MNIInfant is added as a target
        cohort = cohort_by_months("MNIInfant", age_months)
        spaces.add(Reference("MNIInfant", {"cohort": cohort}))

    return spaces


def init_execution_spaces():
    from niworkflows.utils.spaces import Reference, SpatialReferences

    spaces = config.execution.output_spaces or SpatialReferences()
    if not isinstance(spaces, SpatialReferences):
        spaces = SpatialReferences(
            [ref for s in spaces.split(" ") for ref in Reference.from_string(s)]
        )
    return spaces
