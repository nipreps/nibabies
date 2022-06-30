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
# Copyright 2021 The NiPreps Developers <nipreps@gmail.com>
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

import os
import sys
from copy import deepcopy

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from packaging.version import Version

from .. import config
from ..interfaces import DerivativesDataSink
from ..interfaces.reports import AboutSummary, SubjectSummary
from ..utils.bids import group_bolds_ref
from .bold import init_func_preproc_wf


def init_nibabies_wf(participants_table):
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
    participants_table: :obj:`dict`
        Keys of participant labels and values of the sessions to process.
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.bids import BIDSFreeSurferDir

    ver = Version(config.environment.version)
    nibabies_wf = Workflow(name=f"nibabies_{ver.major}_{ver.minor}_wf")
    nibabies_wf.base_dir = config.execution.work_dir

    freesurfer = config.workflow.run_reconall
    if freesurfer:
        fsdir = pe.Node(
            BIDSFreeSurferDir(
                derivatives=config.execution.output_dir,
                freesurfer_home=os.getenv("FREESURFER_HOME"),
                spaces=config.workflow.spaces.get_fs_spaces(),
            ),
            name=f"fsdir_run_{config.execution.run_uuid.replace('-', '_')}",
            run_without_submitting=True,
        )
        if config.execution.fs_subjects_dir is not None:
            fsdir.inputs.subjects_dir = str(config.execution.fs_subjects_dir.absolute())

    for subject_id, sessions in participants_table.items():
        for session_id in sessions:
            single_subject_wf = init_single_subject_wf(subject_id, session_id=session_id)

            bids_level = [f"sub-{subject_id}"]
            if session_id:
                bids_level.append(f"ses-{session_id}")

            log_dir = (
                config.execution.nibabies_dir.joinpath(*bids_level)
                / "log"
                / config.execution.run_uuid
            )

            single_subject_wf.config["execution"]["crashdump_dir"] = str(log_dir)
            for node in single_subject_wf._get_all_nodes():
                node.config = deepcopy(single_subject_wf.config)
            if freesurfer:
                nibabies_wf.connect(
                    fsdir, "subjects_dir", single_subject_wf, "inputnode.subjects_dir"
                )
            else:
                nibabies_wf.add_nodes([single_subject_wf])

            # Dump a copy of the config file into the log directory
            log_dir.mkdir(exist_ok=True, parents=True)
            config.to_filename(log_dir / "nibabies.toml")

    return nibabies_wf


def init_single_subject_wf(subject_id, session_id=None):
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

    from ..utils.misc import fix_multi_source_name
    from .anatomical import init_infant_anat_wf

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
    derivatives = config.execution.derivatives or {}
    anat_modality = "t1w" if subject_data["t1w"] else "t2w"
    spaces = config.workflow.spaces
    # Make sure we always go through these two checks
    if not anat_only and not subject_data["bold"]:
        task_id = config.execution.task_id
        raise RuntimeError(
            "No BOLD images found for participant {} and task {}. "
            "All workflows require BOLD images.".format(
                subject_id, task_id if task_id else "<all>"
            )
        )

    if derivatives:
        from ..utils.bids import collect_precomputed_derivatives

        derivatives = collect_precomputed_derivatives(
            config.execution.layout,
            subject_id,
            derivatives_filters=config.execution.derivatives_filters,
            # session_id=None,  # TODO: Ensure session is visible at workflow level
        )
        config.loggers.workflow.info(f"Found precomputed derivatives: {derivatives}")

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

    # Preprocessing of anatomical (includes registration to UNCInfant)
    anat_preproc_wf = init_infant_anat_wf(
        ants_affine_init=True,
        age_months=config.workflow.age_months,
        anat_modality=anat_modality,
        t1w=subject_data["t1w"],
        t2w=subject_data["t2w"],
        bids_root=config.execution.bids_dir,
        existing_derivatives=derivatives,
        freesurfer=config.workflow.run_reconall,
        longitudinal=config.workflow.longitudinal,
        omp_nthreads=config.nipype.omp_nthreads,
        output_dir=nibabies_dir,
        segmentation_atlases=config.execution.segmentation_atlases_dir,
        skull_strip_mode=config.workflow.skull_strip_t1w,
        skull_strip_template=Reference.from_string(config.workflow.skull_strip_template)[0],
        sloppy=config.execution.sloppy,
        spaces=spaces,
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
            (('t1w', fix_multi_source_name), 'in_file'),
        ]),
        (bidssrc, summary, [
            ('t1w', 't1w'),
            ('t2w', 't2w'),
        ]),
        (bidssrc, ds_report_summary, [
            (('t1w', fix_multi_source_name), 'source_file'),
        ]),
        (bidssrc, ds_report_about, [
            (('t1w', fix_multi_source_name), 'source_file'),
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
    anat_preproc_wf.__postdesc__ = (
        (anat_preproc_wf.__postdesc__ if hasattr(anat_preproc_wf, "__postdesc__") else "")
        + f"""

Functional data preprocessing

: For each of the {len(subject_data['bold'])} BOLD runs found per subject (across all
tasks and sessions), the following preprocessing was performed.
"""
    )

    # calculate reference image(s) for BOLD images
    # group all BOLD files based on same:
    # 1) session
    # 2) PE direction
    # 3) total readout time
    from niworkflows.workflows.epi.refmap import init_epi_reference_wf

    bold_groupings = group_bolds_ref(
        layout=config.execution.layout,
        subject=subject_id,
        sessions=[session_id],
    )

    func_preproc_wfs = []
    has_fieldmap = bool(fmap_estimators)
    for idx, grouping in enumerate(bold_groupings.values()):
        bold_ref_wf = init_epi_reference_wf(
            auto_bold_nss=True,
            name=f"bold_reference_wf{idx}",
            omp_nthreads=config.nipype.omp_nthreads,
        )
        bold_files = grouping.files
        bold_ref_wf.inputs.inputnode.in_files = grouping.files

        if grouping.multiecho_id is not None:
            bold_files = [bold_files]
        for idx, bold_file in enumerate(bold_files):
            func_preproc_wf = init_func_preproc_wf(
                bold_file,
                has_fieldmap=has_fieldmap,
                existing_derivatives=derivatives,
            )
            # fmt: off
            workflow.connect([
                (bold_ref_wf, func_preproc_wf, [
                    ('outputnode.epi_ref_file', 'inputnode.bold_ref'),
                    (
                        ('outputnode.xfm_files', _select_iter_idx, idx),
                        'inputnode.bold_ref_xfm'),
                    (
                        ('outputnode.n_dummy', _select_iter_idx, idx),
                        'inputnode.n_dummy_scans'),
                ]),
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
                    ('outputnode.anat2fsnative_xfm', 'inputnode.anat2fsnative_xfm'),
                    ('outputnode.fsnative2anat_xfm', 'inputnode.fsnative2anat_xfm'),
                ]),
            ])
            # fmt: on
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

                flatten = fmap_wf.get_node(f"wf_{estimator.bids_id}.flatten")
                flatten.inputs.max_trs = config.workflow.topup_max_vols
            else:
                raise NotImplementedError(
                    "Sophisticated PEPOLAR schemes (e.g., using DWI+EPI) are unsupported."
                )

    return workflow


def _prefix(subid):
    return subid if subid.startswith("sub-") else f"sub-{subid}"


def _select_iter_idx(in_list, idx):
    """Returns a specific index of a list/tuple"""
    if isinstance(in_list, (tuple, list)):
        return in_list[idx]
    raise AttributeError(f"Input {in_list} is incompatible type: {type(in_list)}")
