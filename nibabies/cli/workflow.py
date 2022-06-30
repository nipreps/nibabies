"""
The workflow builder factory method.

All the checks and the construction of the workflow are done
inside this function that has pickleable inputs and output
dictionary (``retval``) to allow isolation using a
``multiprocessing.Process`` that allows nibabies to enforce
a hard-limited memory-scope.

"""


def build_workflow(config_file):
    """Create the Nipype Workflow that supports the whole execution graph."""
    from niworkflows.utils.bids import check_pipeline_version, collect_participants
    from niworkflows.utils.misc import check_valid_fs_license

    from .. import config
    from ..reports.core import generate_reports
    from ..utils.misc import check_deps
    from ..workflows.base import init_nibabies_wf

    # initalize config
    config.load(config_file)
    build_logger = config.loggers.workflow

    nibabies_dir = config.execution.nibabies_dir
    version = config.environment.version

    retval = {"return_code": 1, "workflow": None}

    # warn if older results exist: check for dataset_description.json in output folder
    msg = check_pipeline_version(version, nibabies_dir / "dataset_description.json")
    if msg is not None:
        build_logger.warning(msg)

    # Please note this is the input folder's dataset_description.json
    dset_desc_path = config.execution.bids_dir / "dataset_description.json"
    if dset_desc_path.exists():
        from hashlib import sha256

        desc_content = dset_desc_path.read_bytes()
        config.execution.bids_description_hash = sha256(desc_content).hexdigest()

    # First check that bids_dir looks like a BIDS folder
    subject_list = collect_participants(
        config.execution.layout, participant_label=config.execution.participant_label
    )
    subjects_sessions = {
        subject: config.execution.session_id
        or config.execution.layout.get_sessions(scope='raw', subject=subject)
        or [None]
        for subject in subject_list
    }

    # Called with reports only
    if config.execution.reports_only:
        from pkg_resources import resource_filename as pkgrf

        build_logger.log(25, "Running --reports-only on participants %s", ", ".join(subject_list))
        retval["return_code"] = generate_reports(
            subject_list,
            nibabies_dir,
            config.execution.run_uuid,
            config=pkgrf("nibabies", "data/reports-spec.yml"),
            packagename="nibabies",
        )
        return retval

    # Build main workflow
    init_msg = f"""
    Running nibabies version {config.environment.version}:
      * BIDS dataset path: {config.execution.bids_dir}.
      * Participant list: {subject_list}.
      * Run identifier: {config.execution.run_uuid}.
      * Output spaces: {config.execution.output_spaces}."""

    if config.execution.anat_derivatives:
        init_msg += f"""
      * Anatomical derivatives: {config.execution.anat_derivatives}."""

    if config.execution.fs_subjects_dir:
        init_msg += f"""
      * Pre-run FreeSurfer's SUBJECTS_DIR: {config.execution.fs_subjects_dir}."""
    build_logger.log(25, init_msg)

    retval["workflow"] = init_nibabies_wf(subjects_sessions)

    # Check for FS license after building the workflow
    if not check_valid_fs_license():
        build_logger.critical(
            """\
ERROR: a valid license file is required for FreeSurfer to run. nibabies looked for an existing \
license file at several paths, in this order: 1) command line argument ``--fs-license-file``; \
2) ``$FS_LICENSE`` environment variable; and 3) the ``$FREESURFER_HOME/license.txt`` path. Get it \
(for free) by registering at https://surfer.nmr.mgh.harvard.edu/registration.html"""
        )
        retval["return_code"] = 126  # 126 == Command invoked cannot execute.
        return retval

    # Check workflow for missing commands
    missing = check_deps(retval["workflow"])
    if missing:
        build_logger.critical(
            "Cannot run nibabies. Missing dependencies:%s",
            "\n\t* ".join([""] + [f"{cmd} (Interface: {iface})" for iface, cmd in missing]),
        )
        retval["return_code"] = 127  # 127 == command not found.
        return retval

    # config.to_filename(config_file)
    build_logger.info(
        "NiBabies workflow graph with %d nodes built successfully.",
        len(retval["workflow"]._get_all_nodes()),
    )
    retval["return_code"] = 0
    return retval


def build_boilerplate(workflow):
    """Write boilerplate in an isolated process."""
    from .. import config

    logs_path = config.execution.nibabies_dir / "logs"
    boilerplate = workflow.visit_desc()
    citation_files = {ext: logs_path / f"CITATION.{ext}" for ext in ("bib", "tex", "md", "html")}

    if boilerplate:
        # To please git-annex users and also to guarantee consistency
        # among different renderings of the same file, first remove any
        # existing one
        for citation_file in citation_files.values():
            try:
                citation_file.unlink()
            except FileNotFoundError:
                pass

    citation_files["md"].write_text(boilerplate)

    if not config.execution.md_only_boilerplate and citation_files["md"].exists():
        from shutil import copyfile
        from subprocess import CalledProcessError, TimeoutExpired, check_call

        from pkg_resources import resource_filename as pkgrf

        # Generate HTML file resolving citations
        cmd = [
            "pandoc",
            "-s",
            "--bibliography",
            pkgrf("nibabies", "data/boilerplate.bib"),
            "--citeproc",
            "--metadata",
            'pagetitle="nibabies citation boilerplate"',
            str(citation_files["md"]),
            "-o",
            str(citation_files["html"]),
        ]

        config.loggers.cli.info("Generating an HTML version of the citation boilerplate...")
        try:
            check_call(cmd, timeout=10)
        except (FileNotFoundError, CalledProcessError, TimeoutExpired):
            config.loggers.cli.warning("Could not generate CITATION.html file:\n%s", " ".join(cmd))

        # Generate LaTex file resolving citations
        cmd = [
            "pandoc",
            "-s",
            "--bibliography",
            pkgrf("nibabies", "data/boilerplate.bib"),
            "--natbib",
            str(citation_files["md"]),
            "-o",
            str(citation_files["tex"]),
        ]
        config.loggers.cli.info("Generating a LaTeX version of the citation boilerplate...")
        try:
            check_call(cmd, timeout=10)
        except (FileNotFoundError, CalledProcessError, TimeoutExpired):
            config.loggers.cli.warning("Could not generate CITATION.tex file:\n%s", " ".join(cmd))
        else:
            copyfile(pkgrf("nibabies", "data/boilerplate.bib"), citation_files["bib"])
