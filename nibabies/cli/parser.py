# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Parser."""
import sys
from .. import config


def _build_parser():
    """Build parser object off fMRIPrep's."""
    from fmriprep.cli.parser import _build_parser as fmriprep_parser

    parser = fmriprep_parser()
    parser.description = f"NiBabies: Preprocessing workflows for infants v{config.environment.version}"

    # Change a few defaults
    for action in parser._actions:
        avars = vars(action)
        if avars.get('dest') == 'version':
            avars['version']= f"NiBabies v{config.environment.version}"
        elif avars.get('dest') == 'skull_strip_template':
            avars['default'] = 'UNCInfant:cohort-1'

    # Add new options
    g_baby = parser.add_argument_group("NiBabies specific options")
    g_baby.add_argument(
        "--anat-modality",
        default="t1w",
        choices=("t1w", "t2w"),
        help="Preferred modality to use for anatomical processing",
    )
    g_baby.add_argument(
        "--age-months",
        dest="age_months",
        type=int,
        help="Age in months",
    )
    g_baby.add_argument(
        "--segmentation-atlases-dir",
        help="Directory containing precalculated segmentations to use for JointLabelFusion."
    )
    return parser


def parse_args(args=None, namespace=None):
    """Parse args and run further checks on the command line."""
    import logging

    parser = _build_parser()
    opts = parser.parse_args(args, namespace)

    if opts.config_file:
        skip = {} if opts.reports_only else {"execution": ("run_uuid",)}
        config.load(opts.config_file, skip=skip)
        config.loggers.cli.info(f"Loaded previous configuration file {opts.config_file}")

    config.execution.log_level = int(max(25 - 5 * opts.verbose_count, logging.DEBUG))
    config.from_dict(vars(opts))

    # Initialize --output-spaces if not defined
    if config.execution.output_spaces is None:
        from niworkflows.utils.spaces import Reference, SpatialReferences

        config.execution.output_spaces = SpatialReferences(
            [Reference("UNCInfant", {"cohort": "1"})]
        )

    # Retrieve logging level
    build_log = config.loggers.cli

    # Load base plugin_settings from file if --use-plugin
    if opts.use_plugin is not None:
        import yaml

        with open(opts.use_plugin) as f:
            plugin_settings = yaml.load(f, Loader=yaml.FullLoader)
        _plugin = plugin_settings.get("plugin")
        if _plugin:
            config.nipype.plugin = _plugin
            config.nipype.plugin_args = plugin_settings.get("plugin_args", {})
            config.nipype.nprocs = opts.nprocs or config.nipype.plugin_args.get(
                "n_procs", config.nipype.nprocs
            )

    # Resource management options
    # Note that we're making strong assumptions about valid plugin args
    # This may need to be revisited if people try to use batch plugins
    if 1 < config.nipype.nprocs < config.nipype.omp_nthreads:
        build_log.warning(
            f"Per-process threads (--omp-nthreads={config.nipype.omp_nthreads}) exceed "
            f"total threads (--nthreads/--n_cpus={config.nipype.nprocs})"
        )

    # Inform the user about the risk of using brain-extracted images
    if config.workflow.skull_strip_t1w == "auto":
        build_log.warning(
            """\
Option ``--skull-strip-t1w`` was set to 'auto'. A heuristic will be \
applied to determine whether the input T1w image(s) have already been skull-stripped.
If that were the case, brain extraction and INU correction will be skipped for those T1w \
inputs. Please, BEWARE OF THE RISKS TO THE CONSISTENCY of results when using varying \
processing workflows across participants. To determine whether a participant has been run \
through the shortcut pipeline (meaning, brain extraction was skipped), please check the \
citation boilerplate. When reporting results with varying pipelines, please make sure you \
mention this particular variant of NiBabies listing the participants for which it was \
applied."""
        )

    bids_dir = config.execution.bids_dir
    output_dir = config.execution.output_dir
    work_dir = config.execution.work_dir
    version = config.environment.version
    output_layout = config.execution.output_layout

    if config.execution.fs_subjects_dir is None:
        if output_layout == "bids":
            config.execution.fs_subjects_dir = output_dir / "sourcedata" / "freesurfer-infant"
        elif output_layout == "legacy":
            config.execution.fs_subjects_dir = output_dir / "freesurfer-infant"
    if config.execution.fmriprep_dir is None:
        if output_layout == "bids":
            config.execution.fmriprep_dir = output_dir
        elif output_layout == "legacy":
            config.execution.fmriprep_dir = output_dir / "nibabies"

    # Wipe out existing work_dir
    if opts.clean_workdir and work_dir.exists():
        from niworkflows.utils.misc import clean_directory

        build_log.info(f"Clearing previous NiBabies working directory: {work_dir}")
        if not clean_directory(work_dir):
            build_log.warning(
                f"Could not clear all contents of working directory: {work_dir}"
            )

    # Ensure input and output folders are not the same
    if output_dir == bids_dir:
        parser.error(
            "The selected output folder is the same as the input BIDS folder. "
            "Please modify the output path (suggestion: %s)."
            % bids_dir
            / "derivatives"
            / ("nibabies-%s" % version.split("+")[0])
        )

    if bids_dir in work_dir.parents:
        parser.error(
            "The selected working directory is a subdirectory of the input BIDS folder. "
            "Please modify the output path."
        )

    # Validate inputs
    if not opts.skip_bids_validation:
        from ..utils.bids import validate_input_dir

        build_log.info(
            "Making sure the input data is BIDS compliant (warnings can be ignored in most "
            "cases)."
        )
        validate_input_dir(
            config.environment.exec_env, opts.bids_dir, opts.participant_label
        )

    # Setup directories
    config.execution.log_dir = config.execution.fmriprep_dir / "logs"
    # Check and create output and working directories
    config.execution.log_dir.mkdir(exist_ok=True, parents=True)
    work_dir.mkdir(exist_ok=True, parents=True)

    # Force initialization of the BIDSLayout
    config.execution.init()
    all_subjects = config.execution.layout.get_subjects()
    if config.execution.participant_label is None:
        config.execution.participant_label = all_subjects

    participant_label = set(config.execution.participant_label)
    missing_subjects = participant_label - set(all_subjects)
    if missing_subjects:
        parser.error(
            "One or more participant labels were not found in the BIDS directory: "
            "%s." % ", ".join(missing_subjects)
        )

    config.execution.participant_label = sorted(participant_label)
    config.workflow.skull_strip_template = config.workflow.skull_strip_template[0]
