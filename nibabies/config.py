# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright The NiPreps Developers <nipreps@gmail.com>
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
r"""
A Python module to maintain unique, run-wide *NiBabies* settings.

This module implements the memory structures to keep a consistent, singleton config.
Settings are passed across processes via filesystem, and a copy of the settings for
each run and subject is left under
``<nibabies_dir>/sub-<participant_id>/log/<run_unique_id>/nibabies.toml``.
Settings are stored using :abbr:`ToML (Tom's Markup Language)`.
The module has a :py:func:`~nibabies.config.to_filename` function to allow writing out
the settings to hard disk in *ToML* format, which looks like:

.. literalinclude:: ../nibabies/data/tests/config.toml
   :language: toml
   :name: nibabies.toml
   :caption: **Example file representation of nibabies settings**.

This config file is used to pass the settings across processes,
using the :py:func:`~nibabies.config.load` function.

Configuration sections
----------------------
.. autoclass:: environment
   :members:
.. autoclass:: execution
   :members:
.. autoclass:: workflow
   :members:
.. autoclass:: nipype
   :members:

Usage
-----
A config file is used to pass settings and collect information as the execution
graph is built across processes.

.. code-block:: Python

    from nibabies import config
    config_file = config.execution.work_dir / '.nibabies.toml'
    config.to_filename(config_file)
    # Call build_workflow(config_file, retval) in a subprocess
    with Manager() as mgr:
        from .workflow import build_workflow
        retval = mgr.dict()
        p = Process(target=build_workflow, args=(str(config_file), retval))
        p.start()
        p.join()
    config.load(config_file)
    # Access configs from any code section as:
    value = config.section.setting

Logging
-------
.. autoclass:: loggers
   :members:

Other responsibilities
----------------------
The :py:mod:`config` is responsible for other conveniency actions.

  * Switching Python's :obj:`multiprocessing` to *forkserver* mode.
  * Set up a filter for warnings as early as possible.
  * Automated I/O magic operations. Some conversions need to happen in the
    store/load processes (e.g., from/to :obj:`~pathlib.Path` \<-\> :obj:`str`,
    :py:class:`~bids.layout.BIDSLayout`, etc.)

"""

import os
import typing as ty
from multiprocessing import set_start_method

from templateflow.conf import TF_LAYOUT

# Disable NiPype etelemetry always
_disable_et = bool(os.getenv('NO_ET') is not None or os.getenv('NIPYPE_NO_ET') is not None)
os.environ['NIPYPE_NO_ET'] = '1'
os.environ['NO_ET'] = '1'

CONFIG_FILENAME = 'nibabies.toml'

try:
    set_start_method('forkserver')
except RuntimeError:
    pass  # context has been already set
finally:
    # Defer all custom import for after initializing the forkserver and
    # ignoring the most annoying warnings
    import logging
    import random
    import sys
    from pathlib import Path
    from time import strftime
    from uuid import uuid4

    from nipype import __version__ as _nipype_ver
    from templateflow import __version__ as _tf_ver

    from . import __version__

logging.addLevelName(25, 'IMPORTANT')  # Add a new level between INFO and WARNING
logging.addLevelName(15, 'VERBOSE')  # Add a new level between INFO and DEBUG

DEFAULT_MEMORY_MIN_GB = 0.01

# Ping NiPype eTelemetry once if env var was not set
# workers on the pool will have the env variable set from the master process
if not _disable_et:
    # Just get so analytics track one hit
    from contextlib import suppress

    import requests

    with suppress((requests.ConnectionError, requests.ReadTimeout)):
        requests.get('https://rig.mit.edu/et/projects/nipy/nipype', timeout=0.05)

# Execution environment
_exec_env = os.name
_docker_ver = None
# special variable set in the container
if os.getenv('IS_DOCKER_8395080871'):
    _exec_env = 'singularity'
    _cgroup = Path('/proc/1/cgroup')
    if _cgroup.exists() and 'docker' in _cgroup.read_text():
        _docker_ver = os.getenv('DOCKER_VERSION_8395080871')
        _exec_env = 'nibabies-wrapper' if _docker_ver else 'docker'
    del _cgroup

_fs_license = os.getenv('FS_LICENSE')
if not _fs_license and os.getenv('FREESURFER_HOME'):
    _fs_home = os.getenv('FREESURFER_HOME')
    if _fs_home and (Path(_fs_home) / 'license.txt').is_file():
        _fs_license = str(Path(_fs_home) / 'license.txt')
    del _fs_home

_templateflow_home = Path(
    os.getenv('TEMPLATEFLOW_HOME', os.path.join(os.getenv('HOME'), '.cache', 'templateflow'))
)

try:
    from psutil import virtual_memory

    _free_mem_at_start = round(virtual_memory().available / 1024**3, 1)
except ImportError:
    _free_mem_at_start = None

_oc_limit = 'n/a'
_oc_policy = 'n/a'
try:
    # Memory policy may have a large effect on types of errors experienced
    _proc_oc_path = Path('/proc/sys/vm/overcommit_memory')
    if _proc_oc_path.exists():
        _oc_policy = {'0': 'heuristic', '1': 'always', '2': 'never'}.get(
            _proc_oc_path.read_text().strip(), 'unknown'
        )
        if _oc_policy != 'never':
            _proc_oc_kbytes = Path('/proc/sys/vm/overcommit_kbytes')
            if _proc_oc_kbytes.exists():
                _oc_limit = _proc_oc_kbytes.read_text().strip()
            if _oc_limit in ('0', 'n/a') and Path('/proc/sys/vm/overcommit_ratio').exists():
                _oc_limit = '{}%'.format(Path('/proc/sys/vm/overcommit_ratio').read_text().strip())
except Exception:  # noqa: S110, BLE001
    pass

# Debug modes are names that influence the exposure of internal details to
# the user, either through additional derivatives or increased verbosity
DEBUG_MODES = ('compcor', 'registration', 'fieldmaps', 'pdb')


class _Config:
    """An abstract class forbidding instantiation."""

    _paths = ()

    def __init__(self):
        """Avert instantiation."""
        raise RuntimeError('Configuration type is not instantiable.')

    @classmethod
    def load(cls, settings, init=True, ignore=None):
        """Store settings from a dictionary."""
        ignore = ignore or {}
        for k, v in settings.items():
            if k in ignore or v is None:
                continue
            if k in cls._paths:
                if isinstance(v, list | tuple):  # Multiple paths
                    setattr(cls, k, [Path(p).absolute() for p in v])
                elif isinstance(v, dict):
                    setattr(cls, k, {key: Path(val).absolute() for key, val in v.items()})
                else:
                    setattr(cls, k, Path(v).absolute())
            elif hasattr(cls, k):
                setattr(cls, k, v)

        if init:
            try:
                cls.init()
            except AttributeError:
                pass

    @classmethod
    def get(cls):
        """Return defined settings."""
        from niworkflows.utils.spaces import Reference, SpatialReferences

        out = {}
        for k, v in cls.__dict__.items():
            if k.startswith('_') or v is None:
                continue
            if callable(getattr(cls, k)):
                continue
            if k in cls._paths:
                if isinstance(v, list | tuple):  # Multiple paths
                    v = [str(p) for p in v]
                elif isinstance(v, dict):
                    v = {key: str(val) for key, val in v.items()}
                else:
                    v = str(v)
            if isinstance(v, SpatialReferences):
                v = ' '.join([str(s) for s in v.references]) or None
            if isinstance(v, Reference):
                v = str(v) or None
            out[k] = v
        return out


class environment(_Config):
    """
    Read-only options regarding the platform and environment.

    Crawls runtime descriptive settings (e.g., default FreeSurfer license,
    execution environment, nipype and *nibabies* versions, etc.).
    The ``environment`` section is not loaded in from file,
    only written out when settings are exported.
    This config section is useful when reporting issues,
    and these variables are tracked whenever the user does not
    opt-out using the ``--notrack`` argument.

    """

    cpu_count = os.cpu_count()
    """Number of available CPUs."""
    exec_docker_version = _docker_ver
    """Version of Docker Engine."""
    exec_env = _exec_env
    """A string representing the execution platform."""
    free_mem = _free_mem_at_start
    """Free memory at start."""
    overcommit_policy = _oc_policy
    """Linux's kernel virtual memory overcommit policy."""
    overcommit_limit = _oc_limit
    """Linux's kernel virtual memory overcommit limits."""
    nipype_version = _nipype_ver
    """Nipype's current version."""
    templateflow_version = _tf_ver
    """The TemplateFlow client version installed."""
    version = __version__
    """*NiBabies*'s version."""


class nipype(_Config):
    """Nipype settings."""

    crashfile_format = 'txt'
    """The file format for crashfiles, either text or pickle."""
    get_linked_libs = False
    """Run NiPype's tool to enlist linked libraries for every interface."""
    memory_gb = None
    """Estimation in GB of the RAM this workflow can allocate at any given time."""
    nprocs = os.cpu_count()
    """Number of processes (compute tasks) that can be run in parallel (multiprocessing only)."""
    omp_nthreads = None
    """Number of CPUs a single process can access for multithreaded execution."""
    plugin = 'MultiProc'
    """NiPype's execution plugin."""
    plugin_args = {
        'maxtasksperchild': 1,
        'raise_insufficient': False,
    }
    """Settings for NiPype's execution plugin."""
    resource_monitor = False
    """Enable resource monitor."""
    stop_on_first_crash = True
    """Whether the workflow should stop or continue after the first error."""

    @classmethod
    def get_plugin(cls):
        """Format a dictionary for Nipype consumption."""
        out = {
            'plugin': cls.plugin,
            'plugin_args': cls.plugin_args,
        }
        if cls.plugin in ('MultiProc', 'LegacyMultiProc'):
            out['plugin_args']['n_procs'] = int(cls.nprocs)
            if cls.memory_gb:
                out['plugin_args']['memory_gb'] = float(cls.memory_gb)
        return out

    @classmethod
    def init(cls):
        """Set NiPype configurations."""
        from nipype import config as ncfg

        # Configure resource_monitor
        if cls.resource_monitor:
            ncfg.update_config(
                {
                    'monitoring': {
                        'enabled': cls.resource_monitor,
                        'sample_frequency': '0.5',
                        'summary_append': True,
                    }
                }
            )
            ncfg.enable_resource_monitor()

        # Nipype config (logs and execution)
        ncfg.update_config(
            {
                'execution': {
                    'crashdump_dir': str(execution.log_dir),
                    'crashfile_format': cls.crashfile_format,
                    'get_linked_libs': cls.get_linked_libs,
                    'stop_on_first_crash': cls.stop_on_first_crash,
                    'check_version': False,  # disable future telemetry
                }
            }
        )

        if cls.omp_nthreads is None:
            cls.omp_nthreads = min(cls.nprocs - 1 if cls.nprocs > 1 else os.cpu_count(), 8)


class execution(_Config):
    """Configure run-level settings."""

    bids_dir = None
    """An existing path to the dataset, which must be BIDS-compliant."""
    bids_database_dir = None
    """Path to the directory containing SQLite database indices for the input BIDS dataset."""
    bids_description_hash = None
    """Checksum (SHA256) of the ``dataset_description.json`` of the BIDS dataset."""
    bids_filters = None
    """A dictionary of BIDS selection filters."""
    boilerplate_only = False
    """Only generate a boilerplate."""
    copy_derivatives = False
    """Copy any found derivatives into the output directory."""
    sloppy = False
    """Run in sloppy mode (meaning, suboptimal parameters that minimize run-time)."""
    dataset_links = {}
    """A dictionary of dataset links to be used to track Sources in sidecars."""
    debug = []
    """Debug mode(s)."""
    derivatives = {}
    """One or more paths where pre-computed derivatives are found."""
    derivatives_filters = None
    """A dictionary of BIDS selection filters"""
    echo_idx = None
    """Select a particular echo for multi-echo EPI datasets."""
    fs_license_file = _fs_license
    """An existing file containing a FreeSurfer license."""
    fs_subjects_dir = None
    """FreeSurfer's subjects directory."""
    layout = None
    """A :py:class:`~bids.layout.BIDSLayout` object, see :py:func:`init`."""
    log_dir = None
    """The path to a directory that contains execution logs."""
    log_level = 25
    """Output verbosity."""
    low_mem = None
    """Utilize uncompressed NIfTIs and other tricks to minimize memory allocation."""
    mcribs_dir = None
    """M-CRIB-S processing and output directory."""
    md_only_boilerplate = False
    """Do not convert boilerplate from MarkDown to LaTex and HTML."""
    nibabies_dir = None
    """Root of NiBabies BIDS Derivatives dataset. Depends on output_layout."""
    notrack = False
    """Do not monitor *nibabies* using *migas*."""
    output_dir = None
    """Folder where derivatives will be stored."""
    me_output_echos = False
    """Output individual echo time series with slice, motion and susceptibility correction"""
    output_layout = None
    """Layout of derivatives within output_dir."""
    output_spaces = None
    """List of (non)standard spaces designated (with the ``--output-spaces`` flag of
    the command line) as spatial references for outputs."""
    parameters_hash = None
    """Unique hash of the current configuration parameters."""
    reference_anat = None
    """Force usage of this anatomical scan as the structural reference."""
    reports_only = False
    """Only build the reports, based on the reportlets found in a cached working directory."""
    run_uuid = f'{strftime("%Y%m%d-%H%M%S")}_{uuid4()}'
    """Unique identifier of this particular run."""
    segmentation_atlases_dir = None
    """Directory with atlases to use for JLF segmentations"""
    session_id = None
    """List of session identifiers that are to be preprocessed."""
    participant_label = None
    """List of participant identifiers that are to be preprocessed."""
    task_id = None
    """Select a particular task from all available in the dataset."""
    templateflow_home = _templateflow_home
    """The root folder of the TemplateFlow client."""
    unique_labels = None
    """Combinations of subject + session identifiers to be preprocessed."""
    work_dir = Path('work').absolute()
    """Path to a working directory where intermediate results will be available."""
    write_graph = False
    """Write out the computational graph corresponding to the planned preprocessing."""
    _layout = None

    _paths = (
        'anat_derivatives',
        'bids_dir',
        'bids_database_dir',
        'dataset_links',
        'derivatives',
        'fs_license_file',
        'fs_subjects_dir',
        'layout',
        'log_dir',
        'mcribs_dir',
        'nibabies_dir',
        'output_dir',
        'segmentation_atlases_dir',
        'templateflow_home',
        'work_dir',
    )

    @classmethod
    def init(cls):
        """Create a new BIDS Layout accessible with :attr:`~execution.layout`."""
        # Convert string literal None to NoneType
        if cls.unique_labels:
            cls.unique_labels = [
                [sub, ses] if ses != 'None' else [sub, None] for sub, ses in cls.unique_labels
            ]

        if cls.fs_license_file and Path(cls.fs_license_file).is_file():
            os.environ['FS_LICENSE'] = str(cls.fs_license_file)

        if cls._layout is None:
            import re

            from bids.layout import BIDSLayout
            from bids.layout.index import BIDSLayoutIndexer

            _db_path = cls.bids_database_dir or (cls.work_dir / cls.run_uuid / 'bids_db')
            _db_path.mkdir(exist_ok=True, parents=True)

            # Recommended after PyBIDS 12.1
            ignore_patterns = [
                'code',
                'stimuli',
                'sourcedata',
                'models',
                re.compile(r'^\.'),
                re.compile(r'sub-[a-zA-Z0-9]+(/ses-[a-zA-Z0-9]+)?/(beh|dwi|eeg|ieeg|meg|perf)'),
            ]
            if cls.participant_label and cls.bids_database_dir is None:
                # Ignore any subjects who aren't the requested ones.
                # This is only done if the database is written out to a run-specific folder.
                ignore_patterns.append(
                    re.compile(r'sub-(?!(' + '|'.join(cls.participant_label) + r')(\b|_))')
                )

            _indexer = BIDSLayoutIndexer(
                validate=False,
                ignore=ignore_patterns,
            )
            cls._layout = BIDSLayout(
                str(cls.bids_dir),
                database_path=_db_path,
                reset_database=cls.bids_database_dir is None,
                indexer=_indexer,
            )
            cls.bids_database_dir = _db_path
        cls.layout = cls._layout
        if cls.bids_filters:
            from bids.layout import Query

            def _process_value(value):
                """Convert string with "Query" in it to Query object."""
                if isinstance(value, list):
                    return [_process_value(val) for val in value]
                else:
                    return (
                        getattr(Query, value[7:-4])
                        if not isinstance(value, Query) and 'Query' in value
                        else value
                    )

            # unserialize pybids Query enum values
            for acq, filters in cls.bids_filters.items():
                for k, v in filters.items():
                    cls.bids_filters[acq][k] = _process_value(v)

        dataset_links = {
            'raw': cls.bids_dir,
            'templateflow': Path(TF_LAYOUT.root),
        }
        for deriv_name, deriv_path in cls.derivatives.items():
            dataset_links[deriv_name] = deriv_path
        cls.dataset_links = dataset_links

        if 'all' in cls.debug:
            cls.debug = list(DEBUG_MODES)


# These variables are not necessary anymore
del _fs_license
del _exec_env
del _templateflow_home
del _free_mem_at_start
del _oc_limit
del _oc_policy


class workflow(_Config):
    """Configure the particular execution graph of this workflow."""

    age_months = None
    """Age (in months)"""
    analysis_level = 'participant'
    """Level of analysis."""
    anat_only = False
    """Execute the anatomical preprocessing only."""
    bold2anat_dof = None
    """Degrees of freedom of the BOLD-to-anatomical registration steps."""
    bold2anat_init = 'auto'
    """Whether to use standard coregistration ('register') or to initialize coregistration from the
    BOLD image-header ('header')."""
    cifti_output = None
    """Generate HCP Grayordinates, accepts either ``'91k'`` (default) or ``'170k'``."""
    dummy_scans = None
    """Set a number of initial scans to be considered nonsteady states."""
    fd_radius = 45
    """Head radius in mm for framewise displacement calculation"""
    fmap_bspline = None
    """Regularize fieldmaps with a field of B-Spline basis."""
    fmap_demean = None
    """Remove the mean from fieldmaps."""
    force_syn = None
    """Run *fieldmap-less* susceptibility-derived distortions estimation."""
    hires = None
    """Run FreeSurfer ``recon-all`` with the ``-hires`` flag."""
    hmc_bold_frame = 16
    """Frame to start head motion correction estimation on BOLD"""
    ignore = None
    """Ignore particular steps for *nibabies*."""
    level = None
    """Level of preprocessing to complete. One of ['minimal', 'resampling', 'full']."""
    longitudinal = False
    """Run FreeSurfer ``recon-all`` with the ``-logitudinal`` flag."""
    medial_surface_nan = None
    """Fill medial surface with :abbr:`NaNs (not-a-number)` when sampling."""
    multi_step_reg = True
    """Perform multiple registrations (native -> MNIInfant -> template) and concatenate into a
    single transform"""
    norm_csf = False
    """Replace low intensity voxels in CSF mask with average."""
    project_goodvoxels = False
    """Exclude voxels with locally high coefficient of variation from sampling."""
    regressors_all_comps = None
    """Return all CompCor components."""
    regressors_dvars_th = None
    """Threshold for DVARS."""
    regressors_fd_th = None
    """Threshold for :abbr:`FD (frame-wise displacement)`."""
    run_reconall = True
    """Run FreeSurfer's surface reconstruction."""
    skull_strip_fixed_seed = False
    """Fix a seed for skull-stripping."""
    skull_strip_template = 'UNCInfant:cohort-1'
    """Change default brain extraction template."""
    skull_strip_anat = 'force'
    """Skip brain extraction of the anatomical images (default is ``force``, meaning that
    *nibabies* will run brain extraction of the T1w)."""
    slice_time_ref = 0.5
    """The time of the reference slice to correct BOLD values to, as a fraction
    acquisition time. 0 indicates the start, 0.5 the midpoint, and 1 the end
    of acquisition. The alias `start` corresponds to 0, and `middle` to 0.5.
    The default value is 0.5."""
    spaces = None
    """Keeps the :py:class:`~niworkflows.utils.spaces.SpatialReferences`
    instance keeping standard and nonstandard spaces."""
    surface_recon_method = None
    """Method to use for surface reconstruction."""
    use_bbr = False
    """Run boundary-based registration for BOLD-to-T1w registration."""
    use_syn_sdc = None
    """Run *fieldmap-less* susceptibility-derived distortions estimation
    in the absence of any alternatives."""
    me_t2s_fit_method = 'curvefit'
    """The method by which to estimate T2*/S0 for multi-echo data"""


class loggers:
    """Keep loggers easily accessible (see :py:func:`init`)."""

    _fmt = '%(asctime)s,%(msecs)d %(name)-2s %(levelname)-2s:\n\t %(message)s'
    _datefmt = '%y%m%d-%H:%M:%S'

    default = logging.getLogger()
    """The root logger."""
    cli = logging.getLogger('cli')
    """Command-line interface logging."""
    workflow = logging.getLogger('nipype.workflow')
    """NiPype's workflow logger."""
    interface = logging.getLogger('nipype.interface')
    """NiPype's interface logger."""
    utils = logging.getLogger('nipype.utils')
    """NiPype's utils logger."""

    @classmethod
    def init(cls):
        """
        Set the log level, initialize all loggers into :py:class:`loggers`.

            * Add new logger levels (25: IMPORTANT, and 15: VERBOSE).
            * Add a new sub-logger (``cli``).
            * Logger configuration.

        """
        from nipype import config as ncfg

        _handler = logging.StreamHandler(stream=sys.stdout)
        _handler.setFormatter(logging.Formatter(fmt=cls._fmt, datefmt=cls._datefmt))
        cls.cli.addHandler(_handler)
        cls.default.setLevel(execution.log_level)
        cls.cli.setLevel(execution.log_level)
        cls.interface.setLevel(execution.log_level)
        cls.workflow.setLevel(execution.log_level)
        cls.utils.setLevel(execution.log_level)
        ncfg.update_config(
            {'logging': {'log_directory': str(execution.log_dir), 'log_to_file': True}}
        )


class seeds(_Config):
    """Initialize the PRNG and track random seed assignments"""

    _random_seed = None
    master = None
    """Master random seed to initialize the Pseudorandom Number Generator (PRNG)"""
    ants = None
    """Seed used for antsRegistration, antsAI, antsMotionCorr"""
    numpy = None
    """Seed used by NumPy"""

    @classmethod
    def init(cls):
        if cls._random_seed is not None:
            cls.master = cls._random_seed
        if cls.master is None:
            cls.master = random.randint(1, 65536)
        random.seed(cls.master)  # initialize the PRNG
        # functions to set program specific seeds
        cls.ants = _set_ants_seed()
        cls.numpy = _set_numpy_seed()


def _set_ants_seed():
    """Fix random seed for antsRegistration, antsAI, antsMotionCorr"""
    val = random.randint(1, 65536)
    os.environ['ANTS_RANDOM_SEED'] = str(val)
    return val


def _set_numpy_seed():
    """NumPy's random seed is independent from Python's `random` module"""
    import numpy as np

    val = random.randint(1, 65536)
    np.random.seed(val)
    return val


def from_dict(settings, init=True, ignore=None):
    """Read settings from a flat dictionary.

    Arguments
    ---------
    setting : dict
        Settings to apply to any configuration
    init : `bool` or :py:class:`~collections.abc.Container`
        Initialize all, none, or a subset of configurations.
    ignore : :py:class:`~collections.abc.Container`
        Collection of keys in ``setting`` to ignore
    """

    # Accept global True/False or container of configs to initialize
    def initialize(x):
        return init if init in (True, False) else x in init

    nipype.load(settings, init=initialize('nipype'), ignore=ignore)
    execution.load(settings, init=initialize('execution'), ignore=ignore)
    workflow.load(settings, init=initialize('workflow'), ignore=ignore)
    seeds.load(settings, init=initialize('seeds'), ignore=ignore)
    loggers.init()


def load(filename, skip=None, init=True):
    """Load settings from file.

    Arguments
    ---------
    filename : :py:class:`os.PathLike`
        TOML file containing fMRIPrep configuration.
    skip : dict or None
        Sets of values to ignore during load, keyed by section name
    init : `bool` or :py:class:`~collections.abc.Container`
        Initialize all, none, or a subset of configurations.
    """
    from toml import loads

    skip = skip or {}

    # Accept global True/False or container of configs to initialize
    def initialize(x):
        return init if init in (True, False) else x in init

    filename = Path(filename)
    settings = loads(filename.read_text())
    for sectionname, configs in settings.items():
        if sectionname != 'environment':
            section = getattr(sys.modules[__name__], sectionname)
            ignore = skip.get(sectionname)
            section.load(configs, ignore=ignore, init=initialize(sectionname))


def get(flat=False):
    """Get config as a dict."""
    settings = {
        'environment': environment.get(),
        'execution': execution.get(),
        'workflow': workflow.get(),
        'nipype': nipype.get(),
        'seeds': seeds.get(),
    }
    if not flat:
        return settings

    return {
        '.'.join((section, k)): v
        for section, configs in settings.items()
        for k, v in configs.items()
    }


def dumps():
    """Format config into toml."""
    from toml import dumps

    return dumps(get())


def to_filename(filename):
    """Write settings to file."""
    filename = Path(filename)
    filename.write_text(dumps())


def _process_initializer(cwd, omp_nthreads):
    """Initialize the environment of the child process."""
    os.chdir(cwd)
    os.environ['NIPYPE_NO_ET'] = '1'
    os.environ['OMP_NUM_THREADS'] = f'{omp_nthreads}'


def dismiss_entities(entities: list | None = None) -> list:
    """Set entities to dismiss in a DerivativesDataSink."""
    from niworkflows.utils.connections import listify

    entities = set(entities or [])
    echo_idx = execution.echo_idx
    if echo_idx is None or len(listify(echo_idx)) > 2:
        entities.add('echo')
    output_layout = execution.output_layout
    if output_layout != 'multiverse':
        entities.add('hash')
    return list(entities)


DEFAULT_DISMISS_ENTITIES = dismiss_entities()

DEFAULT_CONFIG_HASH_FIELDS = {
    'execution': [
        'sloppy',
        'echo_idx',
        'reference_anat',
    ],
    'workflow': [
        'surface_recon_method',
        'bold2anat_dof',
        'bold2anat_init',
        'dummy_scans',
        'fd_radius',
        'fmap_bspline',
        'fmap_demean',
        'force_syn',
        'hmc_bold_frame',
        'longitudinal',
        'medial_surface_nan',
        'multi_step_reg',
        'norm_csf',
        'project_goodvoxels',
        'regressors_dvars_th',
        'regressors_fd_th',
        'skull_strip_fixed_seed',
        'skull_strip_template',
        'skull_strip_anat',
        'slice_time_ref',
        'surface_recon_method',
        'use_bbr',
        'use_syn_sdc',
        'me_t2s_fit_method',
    ],
}


def hash_config(
    conf: dict[str, ty.Any],
    *,
    fields_required: dict[str, list[str]] = DEFAULT_CONFIG_HASH_FIELDS,
    version: str = None,
    digest_size: int = 4,
) -> str:
    """
    Generate a unique BLAKE2b hash of configuration attributes.

    By default, uses a preselected list of workflow-altering parameters.
    """
    import json
    from hashlib import blake2b

    if version is None:
        from nibabies import __version__ as version

    data = {}
    for level, fields in fields_required.items():
        for f in fields:
            data[f] = conf[level].get(f, None)

    datab = json.dumps(data, sort_keys=True).encode()
    return blake2b(datab, digest_size=digest_size).hexdigest()
