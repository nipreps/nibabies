## *NiBabies*: A robust preprocessing workflow tailored for neonate and infant MRI

![nibabies](https://github.com/nipreps/nibabies/actions/workflows/pytest.yml/badge.svg)
[![DOI](https://zenodo.org/badge/264223087.svg)](https://zenodo.org/badge/latestdoi/264223087)


<p align="center">
    <img src="./docs/_static/nibabies_anat.png" width="600" height="600" alt="nibabies-anat">
</p>

*NiBabies* is an extension of [fMRIPrep](https://fmriprep.org/en/stable/) designed and tested for infants 0-2 years old. *NiBabies* offers structural and functional MRI preprocessing.

----
### Requirements

Given its extensive dependencies, the easiest way to get up and running with *NiBabies* is by using the available [Docker](https://hub.docker.com/r/mgxd/nibabies/tags?page=1&ordering=last_updated) or [Singularity](https://cloud.sylabs.io/library/mathiasg/default/nibabies) containers.

If you insist on installing this tool locally, you can use the [Dockerfile](./Dockerfile) as a guide.


### Usage

*NiBabies* follow the [BIDS App Specifications](http://bids-apps.neuroimaging.io/about/), meaning you only need to provide three positional arguments:

- **bids_dir** - the root folder of a BIDS valid dataset.
- **output_dir** - folder to store outputs and reports.
- **level** - processing stage to be run, currently can only be `participant`.


However, as infant brains can vastly differ depending on age, providing the following arguments is highly recommended:

- **--age-months** - participant age in months
- **--segmentation-dir** - directory containing pre-labeled segmentations to use for Joint Label Fusion.

<details>
<summary>Extensive argument list</summary>

```
usage: nibabies [-h] [--version] [--skip_bids_validation]
                [--participant-label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]]
                [-t TASK_ID] [--echo-idx ECHO_IDX] [--bids-filter-file FILE]
                [--anat-derivatives PATH] [--bids-database-dir PATH]
                [--nprocs NPROCS] [--omp-nthreads OMP_NTHREADS]
                [--mem MEMORY_GB] [--low-mem] [--use-plugin FILE]
                [--anat-only] [--boilerplate_only] [--md-only-boilerplate]
                [--error-on-aroma-warnings] [-v]
                [--ignore {fieldmaps,slicetiming,sbref,t2w,flair} [{fieldmaps,slicetiming,sbref,t2w,flair} ...]]
                [--longitudinal]
                [--output-spaces [OUTPUT_SPACES [OUTPUT_SPACES ...]]]
                [--bold2t1w-init {register,header}] [--bold2t1w-dof {6,9,12}]
                [--force-bbr] [--force-no-bbr] [--medial-surface-nan]
                [--dummy-scans DUMMY_SCANS] [--random-seed _RANDOM_SEED]
                [--use-aroma]
                [--aroma-melodic-dimensionality AROMA_MELODIC_DIM]
                [--return-all-components]
                [--fd-spike-threshold REGRESSORS_FD_TH]
                [--dvars-spike-threshold REGRESSORS_DVARS_TH]
                [--skull-strip-template SKULL_STRIP_TEMPLATE]
                [--skull-strip-fixed-seed]
                [--skull-strip-t1w {auto,skip,force}] [--fmap-bspline]
                [--fmap-no-demean] [--use-syn-sdc] [--force-syn]
                [--fs-license-file FILE] [--fs-subjects-dir PATH]
                [--no-submm-recon] [--cifti-output [{91k,170k}] |
                --fs-no-reconall] [--output-layout {bids,legacy}]
                [-w WORK_DIR] [--clean-workdir] [--resource-monitor]
                [--reports-only] [--config-file FILE] [--write-graph]
                [--stop-on-first-crash] [--notrack]
                [--debug {compcor,registration,fieldmaps,all} [{compcor,registration,fieldmaps,all} ...]]
                [--sloppy]
                [--age-months AGE_MONTHS]
                [--segmentation-atlases-dir SEGMENTATION_ATLASES_DIR]
                [--ants-affine-init {random,search}]
                bids_dir output_dir {participant}

NiBabies: Preprocessing workflows for infants (version 0.1.0)

positional arguments:
  bids_dir              the root folder of a BIDS valid dataset (sub-XXXXX
                        folders should be found at the top level in this
                        folder).
  output_dir            the output path for the outcomes of preprocessing and
                        visual reports
  {participant}         processing stage to be run, only "participant" in the
                        case of NiBabies (see BIDS-Apps specification).

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit

Options for filtering BIDS queries:
  --skip_bids_validation, --skip-bids-validation
                        assume the input dataset is BIDS compliant and skip
                        the validation (default: False)
  --participant-label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...], --participant_label PARTICIPANT_LABEL [PARTICIPANT_LABEL ...]
                        a space delimited list of participant identifiers or a
                        single identifier (the sub- prefix can be removed)
                        (default: None)
  -t TASK_ID, --task-id TASK_ID
                        select a specific task to be processed (default: None)
  --echo-idx ECHO_IDX   select a specific echo to be processed in a multiecho
                        series (default: None)
  --bids-filter-file FILE
                        a JSON file describing custom BIDS input filters using
                        PyBIDS. For further details, please check out
                        https://fmriprep.readthedocs.io/en/0.0.3/faq.html#how-
                        do-I-select-only-certain-files-to-be-input-to-fMRIPrep
                        (default: None)
  --anat-derivatives PATH
                        Reuse the anatomical derivatives from another NiBabies
                        run or calculated with an alternative processing tool
                        (NOT RECOMMENDED). (default: None)
  --bids-database-dir PATH
                        Path to an existing PyBIDS database folder, for faster
                        indexing (especially useful for large datasets).
                        (default: None)

Options to handle performance:
  --nprocs NPROCS, --nthreads NPROCS, --n_cpus NPROCS, --n-cpus NPROCS
                        maximum number of threads across all processes
                        (default: None)
  --omp-nthreads OMP_NTHREADS
                        maximum number of threads per-process (default: None)
  --mem MEMORY_GB, --mem_mb MEMORY_GB, --mem-mb MEMORY_GB
                        upper bound memory limit for NiBabies processes
                        (default: None)
  --low-mem             attempt to reduce memory usage (will increase disk
                        usage in working directory) (default: False)
  --use-plugin FILE, --nipype-plugin-file FILE
                        nipype plugin configuration file (default: None)
  --anat-only           run anatomical workflows only (default: False)
  --boilerplate_only    generate boilerplate only (default: False)
  --md-only-boilerplate
                        skip generation of HTML and LaTeX formatted citation
                        with pandoc (default: False)
  --error-on-aroma-warnings
                        Raise an error if ICA_AROMA does not produce sensible
                        output (e.g., if all the components are classified as
                        signal or noise) (default: False)
  -v, --verbose         increases log verbosity for each occurence, debug
                        level is -vvv (default: 0)

Workflow configuration:
  --ignore {fieldmaps,slicetiming,sbref,t2w,flair} [{fieldmaps,slicetiming,sbref,t2w,flair} ...]
                        ignore selected aspects of the input dataset to
                        disable corresponding parts of the workflow (a space
                        delimited list) (default: [])
  --longitudinal        treat dataset as longitudinal - may increase runtime
                        (default: False)
  --output-spaces [OUTPUT_SPACES [OUTPUT_SPACES ...]]
                        Standard and non-standard spaces to resample
                        anatomical and functional images to. Standard spaces
                        may be specified by the form
                        ``<SPACE>[:cohort-<label>][:res-<resolution>][...]``,
                        where ``<SPACE>`` is a keyword designating a spatial
                        reference, and may be followed by optional, colon-
                        separated parameters. Non-standard spaces imply
                        specific orientations and sampling grids. Important to
                        note, the ``res-*`` modifier does not define the
                        resolution used for the spatial normalization. To
                        generate no BOLD outputs, use this option without
                        specifying any spatial references. For further
                        details, please check out
                        https://fmriprep.readthedocs.io/en/0.0.3/spaces.html
                        (default: None)
  --bold2t1w-init {register,header}
                        Either "register" (the default) to initialize volumes
                        at center or "header" to use the header information
                        when coregistering BOLD to T1w images. (default:
                        register)
  --bold2t1w-dof {6,9,12}
                        Degrees of freedom when registering BOLD to T1w
                        images. 6 degrees (rotation and translation) are used
                        by default. (default: 6)
  --force-bbr           Always use boundary-based registration (no goodness-
                        of-fit checks) (default: None)
  --force-no-bbr        Do not use boundary-based registration (no goodness-
                        of-fit checks) (default: None)
  --medial-surface-nan  Replace medial wall values with NaNs on functional
                        GIFTI files. Only performed for GIFTI files mapped to
                        a freesurfer subject (fsaverage or fsnative).
                        (default: False)
  --dummy-scans DUMMY_SCANS
                        Number of non steady state volumes. (default: None)
  --random-seed _RANDOM_SEED
                        Initialize the random seed for the workflow (default:
                        None)

Specific options for running ICA_AROMA:
  --use-aroma           add ICA_AROMA to your preprocessing stream (default:
                        False)
  --aroma-melodic-dimensionality AROMA_MELODIC_DIM
                        Exact or maximum number of MELODIC components to
                        estimate (positive = exact, negative = maximum)
                        (default: -200)

Specific options for estimating confounds:
  --return-all-components
                        Include all components estimated in CompCor
                        decomposition in the confounds file instead of only
                        the components sufficient to explain 50 percent of
                        BOLD variance in each CompCor mask (default: False)
  --fd-spike-threshold REGRESSORS_FD_TH
                        Threshold for flagging a frame as an outlier on the
                        basis of framewise displacement (default: 0.5)
  --dvars-spike-threshold REGRESSORS_DVARS_TH
                        Threshold for flagging a frame as an outlier on the
                        basis of standardised DVARS (default: 1.5)

Specific options for ANTs registrations:
  --skull-strip-template SKULL_STRIP_TEMPLATE
                        select a template for skull-stripping with
                        antsBrainExtraction (default: UNCInfant:cohort-1)
  --skull-strip-fixed-seed
                        do not use a random seed for skull-stripping - will
                        ensure run-to-run replicability when used with --omp-
                        nthreads 1 and matching --random-seed <int> (default:
                        False)
  --skull-strip-t1w {auto,skip,force}
                        determiner for T1-weighted skull stripping ('force'
                        ensures skull stripping, 'skip' ignores skull
                        stripping, and 'auto' applies brain extraction based
                        on the outcome of a heuristic to check whether the
                        brain is already masked). (default: force)

Specific options for handling fieldmaps:
  --fmap-bspline        fit a B-Spline field using least-squares
                        (experimental) (default: False)
  --fmap-no-demean      do not remove median (within mask) from fieldmap
                        (default: True)

Specific options for SyN distortion correction:
  --use-syn-sdc         EXPERIMENTAL: Use fieldmap-free distortion correction
                        (default: False)
  --force-syn           EXPERIMENTAL/TEMPORARY: Use SyN correction in addition
                        to fieldmap correction, if available (default: False)

Specific options for FreeSurfer preprocessing:
  --fs-license-file FILE
                        Path to FreeSurfer license key file. Get it (for free)
                        by registering at
                        https://surfer.nmr.mgh.harvard.edu/registration.html
                        (default: None)
  --fs-subjects-dir PATH
                        Path to existing FreeSurfer subjects directory to
                        reuse. (default: OUTPUT_DIR/freesurfer) (default:
                        None)

Surface preprocessing options:
  --no-submm-recon      disable sub-millimeter (hires) reconstruction
                        (default: True)
  --cifti-output [{91k,170k}]
                        output preprocessed BOLD as a CIFTI dense timeseries.
                        Optionally, the number of grayordinate can be
                        specified (default is 91k, which equates to 2mm
                        resolution) (default: False)
  --fs-no-reconall      disable FreeSurfer surface preprocessing. (default:
                        True)

Other options:
  --output-layout {bids,legacy}
                        Organization of outputs. legacy (default) creates
                        derivative datasets as subdirectories of outputs. bids
                        places NiBabies derivatives directly in the output
                        directory, and defaults to placing FreeSurfer
                        derivatives in <output-dir>/sourcedata/freesurfer.
                        (default: legacy)
  -w WORK_DIR, --work-dir WORK_DIR
                        path where intermediate results should be stored
                        (default: /tmp/work)
  --clean-workdir       Clears working directory of contents. Use of this flag
                        is notrecommended when running concurrent processes of
                        NiBabies. (default: False)
  --resource-monitor    enable Nipype's resource monitoring to keep track of
                        memory and CPU usage (default: False)
  --reports-only        only generate reports, don't run workflows. This will
                        only rerun report aggregation, not reportlet
                        generation for specific nodes. (default: False)
  --config-file FILE    Use pre-generated configuration file. Values in file
                        will be overridden by command-line arguments.
                        (default: None)
  --write-graph         Write workflow graph. (default: False)
  --stop-on-first-crash
                        Force stopping on first crash, even if a work
                        directory was specified. (default: False)
  --notrack             Opt-out of sending tracking information of this run to
                        the NiBabies developers. This information helps to
                        improve NiBabies and provides an indicator of real
                        world usage crucial for obtaining funding. (default:
                        False)
  --debug {compcor,registration,fieldmaps,all} [{compcor,registration,fieldmaps,all} ...]
                        Debug mode(s) to enable. 'all' is alias for all
                        available modes. (default: None)
  --sloppy              Use low-quality tools for speed - TESTING ONLY
                        (default: False)

NiBabies specific options:
  --age-months AGE_MONTHS
                        Age in months (default: None)
  --segmentation-atlases-dir SEGMENTATION_ATLASES_DIR
                        Directory containing precalculated segmentations to
                        use for JointLabelFusion. (default: None)
```

</details>

### Outputs

TODO - Refer to [fMRIPrep's outputs](https://fmriprep.org/en/20.2.1/outputs.html) for now.
