# Outputs of *NiBabies*

*NiBabies* outputs conform to the {abbr}`BIDS (brain imaging data structure)`
Derivatives specification (see `BIDS Derivatives`_, along with the
upcoming `BEP 011`_ and `BEP 012`_).
*NiBabies* generates three broad classes of outcomes:

1. **Visual QA (quality assessment) reports**:
   one {abbr}`HTML (hypertext markup language)` per subject,
   that allows the user a thorough visual assessment of the quality
   of processing and ensures the transparency of *NiBabies* operations.

2. **Derivatives (preprocessed data)** the input fMRI data ready for
   analysis, i.e., after the various preparation procedures
   have been applied.
   For example, {abbr}`INU (intensity non-uniformity)`-corrected versions
   of the T1-weighted image (per subject), the brain mask,
   or {abbr}`BOLD (blood-oxygen level dependent)`
   images after head-motion correction, slice-timing correction and aligned into
   the same-subject's anatomical space or in some standard space.

3. **Confounds**: this is a special family of derivatives that can be utilized
   to inform subsequent denoising steps.

:::{warning}
Confounds, particularly a/t CompCor, are still in alpha and require additional testing.
:::

:::{important}
In order to remain agnostic to any possible subsequent analysis,
*NiBabies* does not perform any denoising (e.g., spatial smoothing) itself.
There are exceptions to this principle (described in its corresponding
section below):

- CompCor regressors, which are calculated after temporal high-pass filtering.
:::

## Layout

Assuming NiBabies is invoked with:

```
nibabies <input_dir>/ <output_dir>/ participant [OPTIONS]
```

The outputs will be a [BIDS Derivatives](https://bids-specification.readthedocs.io/en/stable/derivatives/introduction.html) dataset of the form:

```
<output_dir>/
├── logs/
├── sub-<label>/
├── sub-<label>[_ses-<label>].html
├── dataset_description.json
└── .bidsignore
```

For each participant in the dataset,
a directory of derivatives (`sub-<label>/`)
and a visual report (`sub-<label>.html`) are generated.
The log directory contains `citation boilerplate`_ text.
`dataset_description.json` is a metadata file in which NiBabies
records metadata recommended by the BIDS standard.

:::{admonition} fMRIPrep deviation - multi-session processing
:class: important

Due to the potential for significant development across sessions,
NiBabies processes at the session level. Each unique session will be nested
within a subject's output directory, and an individual report will be created
(i.e. `sub-<label>_ses-<label>.html`).
:::

## Processing level

As of version 24.0.0, NiBabies supports three levels of derivatives:

* `--level minimal`: This processing mode aims to produce the smallest
  working directory and output dataset possible, while enabling all further
  processing results to be deterministically generated. Most components of
  the `visual reports`_ can be generated at this level, so the quality of
  preprocessing can be assessed. Because no resampling is done, confounds
  and carpetplots will be missing from the reports.
* `--level resampling`: This processing mode aims to produce additional
  derivatives that enable third-party resampling, resampling BOLD series
  in the working directory as needed, but these are not saved to the output
  directory.
  The `--me-output-echos` flag will be enabled at this level, in which
  case the individual echos will be saved to the working directory after
  slice-timing correction, head-motion correction, and susceptibility
  distortion correction.
* `--level full`: This processing mode aims to produce all derivatives
  that have previously been a part of the NiBabies output dataset.
  This is the default processing level.

## Visual Reports

*NiBabies* outputs summary reports, written to `<output dir>/nibabies/sub-<subject_label>[_ses-<session_label>].html`.
These reports provide a quick way to make visual inspection of the results easy.
`View a sample report. <_static/SampleReport/sample_report.html>`_

## Derivatives of *NiBabies* (preprocessed data)

Preprocessed, or derivative, data are written to
`<output dir>/sub-<subject_label>/`. If the data is composed of sessions, each session will
have a separate output.
The `BIDS Derivatives`_ specification describes the naming and metadata conventions we follow.

### Anatomical derivatives

Anatomical derivatives are placed in each subject's `anat` subfolder::

```
sub-<subject_label>/[ses-<session_label>/]
  anat/
    sub-<subject_label>[_space-<space_label>]_desc-preproc_T1w.nii.gz
    sub-<subject_label>[_space-<space_label>]_desc-preproc_T2w.nii.gz
    sub-<subject_label>[_space-<space_label>]_desc-brain_mask.nii.gz
    sub-<subject_label>[_space-<space_label>]_dseg.nii.gz
    sub-<subject_label>[_space-<space_label>]_label-CSF_probseg.nii.gz
    sub-<subject_label>[_space-<space_label>]_label-GM_probseg.nii.gz
    sub-<subject_label>[_space-<space_label>]_label-WM_probseg.nii.gz
```

Spatially-standardized derivatives are denoted with a space label,
such as `MNI152NLin2009cAsym`.

:::{versionchanged} 24.0.0
The anatomical output reference can now be either in T1w or T2w space, and will always
be denoted with the `space-{T1w/T2w}` entity.
:::

Additionally, the following transforms are saved

```
sub-<subject_label>/[ses-<session_label>/]
  anat/
    sub-<subject_label>_from-MNIInfant+<cohort>_to-<anat>_mode-image_xfm.h5
    sub-<subject_label>_from-<anat>_to-MNIInfant+<cohort>_mode-image_xfm.h5
```

If surface reconstructions are used, the following surface files are generated

```
sub-<subject_label>/[ses-<session_label>/]
  anat/
    sub-<subject_label>_hemi-[LR]_white.surf.gii
    sub-<subject_label>_hemi-[LR]_midthickness.surf.gii
    sub-<subject_label>_hemi-[LR]_pial.surf.gii
    sub-<subject_label>_hemi-[LR]_desc-reg_sphere.surf.gii
    sub-<subject_label>_hemi-[LR]_space-fsLR_desc-reg_sphere.surf.gii
```

The registration spheres target `fsaverage` and `fsLR` spaces.

:::{warning}
Unlike fMRIPrep, MSMSulc support is not available at the moment.
:::

And the affine translation (and inverse) between the anatomical reference sampling and
FreeSurfer's conformed space for surface reconstruction (`fsnative`) is stored in::

```
sub-<subject_label>/[ses-<session_label>/]
  anat/
    sub-<subject_label>_from-fsnative_to-anat_mode-image_xfm.txt
    sub-<subject_label>_from-anat_to-fsnative_mode-image_xfm.txt
```

Finally, cortical thickness, curvature, and sulcal depth maps are converted to GIFTI
and CIFTI-2::

```
sub-<subject_label>/[ses-<session_label>/]
  anat/
    sub-<subject_label>_hemi-[LR]_thickness.shape.gii
    sub-<subject_label>_hemi-[LR]_curv.shape.gii
    sub-<subject_label>_hemi-[LR]_sulc.shape.gii
    sub-<subject_label>_space-fsLR_den-32k_thickness.dscalar.nii
    sub-<subject_label>_space-fsLR_den-32k_curv.dscalar.nii
    sub-<subject_label>_space-fsLR_den-32k_sulc.dscalar.nii
```

:::{warning}

GIFTI metric files follow the FreeSurfer conventions and are not modified
by *NiBabies* in any way.

The Human Connectome Project (HCP) inverts the sign of the curvature and
sulcal depth maps. For consistency with HCP, *NiBabies* follows these
conventions and masks the medial wall of CIFTI-2 dscalar files.
:::

TODO: FSDERIVS? .. _fsderivs:

### Surface Reconstruction

If any of the surface reconstruction methods (M-CRIB-S, infant FreeSurfer, or FreeSurfer) are enabled,
then a FreeSurfer-like subjects directory is created in
`<output dir>/sourcedata/freesurfer` or the directory indicated with the
`--fs-subjects-dir` flag.

Additionally, FreeSurfer segmentations are resampled into the BOLD space,
and lookup tables are provided.

```
<output_dir>/
  sourcedata/
    freesurfer/
      fsaverage{,5,6}/
          mri/
          surf/
          ...
      sub-<label>/
          mri/
          surf/
          ...
      ...
  desc-aparc_dseg.tsv
  desc-aparcaseg_dseg.tsv
```

Copies of the `fsaverage` subjects distributed with the running version of
FreeSurfer are copied into this subjects directory, if any functional data are
sampled to those subject spaces.

Note that the use of `sourcedata/` recognizes FreeSurfer derivatives as an input to
the NiBabies workflow.
This is strictly true when pre-computed FreeSurfer derivatives are provided either in
the `sourcedata/` directory or passed via the `--fs-subjects-dir` flag;
if NiBabies runs FreeSurfer, then there is a mutual dependency.

 **If M-CRIB-S is used** for surface reconstruction, an additional folder is created in `<output dir>/sourcedata/mcribs`. M-CRIB-S converts and organizes its output into a FreeSurfer-compatible format, i.e. following the structure of FreeSurfer’s [recon-all](https://surfer.nmr.mgh.harvard.edu/fswiki/ReconAllOutputFiles), to `mcribs/sub-<label>_ses-<label>/freesurfer/`. This is then mapped by Infant fMRIPrep into the `sourcedata/freesurfer/` folder for further processing.

### Functional derivatives

Functional derivatives are stored in the `func/` subfolder.
All derivatives contain `task-<task_label>` (mandatory) and `run-<run_index>` (optional), and
these will be indicated with `[specifiers]`

```
  sub-<subject_label>/[ses-<session_label>/]
    func/
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-brain_mask.nii.gz
      sub-<subject_label>_[specifiers]_space-<space_label>_desc-preproc_bold.nii.gz
```

:::{note}
The mask file is part of the *minimal* processing level. The BOLD series
is only generated at the *full* processing level.
:::

#### Motion correction outputs

Head-motion correction (HMC) produces a reference image to which all volumes
are aligned, and a corresponding transform that maps the original BOLD series
to the reference image::

```
sub-<subject_label>/[ses-<session_label>/]
  func/
    sub-<subject_label>_[specifiers]_desc-hmc_boldref.nii.gz
    sub-<subject_label>_[specifiers]_from-orig_to_boldref_mode-image_desc-hmc_xfm.nii.gz
```

:::{note}
Motion correction outputs are part of the *minimal* processing level.
:::

#### Coregistration outputs

Registration of the BOLD series to the anatomical reference image generates a further reference
image and affine transform

```
sub-<subject_label>/[ses-<session_label>/]
  func/
    sub-<subject_label>_[specifiers]_desc-coreg_boldref.nii.gz
    sub-<subject_label>_[specifiers]_from-boldref_to-<anat>_mode-image_desc-coreg_xfm.txt
```

:::{note}
Coregistration outputs are part of the *minimal* processing level.
:::

#### Fieldmap registration

If a fieldmap is used for the correction of a BOLD series, then a registration
is calculated between the BOLD series and the fieldmap. If, for example, the fieldmap
is identified with `"B0FieldIdentifier": "TOPUP"`, the generated transform will be named

```
sub-<subject_label>/[ses-<session_label>/]
  func/
    sub-<subject_label>_[specifiers]_from-boldref_to-TOPUP_mode-image_xfm.nii.gz
```

If the association is discovered through the `IntendedFor` field of the
fieldmap metadata, then the transform will be given an auto-generated name

```
sub-<subject_label>/[ses-<session_label>/]
  func/
    sub-<subject_label>_[specifiers]_from-boldref_to-auto000XX_mode-image_xfm.txt
```

:::{note}
Fieldmap registration outputs are part of the *minimal* processing level.
:::

#### Regularly gridded outputs (images)

Volumetric output spaces labels (`<space_label>` above, and in the following) include
`T1w` and `MNI152NLin2009cAsym` (default).

#### Surfaces, segmentations and parcellations from FreeSurfer

If FreeSurfer reconstructions are used, the `(aparc+)aseg` segmentations are aligned to the
subject's anatomical reference space and resampled to the BOLD grid, and the BOLD series are
resampled to the mid-thickness surface mesh

```
sub-<subject_label>/[ses-<session_label>/]
  func/
    sub-<subject_label>_[specifiers]_space-<anat>_desc-aparcaseg_dseg.nii.gz
    sub-<subject_label>_[specifiers]_space-<anat>_desc-aseg_dseg.nii.gz
    sub-<subject_label>_[specifiers]_hemi-[LR]_space-<space_label>_bold.func.gii
```

Surface output spaces include `fsnative` (full density subject-specific mesh),
`fsaverage` and the down-sampled meshes `fsaverage6` (41k vertices) and
`fsaverage5` (10k vertices, default).

#### Grayordinates files

[CIFTI](https://www.nitrc.org/forum/attachment.php?attachid=333&group_id=454&forum_id=1955) is
a container format that holds both volumetric (regularly sampled in a grid) and surface
(sampled on a triangular mesh) samples.
Sub-cortical time series are sampled on a regular grid derived from one MNI template, while
cortical time series are sampled on surfaces projected from the [^Glasser2016] template.
If CIFTI outputs are requested (with the `--cifti-outputs` argument), the BOLD series are also
saved as `dtseries.nii` CIFTI2 files

```bash
sub-<subject_label>/[ses-<session_label>/]
  func/
    sub-<subject_label>_[specifiers]_bold.dtseries.nii
```

CIFTI output resolution can be specified as an optional parameter after `--cifti-output`.
By default, '91k' outputs are produced and match up to the standard `HCP Pipelines`_ CIFTI
output (91282 grayordinates @ 2mm). However, '170k' outputs are also possible, and produce
higher resolution CIFTI output (170494 grayordinates @ 1.6mm).

#### Extracted confounding time series

For each {abbr}`BOLD (blood-oxygen level dependent)` run processed with *NiBabies*, an
accompanying *confounds* file will be generated.
Confounds_ are saved as a {abbr}`TSV (tab-separated value)` file::

```
sub-<subject_label>/[ses-<session_label>/]
  func/
    sub-<subject_label>_[specifiers]_desc-confounds_timeseries.tsv
    sub-<subject_label>_[specifiers]_desc-confounds_timeseries.json
```

These {abbr}`TSV (tab-separated values)` tables look like the example below,
where each row of the file corresponds to one time point found in the
corresponding {abbr}`BOLD (blood-oxygen level dependent)` time series

| csf | white_matter | global_signal | std_dvars | dvars | framewise_displacement | t_comp_cor_00 | t_comp_cor_01 | t_comp_cor_02 | t_comp_cor_03 | t_comp_cor_04 | t_comp_cor_05 | a_comp_cor_00 | a_comp_cor_01 | a_comp_cor_02 | a_comp_cor_03 | a_comp_cor_04 | a_comp_cor_05 | non_steady_state_outlier00 | trans_x | trans_y | trans_z | rot_x | rot_y | rot_z |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 682.75275 | 0.0 | 491.64752000000004 | n/a | n/a | n/a | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 | 1.0 | 0.0 | 0.0 | 0.0 | -0.00017029 | -0.0 | 0.0 |
| 669.14166 | 0.0 | 489.4421 | 1.168398 | 17.575331 | 0.07211929999999998 | -0.4506846719 | 0.1191909139 | -0.0945884724 | 0.1542023065 | -0.2302324641 | 0.0838194238 | -0.032426848599999995 | 0.4284323184 | -0.5809158299 | 0.1382414008 | -0.1203486637 | 0.3783661265 | 0.0 | 0.0 | 0.0207752 | 0.0463124  | -0.000270924 | -0.0 | 0.0 |

#### Multi-echo derivatives

For multi-echo datasets, the output `_bold` series are "optimally combined" by
`tedana`_ to better estimate the BOLD signal.
This process also generates a T2\* map, which is resampled into every requested output
space.

```bash
sub-<subject_label>/[ses-<session_label>/]
  func/
    sub-<subject_label>_[specifiers]_T2starmap.nii.gz
```

If the `--me-output-echos` flag is specified, then the distortion-corrected (STC, HMC, SDC)
per-echo time series are output. For example, if the inputs are of the form::

```bash
sub-01/
  func/
    sub-01_task-rest_echo-1_bold.nii.gz
    sub-01_task-rest_echo-2_bold.nii.gz
    sub-01_task-rest_echo-3_bold.nii.gz
```

Then the output will include::

```bash
sub-01/
  func/
    sub-01_task-rest_boldref.nii.gz
    sub-01_task-rest_desc-brain_mask.nii.gz
    sub-01_task-rest_echo-1_desc-preproc_bold.nii.gz
    sub-01_task-rest_echo-2_desc-preproc_bold.nii.gz
    sub-01_task-rest_echo-3_desc-preproc_bold.nii.gz
```

These may then be used independently with multi-echo tools, such as `tedana`_,
to perform more advanced denoising or alternative combination strategies.

:::{danger}
Slice timing correction in *NiBabies* is referenced to the middle slice by default,
which leads to a time shift in the volume onsets by 0.5 TR (repetition time).
For example, assuming a TR of 2s, original onsets of 0, 2, and 4s would be shifted
to 1, 3, and 5s, respectively.
In case you did execute slice timing correction, you must check that subsequent
analyses (e.g., general linear modeling) consider the right onset shifts.
For example, when specifying a first-level model, you should set parameters in your
software package or first-level model function accordingly (e.g., select the middle
slice as reference).
Alternatively, you could manually adjust the volume onsets (e.g. as mentioned in
the example above from [0, 2, 4] to [1, 3, 5]) or the event onsets accordingly.
In contrast to volume onsets, event onsets need to be shifted *backward* by half a TR,
for example, from [5, 10, 15] to [4, 9, 14].

Further information on this issue is found at
[this blog post (with thanks to Russell Poldrack and Jeanette Mumford)](https://reproducibility.stanford.edu/slice-timing-correction-in-fmriprep-and-linear-modeling/).
:::

## Confounds

The {abbr}`BOLD (blood-oxygen level dependent)` signal measured with fMRI is a mixture of fluctuations
of both neuronal and non-neuronal origin.
Neuronal signals are measured indirectly as changes in the local concentration of oxygenated hemoglobin.
Non-neuronal fluctuations in fMRI data may appear as a result of head motion, scanner noise,
or physiological fluctuations (related to cardiac or respiratory effects).
For a detailed review of the possible sources of noise in the BOLD signal, refer to [^Greve2013].

*Confounds* (or nuisance regressors) are variables representing fluctuations with a potential
non-neuronal origin.
Such non-neuronal fluctuations may drive spurious results in fMRI data analysis,
including standard activation {abbr}`GLM (General Linear Model)` and functional connectivity analyses.
It is possible to minimize confounding effects of non-neuronal signals by including
them as nuisance regressors in the GLM design matrix or regressing them out from
the fMRI data - a procedure known as *denoising*.
There is currently no consensus on an optimal denoising strategy in the fMRI community.
Rather, different strategies have been proposed, which achieve different compromises between
how much of the non-neuronal fluctuations are effectively removed, and how much of neuronal fluctuations
are damaged in the process.
The *NiBabies* pipeline generates a large array of possible confounds.

The most well established confounding variables in neuroimaging are the six head-motion parameters
(three rotations and three translations) - the common output of the head-motion correction
(also known as *realignment*) of popular fMRI preprocessing software
such as SPM_ or FSL_.
Beyond the standard head-motion parameters, the NiBabies pipeline generates a large array
of possible confounds, which enable researchers to choose the most suitable denoising
strategy for their downstream analyses.

Confounding variables calculated in *NiBabies* are stored separately for each subject,
session and run in {abbr}`TSV (tab-separated value)` files - one column for each confound variable.
Such tabular files may include over 100 columns of potential confound regressors.

:::{danger}
Do not include all columns of `~_desc-confounds_timeseries.tsv` table
into your design matrix or denoising procedure.
Filter the table first, to include only the confounds (or components thereof)
you want to remove from your fMRI signal.
The choice of confounding variables may depend on the analysis you want to perform,
and may be not straightforward as no gold standard procedure exists.
For a detailed description of various denoising strategies and their performance,
see [Parkes2018] and [Ciric2017].
:::

### Confound regressors description

**Basic confounds**. The most commonly used confounding time series:

- Estimated head-motion parameters:
  `trans_x`, `trans_y`, `trans_z`, `rot_x`, `rot_y`, `rot_z` - the 6 rigid-body motion
  parameters (3 translations and 3 rotation), estimated relative to a reference image;

- Global signals:
  - `csf` - the average signal within anatomically-derived eroded {abbr}`CSF (cerebro-spinal fluid)` mask;
  - `white_matter` - the average signal within  the anatomically-derived eroded {abbr}`WM (white matter)` masks;
  - `global_signal` -  the average signal within the brain mask.

#### Parameter expansion of basic confounds

The standard six-motion parameters may not account for all the variance related
to head-motion.
[Friston1996] and [Satterthwaite2013] proposed an expansion of the six fundamental
head-motion parameters.
To make this technique more accessible, *NiBabies* automatically calculates motion parameter
expansion [Satterthwaite2013], providing time series corresponding to the first
*temporal derivatives* of the six base motion parameters, together with their
*quadratic terms*, resulting in the total 24 head motion parameters
(six base motion parameters + six temporal derivatives of six motion parameters +
12 quadratic terms of six motion parameters and their six temporal derivatives).
Additionally, *NiBabies* returns temporal derivatives and quadratic terms for the
three global signals (`csf`, `white_matter` and `global_signal`)
to enable applying the 36-parameter denoising strategy proposed by [Satterthwaite2013].

Derivatives and quadratic terms are stored under column names with
suffixes: `_derivative1` and powers `_power2`.
These are calculated for head-motion estimates (`trans_` and `rot_`) and global signals
(`white_matter`, `csf`, and `global_signal`).

#### Outlier detection

These confounds can be used to detect potential outlier time points -
frames with sudden and large motion or intensity spikes.

- `framewise_displacement` - is a quantification of the estimated bulk-head motion calculated using
  formula proposed by [^Power2012];
- `rmsd` - is a quantification of the estimated relative (frame-to-frame) bulk head motion
  calculated using the {abbr}`RMS (root mean square)` approach of [^Jenkinson2002];
- `dvars` - the derivative of RMS variance over voxels (or {abbr}`DVARS (derivative of
  RMS variance over voxels)`) [^Power2012];
- `std_dvars` - standardized {abbr}`DVARS (derivative of RMS variance over voxels)`;
- `non_steady_state_outlier_XX` - columns indicate non-steady state volumes with a single
  `1` value and `0` elsewhere (*i.e.*, there is one `non_steady_state_outlier_XX` column per
  outlier/volume).

Detected outliers can be further removed from time series using methods such as:
volume *censoring* - entirely discarding problematic time points [^Power2012],
regressing signal from outlier points in denoising procedure, or
including outlier points in the subsequent first-level analysis when building
the design matrix.
Averaged value of confound (for example, mean `framewise_displacement`)
can also be added as regressors in group level analysis [^Yan2013].
*Regressors of motion spikes* for outlier censoring are generated from within *NiBabies*,
and their calculation may be adjusted with the command line options `--fd-spike-threshold`
and `--dvars-spike-threshold` (defaults are FD > 0.5 mm or DVARS > 1.5).
Regressors of motion spikes are stored in separate `motion_outlier_XX` columns.

#### Discrete cosine-basis regressors

Physiological and instrumental (scanner) noise sources are generally present in fMRI
data, typically taking the form of low-frequency signal drifts.
To account for these drifts, temporal high-pass filtering is the immediate option.
Alternatively, low-frequency regressors can be included in the statistical model to account
for these confounding signals.
Using the {abbr}`DCT (discrete cosine transform)` basis functions, *NiBabies* generates
these low-frequency predictors:

- `cosine_XX` - DCT-basis regressors.

One characteristic of the cosine regressors is that they are identical for two different
datasets with the same {abbr}`TR (repetition time)` and the same *effective* number of
time points (*effective* length).
It is relevant to mention *effective* because initial time points identified as *nonsteady
states* are removed before generating the cosine regressors.

:::{caution}
If your analysis includes separate high-pass filtering, do not include
`cosine_XX` regressors in your design matrix.
:::

:::{seealso}
- A detailed explanation about temporal high-pass filtering is provided with
  the [BrainVoyager User Guide](https://www.brainvoyager.com/bvqx/doc/UsersGuide/Preprocessing/TemporalHighPassFiltering.html).
- [This comment](https://github.com/nipreps/fmriprep/issues/1899#issuecomment-561687460)
  on an issue regarding CompCor regressors.
:::

#### CompCor confounds

{abbr}`CompCor (Component Based Noise Correction)` is a {abbr}`PCA (principal component analysis)`,
hence component-based, noise pattern recognition method.
In the method, principal components are calculated within an {abbr}`ROI (Region of Interest)`
that is unlikely to include signal related to neuronal activity, such as {abbr}`CSF (cerebro-spinal fluid)`
and {abbr}`WM (white matter)` masks.
Signals extracted from CompCor components can be further regressed out from the fMRI data with a
denoising procedure [^Behzadi2007].

- `a_comp_cor_XX` - additional noise components are calculated using anatomical {abbr}`CompCor
  (Component Based Noise Correction)`;
- `t_comp_cor_XX` - additional noise components are calculated using temporal {abbr}`CompCor
  (Component Based Noise Correction)`.

Four separate CompCor decompositions are performed to compute noise components: one temporal
decomposition (`t_comp_cor_XX`) and three anatomical decompositions (`a_comp_cor_XX`) across
three different noise ROIs: an eroded white matter mask, an eroded CSF mask, and a combined mask derived
from the union of these.

Each confounds data file will also have a corresponding metadata file
(`~desc-confounds_regressors.json`).
Metadata files contain additional information about columns in the confounds TSV file:

```json

{
  "a_comp_cor_00": {
    "CumulativeVarianceExplained": 0.1081970825,
    "Mask": "combined",
    "Method": "aCompCor",
    "Retained": true,
    "SingularValue": 25.8270209974,
    "VarianceExplained": 0.1081970825
  },
  "dropped_0": {
    "CumulativeVarianceExplained": 0.5965809597,
    "Mask": "combined",
    "Method": "aCompCor",
    "Retained": false,
    "SingularValue": 20.7955177198,
    "VarianceExplained": 0.0701465624
  }
}
```

For CompCor decompositions, entries include:

- `Method`: anatomical or temporal CompCor.
- `Mask`: denotes the {abbr}`ROI (region of interest)` where the decomposition that generated
  the component was performed: `CSF`, `WM`, or `combined` for anatomical CompCor.
- `SingularValue`: singular value of the component.
- `VarianceExplained`: the fraction of variance explained by the component across the decomposition ROI mask.
- `CumulativeVarianceExplained`: the total fraction of variance explained by this particular component
  and all preceding components.
- `Retained`: Indicates whether the component was saved in `desc-confounds_timeseries.tsv`
  for use in denoising.
  Entries that are not saved in the data file for denoising are still stored in metadata with the
  `dropped` prefix.

:::{caution}
Only a subset of these CompCor decompositions should be used for further denoising.
The original Behzadi aCompCor implementation [^Behzadi2007] can be applied using
components from the combined masks, while the more recent Muschelli implementation
[^Muschelli2014] can be applied using
the {abbr}`WM (white matter)` and {abbr}`CSF (cerebro-spinal fluid)` masks.
To determine the provenance of each component, consult the metadata file (described above).

There are many valid ways of selecting CompCor components for further denoising.
In general, the components with the largest singular values (i.e., those that
explain the largest fraction of variance in the data) should be selected.
*NiBabies* outputs components in descending order of singular value.
Common approaches include selecting a fixed number of components (e.g., the
first 5 or 6), using a quantitative or qualitative criterion (e.g., elbow, broken
stick, or condition number), or using sufficiently many components that a minimum
cumulative fraction of variance is explained (e.g., 50%).
:::

:::{caution}
Similarly, if you are using anatomical or temporal CompCor it may not make sense
to use the `csf`, or `white_matter` global regressors -
see [fmriprep#1049](https://github.com/nipreps/fmriprep/issues/1049).
Conversely, using the overall `global_signal` confound in addition to CompCor's
regressors can be beneficial (see [^Parkes2018]).
:::

:::{danger}
*NiBabies* does high-pass filtering before running anatomical or temporal CompCor.
Therefore, when using CompCor regressors, the corresponding `cosine_XX` regressors
should also be included in the design matrix.
:::

:::{seealso}
This didactic [discussion on NeuroStars.org](https://neurostars.org/t/fmrirep-outputs-very-high-number-of-acompcors-up-to-1000/5451)
where Patrick Sadil gets into details about PCA and how that base technique applies
to CompCor in general and *fMRIPrep*'s implementation in particular.
:::

#### Confounds estimated from the brain's outer edge

Reusing the implementation of aCompCor, *NiBabies* generates regressors corresponding to the
24 first principal components extracted with PCA using the voxel time-series delineated by
the brain's outer edge (*crown*) mask.
The procedure essentially follows the initial proposal of the approach by Patriat et al.
[^Patriat2017] and is described in our ISMRM abstract [^Provins2022].

#### Confounds and "carpet"-plot on the visual reports

The visual reports provide several sections per task and run to aid designing
a denoising strategy for subsequent analysis.
Some of the estimated confounds are plotted with a "carpet" visualization of the
{abbr}`BOLD (blood-oxygen level-dependent)` time series [^Power2016].

Noise components computed during each CompCor decomposition are evaluated according
to the fraction of variance that they explain across the nuisance ROI.
This is used by *NiBabies* to determine whether each component should be saved for
use in denoising operations: a component is saved if it contributes to explaining
the top 50 percent of variance in the nuisance ROI.
*NiBabies* can be configured to save all components instead using the command line
option `--return-all-components`.
*NiBabies* reports include a plot of the cumulative variance explained by each
component, ordered by descending singular value.

Also included is a plot of correlations among confound regressors.
This can be used to guide selection of a confound model or to assess the extent
to which tissue-specific regressors correlate with global signal.

See implementation on <python#nibabies.workflows.bold.confounds.init_bold_confs_wf>.

## References

[^Behzadi2007]: Behzadi Y, Restom K, Liau J, Liu TT.
A component-based noise correction method (CompCor) for BOLD and perfusion-based fMRI. NeuroImage. 2007.
doi:[10.1016/j.neuroimage.2007.04.042](https://doi.org/10.1016/j.neuroimage.2007.04.042)

[^Ciric2017]: Ciric R, Wolf DH, Power JD, Roalf DR, Baum GL, Ruparel K, Shinohara RT, Elliott MA, Eickhoff SB, Davatzikos C., Gur RC, Gur RE, Bassett DS, Satterthwaite TD.
Benchmarking of participant-level confound regression strategies for the control of motion artifact in studies of functional connectivity. 
Neuroimage. 2017.
doi:[10.1016/j.neuroimage.2017.03.020](https://doi.org/10.1016/j.neuroimage.2017.03.020)

[^Greve2013]: Greve DN, Brown GG, Mueller BA, Glover G, Liu TT.
A Survey of the Sources of Noise in fMRI.
Psychometrika. 2013.
doi:[10.1007/s11336-013-9344-2](https://doi.org/10.1007/s11336-013-9344-2)

[^Friston1996]: Friston KJ1, Williams S, Howard R, Frackowiak RS, Turner R.
Movement‐Related effects in fMRI time‐series.
Magnetic Resonance in Medicine. 1996.
doi:[10.1002/mrm.191035031](https://doi.org/10.1002/mrm.1910350312)

[^Glasser2016]: Glasser MF, Coalson TS Robinson EC, Hacker CD, Harwell J, Yacoub E, Ugurbil K, Andersson J, Beckmann CF, Jenkinson M, Smith SM, Van Essen DC.
A multi-modal parcellation of human cerebral cortex.
Nature. 2016.
doi:[10.1038/nature18933](https://doi.org/10.1038/nature18933)

[^Jenkinson2002]: Jenkinson M, Bannister P, Brady M, Smith S.
Improved optimization for the robust and accurate linear registration and motion correction of brain images.
Neuroimage. 2002.
doi:[10.1016/s1053-8119(02)91132-8](https://doi.org/10.1016/s1053-8119(02)91132-8)

[^Muschelli2014]: Muschelli J, Nebel MB, Caffo BS, Barber AD, Pekar JJ, Mostofsky SH.
Reduction of motion-related artifacts in resting state fMRI using aCompCor.
NeuroImage. 2014.
doi:[10.1016/j.neuroimage.2014.03.028](https://doi.org/10.1016/j.neuroimage.2014.03.028)

[^Parkes2018]: Parkes L, Fulcher B, Yücel M, Fornito A.
An evaluation of the efficacy, reliability, and sensitivity of motion correction strategies for resting-state functional MRI.
NeuroImage. 2018.
doi:[10.1016/j.neuroimage.2017.12.073](https://doi.org/10.1016/j.neuroimage.2017.12.073)

[^Patriat2015]: Patriat R, EK Molloy, RM Birn, T. Guitchev, and A. Popov.
Using Edge Voxel Information to Improve Motion Regression for Rs-FMRI Connectivity Studies.
Brain Connectivity. 2015.
doi:[10.1089/brain.2014.0321](https://doi.org/10.1089/brain.2014.0321)

[^Patriat2017]: Patriat R, Reynolds RC, Birn RM,
An improved model of motion-related signal changes in fMRI.
NeuroImage. 2017.
doi:[10.1016/j.neuroimage.2016.08.051](https://doi.org/10.1016/j.neuroimage.2016.08.051)

[^Power2012]: Power JD, Barnes KA, Snyder AZ, Schlaggar BL, Petersen, SA.
Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion.
NeuroImage. 2012.
doi:[10.1016/j.neuroimage.2011.10.018](https://doi.org/10.1016/j.neuroimage.2011.10.018=[]

[^Power2016]: Power JD.
A simple but useful way to assess fMRI scan qualities.
NeuroImage. 2016.
doi:[10.1016/j.neuroimage.2016.08.009](https://doi.org/10.1016/j.neuroimage.2016.08.009)

[^Provins2022]: Provins C et al.
Quality control and nuisance regression of fMRI, looking out where signal should not be found.
Proc. Intl. Soc. Mag. Reson. Med. 31, London (UK). 2022
doi:[10.31219/osf.io/hz52v](https://doi.org/10.31219/osf.io/hz52v)

[^Satterthwaite2013]: Satterthwaite TD, Elliott MA, Gerraty RT, Ruparel K, Loughead J, Calkins ME, Eickhoff SB, Hakonarson H, Gur RC, Gur RE, Wolf DH.
An improved framework for confound regression and filtering for control of motion artifact in the preprocessing of resting-state functional connectivity data.
NeuroImage. 2013.
doi:[10.1016/j.neuroimage.2012.08.052](https://doi.org/10.1016/j.neuroimage.2012.08.052)

[^Yan2013]: Yan CG, Cheung B, Kelly C, Colcombe S, Craddock RC, Di Martino A, Li Q, Zuo XN, Castellanos FX, Milham MP.
A comprehensive assessment of regional variation in the impact of head micromovements on functional connectomics.
NeuroImage. 2013.
doi:[10.1016/j.neuroimage.2013.03.004](https://doi.org/10.1016/j.neuroimage.2013.03.004)
