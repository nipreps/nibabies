25.2.0 (September 22, 2025)
===========================
The first release of the 25.2.x series.

This release synchronizes with downstream dependencies, as well as
leverages TemplateFlow to retrieve intermediate xfms used during multi-step registration.

Thanks to @LuciMoore for the contribution!

### Documentation
  * DOC: MCRIBS and surf recon methods chosen based on age (#485)

### Enhancements
  * ENH: Retrieve transforms with templateflow (#486)
  * ENH: Save cortex mask (#487)
  * ENH: Add multiverse output layout (#481)

### Bug Fixes
  * FIX: Default to grid-constant when resampling BOLD to output spaces

25.1.2 (August 28, 2025)
========================
A patch release in the 25.1.x series.

Fixes compatibility with Python 3.12 when running MCRIBS.

### Bug Fixes
  * FIX: MCRIBS compatibility (#483)

### Internals / Maintenance
  * MAINT: Fix outputs check (#479)

25.1.1 (July 2, 2025)
=====================
A patch release in the 25.1.x series.

Includes a fix to the fieldmap reference orientation to match the B-spline coefficients.

Thanks to @joey-scanga for the contribution!

### Bug Fixes
  * FIX: Orient fieldmap before checking spline fit (#475)


25.1.0 (June 24, 2025)
======================
First release of the 25.1.x series. A few key changes include a new workflow
for derivatives compatibility when an anatomical template is not present, and performing
a two step registration (currently only for MNI152NLin6Asym) is now the default behavior.
This can be disabled by adding `--no-multi-step-reg` to the command.

### Enhancements
  * ENH: Verify derivatives are compatible with anatomical reference (#459)
  * ENH: Make `--multi-step-reg` a boolean action, enable by default (#470)
  * ENH: Convert pooch retrieval to interface, allow setting cache dir (#467)

### Internals / Maintenance
  * MAINT: bump dependencies, test on python 3.13 (#468)
  * MAINT: deprecate `--clean-workdir` (#473)


25.0.2 (May 19, 2025)
---------------------
A patch release including a couple key bug fixes:
 - Adds missing dependency to Docker image when BOLD coregistration falls back to FSL
 - Fixes connections to allow fieldmap-less SDC

Thanks to @joey-scanga for the contribution!

### Documentation
  * DOC: Changed "B0Identifier" -> "B0FieldIdentifier" (#456)

### Enhancements
  * ENH: Add message to show which files are to be used (#452)
  * ENH: Tag anatomical workflows (#463)

### Bug Fixes
  * FIX: SyN workflow (#453)
  * FIX: Add Convert3d to Docker image (#461)

### Internals / Maintenance
  * MAINT: Support sdcflows@main and smriprep >=0.18 (#462)


25.0.1 (February 5, 2025)
-------------------------
A bug-fix release to address an issue when using a precomputed mask for MCRIBS surface reconstruction without an anatomical template.

  * FIX: Match mask header prior to n4 correction (#443)


25.0.0 (January 28, 2025)
=========================
A new minor release with improvements to anatomical to template spatial normalization.

Registration will now prioritize the same modality as the anatomical template, if available.

A new flag `--norm-csf` performs CSF normalization on the anatomical template prior to template registration.

A new flag `--multi-step-reg` adds an intermediate step when registering to MNI152NLin6Asym, first performing anatomical -> MNIInfant:cohort-X (age matched by default), and then concatenates the transform with an already computed MNIInfant -> MNI152NLin6Asym.

Both of the new flags above are disabled by default, but have shown promise and may become defaults in the next release. Please experiment with your data, and any feedback on the results would be greatly appreciated!


### Enhancements
  * ENH: Output anatomical coregistration transform + report (#437)
  * ENH: Minimize clipping prior to surface reconstruction with MCRIBS (#436)
  * ENH: Output fsLR meshes on subject surfaces (#427)
  * ENH: Add flag for multi-step registration to adult templates (#415) (#425) (#430) (#433)
  * ENH: Option to normalize CSF prior to template registration (#419)
  * ENH: Expand template registration to use either anatomical modality (#418)

### Bug Fixes
  * FIX: Reduce range that --surface-recon-method auto recommends MCRIBS (#438)
  * FIX: Allow T2 only without the use of --derivatives
  * FIX: New styling catches (#417)
  * FIX: Default surface recon method should be None (#416)

### Internals / Maintenance
  * TST: Build workflow across different conditions (#409)
  * MAINT: Remove deprecated parser arguments (#407)


24.1.0 (October 02, 2024)
=========================
This new minor release includes a few bug fixes, such as excluding MCRIBS from surface reconstruction without a precomputed segmentation and ensuring generated derivatives are not masked, as well as improvements to reporting.

### Enhancements
  * ENH: Add boilerplate, errors to report (#403)
  * ENH: Add age to session report (#402)
  * ENH: Improvements to age parsing (#395, #398)

### Bug Fixes
  * FIX: MCRIBS auto surface reconstruction logic (#399)
  * FIX: Do not force masking of anatomicals when using `--derivatives` (#400)

### Internals / Maintenance
  * MAINT: Revisit warnings filter (#396)
  * MAINT: Automate testing with tox (#404)
  * MAINT: Port over parser arguments and tests from fmriprep (#401)


24.0.1 (August 31, 2024)
========================
A patch release with a fix for the BOLD T2\* workflow. The command line argument  `--me-t2s-fit-method` was added for finer control when processing multi-echo datasets.

* FIX: Add missing me-t2s-fit-method option (#385)
* DOC: Reformat abbreviations (#386)


24.0.0 (August 29, 2024)
========================
This major release includes a substantial refactoring of the pipeline.

One key addition is the addition of the `--level` flag, which can take the arguments minimal, resampling or full. The default is full, which should produce nearly the same results as previous versions. minimal will produce only the minimum necessary to deterministically generate the remaining derivatives. resampling will produce some additional derivatives, intended to simplify resampling with other tools.

The `--derivatives` flag was altered to take arguments in the form `name=/path/to/dir`.
For each directory provided, if a derivative is found - it will be used instead of computing it from scratch. If a derivative is not found, NiBabies will compute it and proceed as usual.

Taken together, these features can allow a dataset provider to run a minimal NiBabies run, targeting many output spaces, while a user can then run a `--derivatives` run to generate additional derivatives in only the output spaces they need. Another use case is to provide an precomputed derivative to override the default NiBabies behavior, enabling easier workarounds for bugs or experimentation with alternatives.

Another new feature is a dynamic anatomical reference, which is set based on surface reconstruction method or through the `--reference-anatomical` flag. Previously, T1w was the default output space. Now, the reference anatomical is determined based on the surface reconstruction method.

Additionally, minor adjustments have been made to MCRIBS surface reconstruction to address failure rates. This is still an on-going investigation, but preliminary results look promising.

This release resolves a number of issues with fieldmaps inducing distortions during correction. Phase difference and direct fieldmaps are now masked correctly, preventing the overestimation of distortions outside the brain. Additionally, we now implement Jacobian weighting during unwarping, which corrects for compression and expansion effects on signal intensity. To disable Jacobian weighting, add `fmap-jacobian` to the `--ignore` argument.

Finally, a new resampling method has been added, to better account for susceptibility distortion and motion in a single shot resampling to a volumetric target space. We anticipate extending this to surface targets in the future.

  * FIX: nest pathlib import in fix_multi_source_name (#365)
  * FIX: Avoid retrieving multiple templates from latest TF (#353)
  * FIX: Raise informative error if no t1w or t2w found (#347)
  * FIX: Easier pyenv usage (#342)
  * FIX: Catch nonexistent derivatives, clean up subworkflow logic (#336)
  * FIX: Use fsLR reg sphere for MCRIBS morphometrics resampling (#334)
  * FIX: T2star map MNI scaling (#320)

  * ENH: Alter outputs when MCRIBS reconstruction is used (#329)
  * ENH: Use nireports for Report generation + add reportlet per reconstruction (#328)
  * ENH: better repr for Derivatives class (#351)

  * RF: Move to fit/apply workflow (#360)
  * RF: Replace `resource_filename` with `load_data` (#345)

  * MAINT: Bump urllib3 from 2.0.3 to 2.0.7 (#319)
  * MAINT: Raise minimum to 3.10, bump actions (#337)
  * MAINT: Bump pillow from 9.5.0 to 10.0.1 (#317)
  * MAINT: Update to latest migas API (#326)

  * DOC: Use correct argument flag (#338)
  * DOC: Move to new theme, add outputs description (#383)


23.1.0 (November 22, 2023)
============
The next minor release of *NiBabies*, this release includes a number of new goodies, including:

### New surface reconstruction option
M-CRIB-S (Adamson et al., https://www.nature.com/articles/s41598-020-61326-2), has shown to improve performance in participants under 9 months. If you would like to try this method, add the following to your command: `--surface-recon-method mcribs`.

Note: Currently, a T2w image and pre-computed segmentation derivative must be provided to run mcribs.

### Improved batch processing
*NiBabies* now automatically parses the BIDS directory for participant ages, first searching in the
participant's `session.tsv`, and falling back to `participants.tsv`. This simplifies batch submissions including multiple subjects & sessions. As a result, the `--age-months` flag has been deprecated, and will be removed in a later release.

### Goodvoxels projection
An option to determine and exclude high-variance voxels from being projected to the surface when creating CIFTI files. To enable this, add `--project-goodvoxels` to your command.

### Single anatomical processing
Running *NiBabies* is now less restrictive, and will still process data missing either a T1w / T2w image. However, for best results, it is recommended to collect and include both for processing.

### Anat-specific derivatives inputs
Previous, *NiBabies* expected input from the `--derivatives` flag to be in T1w space, using the entity `space-orig`. This has now been changed to support derivatives in either T1w or T2w space. For more information, please see https://nibabies.readthedocs.io/en/23.1.0/faqs.html#leveraging-precomputed-results


## Full Changelog
  * CI: Purge codecov python package (#282)
  * DKR: Upgrade Docker base, c3d (#275)
  * DKR: Add M-CRIB-S to Docker container (#283)
  * DKR: Update dependencies, split into multi-stage build
  * ENH: Add option to exclude projecting high variance voxels to surface (#278)
  * ENH: Resample morphometrics to fsLR CIFTI-2 files when outputting CIFTIs (#279)
  * ENH: Add MCRIBReconAll as alternative surface reconstruction method (#283)
  * ENH: Reorder anatomical processing, run ANTs DenoiseImage on anatomicals (#286)
  * ENH: Extract participant ages from BIDS sources, deprecate `--age-months` (#287)
  * ENH: Dilate BOLD mask by 2 voxels to prevent over-aggressive masking degrading T2star map estimation (#296)
  * ENH: Allow precomputed derivatives in T1w or T2w space (#305)
  * ENH: Add separate workflow for single anatomical processing (#316)
  * FIX: Improve free memory estimation (#284)
  * FIX: Ensure age is extracted from sessions file (#291)
  * FIX: Restore CIFTI medial wall masking, subcortical volume LAS reorientation (#298)
  * FIX: Recify "goodvoxels" surface projection (#301)
  * FIX: Connect derivatives mask to mcribs recon (#323)
  * MAINT: Drop TemplateFlowSelect patches (#290)

23.0.0 (January 23, 2023)
=========================
New year, new *NiBabies* minor series!
Some of the highlights of this release include:
- New run-wise BOLD reference generation, prioritizing single-band references if available, unless avoided with the `--ignore sbrefs` flag.
- New output: Preprocessed T2w in T1w space.

A full list of changes can be found below.

## Full Changelog
  * ENH: Runwise bold reference generation (#268)
  * ENH: Add preprocessed T2w volume to outputs (#271)
  * MAINT: Drop versioneer for hatch backend, fully embrace pyproject.toml (#265)
  * MAINT: Rotate CircleCI secrets and setup up org-level context (#266)
  * CI: Bump convenience images, limit datalad (#267)
  * FIX: Remove legacy CIFTI variant support (#264)


22.2.0 (December 13, 2022)
==========================
The final *NiBabies* minor series of 2022!
Some highlights of the new additions in this release series includes:
- surface morphometrics outputs, including cortical thickness
- T2star maps for multiecho data, projected to target output spaces

This series will be the last to support Python 3.7.

A full list of changes can be found below.

## Full Changelog
  * FIX: Remove cortex masking during vol2surf sampling (#258)
  * ENH: Improve migas telemetry (#257)
  * CI: GitHub actions update (#256)
  * ENH: Add morphometric outputs (#255)
  * ENH: Output T2star maps for multiecho data (#252)
  * FIX: Use the binarized output from the brain extraction (#246)
  * DOC: Add long description including background/significance (#243)
  * CI: Fix docker credential error (#244)
  * DOC: Advertise nipreps community pages, add section on contributions (#242)

22.1.3 (September 12, 2022)
===========================
This patch release includes a vital fix for susceptibility distortion correction on multi-echo data.

* FIX: Field name for multi-echo fieldmap correction (#233)

22.1.2 (August 22, 2022)
========================
This patch release includes a fix to FreeSurfer version detection, which was causing `recon-all` to use outdated flags.

22.1.1 (August 15, 2022)
========================
A bugfix release that includes missing files needed to run `infant_recon_all`.

## Full Changelog
* FIX: Add missing shared object for `infant_recon_all` (#231)
* RF: `migas` reporting (#230)

22.1.0 (August 2, 2022)
=======================
A new minor release! The 22.1.x series of *NiBabies* includes:

- Improved alignment between FreeSurfer outputs and processed anatomical.
- Decreased memory usage while running across multiple processes (default).
- Fix to multi-echo processing in cases where an optimally combined file of all echoes was missing.
- Fix to the subcortical CIFTI to be in *LAS* orientation.

## Full Changelog
* FIX: Correct fsnative <-> anatomical transforms (#223)
* FIX: Vastly improve multi-echo handling (#220)
* ENH: Add migas telemetry to nibabies (#226)
* ENH: Add interface for reorienting images (#229)
* DOCKER: Bump Python to 3.9 (#221)
* RF/ENH: Rework workflow generation (#219)


22.0.2 (May 11, 2022)
---------------------
A bug-fix release in the 22.0.x series.

This release includes a fix to a problem where `--cifti-output` was not
producing any outputs.

## Full Changelog
* CI: Force all git-annex dependencies to be installed (#217)
* CI: Simplify config with anchors (#209)
* FIX: Remedy missing CIFTI outputs (#212)
* MAINT: Set maximum scipy for Python 3.7.x (#216)

22.0.1 (April 6, 2022)
----------------------
A bug-fix release in the 22.0.x series.

This release includes a fix for when using `UNCInfant` as an output space,
as well as a few improvements to susceptibility distortion correction (SDC).
These include a new flag (`--topup-max-vols`) for controlling the number of
volumes used by TOPUP, and support for SDC in the case where single phase-encoding
fieldmap is used to correct opposite phase-encoding BOLD/EPI runs.

### Warning
This release includes a new version of PyBIDS, which now preserves any
zero-padding within the `run` entity. As a result, NiBabies output naming
may slightly differ from previous versions.

## Full Changelog
* CI: Migrate to token auth when uploading to pypi (#203)
* ENH: Improve fieldmap support (#205)
* MAINT: Bump niworkflows (#208)
* STY: Bump style dependencies, run isort on repo (#206)

22.0.0 (March 28, 2022)
-----------------------
A new `NiBabies` minor series!

This release includes a number of new features, as well as various bug fixes. Some of the biggest changes include:
- [Ability to pass in anatomical derivatives](https://nibabies.readthedocs.io/en/latest/faqs.html#leveraging-precomputed-results).
Users can now leverage a precomputed brain mask and/or discrete anatomical segmentations.
- A new flag `--me-output-echoes` to output individual corrected echo time series.
This is useful when doing additional multi-echo processing.

Thank you for using *NiBabies*!
If you encounter any issues with this release, please let us know
by posting an issue on our GitHub page!

A full list of changes can be found below.

## Full Changelog
* CI: Add workflow smoke tests (#100)
* DOC: Add FAQs page (#164)
* DOCKER: Upgrade to FSL6, use niprep miniconda base layer (#191)
* ENH: Add major/minor version prefix to base workflow name (#202)
* ENH: Add `--me-output-echos` CLI flag (#174)
* ENH: Precomputed derivatives (#173)
* ENH: Validate files passed with `--derivatives` (#182)
* FIX: Clean up generated boilerplate (#200)
* FIX: Various Configuration module touch-ups (#197)
* FIX: Clean up default output space handling (#196)
* FIX: Pandoc citeproc API incompatibility (#195)
* FIX: Check if segmentation directory exists (#165)
* FIX: Update report spec to reflect `infant_recon_all` (#167)
* FIX: ICA Aroma imports (#170)
* FIX: Relabel sub-structures before discarding (#186)
* FIX: Use precomputed aseg within `infant_recon_all` (#184)
* FIX: Remove excess arguments on wrapper tests (#181)
* MAINT: Update versioneer, allow static versioning (#190)
* MAINT: Ensure version is written to version file (#189)
* MAINT: Add missing toml dependency
* MAINT: Add pre-commit checks (#178)
* MAINT: Add RTD config (#173)
* MAINT: Freeze `black` version (#185)
* RF: Wrapper usage logic (#183)
* STY/TEST: Set global style and doctest options (#162)

21.0.2 (November 29, 2021)
--------------------------
A patch release in the 21.0.x series.
This release removes the 24 month age cap for infant recon all processing, as well as includes various small maintenance fixes.

## Full Changelog
* DOC: Use dynamic versioning for examples (#151)
* ENH: Remove infant recon all age cap (#154)
* FIX: Generate boilerplate (#157)
* FIX: Avoid requiring service when checking wrapper version (#159)
* MAINT: Fix dirty version on release (#144)
* MAINT: Only alter pybids config in legacy versions (#152)
* MAINT: Prefetch neonate MNIInfant templates (#159)
* MAINT: Update git-annex version (#159)
* RF: Initialize BIDSLayout with dedicated indexer (#146)

21.0.1 (October 18, 2021)
-------------------------
A patch release in the 21.0.x series.

This patch release is for all **Docker/Singularity** users: `infant_recon_all` did not have all the available templates, which would cause failures for certain ages.
Upgrading to this release will ensure you have all the necessary templates.

## Full Changelog
* DOCKER: Add missing `infant_recon_all` templates (130dcf3)

21.0.0 (October 15, 2021)
-------------------------
The first major release series of 2021.

This release includes enhancements, such as:
 - Fine-grain subcortical alignment during CIFTI generation
 - Improved functional registration to template space
 - Greatly minimized container environment

Additionally, a plethora of bug-fixes are included, and documentation has been improved.

### Caution!
As with all minor version increments, working directories from previous versions **should not be reused**.

### Thank you for using *NiBabies*!
If you encounter any issues with this release, please let us know by posting an issue on our GitHub page!


## Full Changelog
* DOC: Set up external readthedocs documentation (#119) (#126) (#128)
* DOCKER: Reduce container image size (#105) (#133)
* DOCKER: Strip ABI tag from libQt5Core.so.5 (#109)
* DOCKER: Modernize Dockerfile (#85)
* ENH: Port slice timing correction enhancements from fMRIPrep (#137)
* ENH: Change default `--output-layout` to bids (#130)
* ENH: Subcortical alignment workflow (#72)
* ENH: Framewise displacement head radius flag (#104)
* ENH: Incorporate subcortical CIFTI alignment to functional processing (#102)
* ENH: Do not run infant_recon_all if already completed (#101)
* FIX: Handle sessions when grouping BOLDs (#139)
* FIX: Ensure MNIInfant is added if no `--output-spaces` are used (#136)
* FIX: Ensure `nibabies-wrapper` patches are correctly bound (#113)
* FIX: BOLD to template normalization (#99)
* FIX: SDC fieldwarp application (#98)
* FIX: Avoid running BBReg under certain conditions (#95)
* FIX: Standard output spaces (#92)
* FIX: Small Docker environment fixes (#86)
* FIX: Feed NiTransforms with LTAs of type RAS2RAS (#84)
* MAINT: Rename default `infant_recon_all` output directory (#129)
* MAINT: Bump SDCFlows to latest bugfix version (d799fee)
* MAINT: Attempt to pull most recent dev version (#94)
* MAINT: Initial CircleCI workflow (#93)
* STY: `black` nibabies module (#118)

0.1.2 (June 30, 2021)
---------------------
- FIX: BOLD file duplication error

0.1.1 (June 21, 2021)
---------------------
- FIX: BIDS validation error
- ENH: Additional wrapper script (``nibabies-wrapper``) to facilitate running through Docker/Singularity
- DOC: Updates anatomical processing figure, adds functional processing figure
- DOC: Fixes typo in README usage

0.1.0 (May 26, 2021)
--------------------
- Implements susceptibility distortion correction for BOLDs
- Improved documentation
- Adds estimated fieldmap report
- Removed unused commandline arguments


0.0.3 (March 17, 2021)
-----------------------
- Adds functional processing to the workflow
  - minus SDC and fine-grain subcortical CIFTI generation
- Revision of T1w/T2w templates when multiple runs are detected
- Modularization of T1w-T2w coregistration, with increased robustness
- Removed of UNCInfant for a default output space

0.0.2 (Feb 18, 2021)
--------------------
- Improved robustness of structural coregistration
- Removed fMRIPrep dependency
- Bumped Infant FreeSurfer version to latest developmental release
- Added figure overviewing structural workflow

0.0.1 (Nov 13, 2020)
--------------------
- Initial implementation of anatomical pipeline
