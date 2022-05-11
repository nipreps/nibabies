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
