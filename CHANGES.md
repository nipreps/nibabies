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