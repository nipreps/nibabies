21.0.0 (TBD)
============

  21.0.0rc1 (October 1, 2021)
  ---------------------------

  * DOCKER: Strip ABI tag from libQt5Core.so.5 (#109)
  * MAINT: Bump SDCFlows to latest bugfix version (d799fee)

  21.0.0rc0 (September 29, 2021)
  ------------------------------

  * DOCKER: Reduce image size (#105)
  * ENH: Subcortical alignment workflow (#72)
  * ENH: Framewise displacement head radius flag (#104)
  * ENH: Incorporate subcortical CIFTI alignment to functional processing (#102)
  * ENH: Do not run infant_recon_all if already completed (#101)
  * ENH: Modernize Dockerfile (#85)
  * FIX: BOLD to template normalization (#99)
  * FIX: SDC fieldwarp application (#98)
  * FIX: Avoid running BBReg under certain conditions (#95)
  * FIX: Standard output spaces (#92)
  * FIX: Small Docker environment fixes (#86)
  * FIX: Feed NiTransforms with LTAs of type RAS2RAS (#84)
  * MAINT: Attempt to pull most recent dev version (#94)
  * MAINT: Initial CircleCI workflow (#93)

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