## Changelog

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