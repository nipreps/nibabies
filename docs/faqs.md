# Tips and FAQs

## Leveraging precomputed results

Whether to allow for manual intervention for tough cases, or simply to reduce processing time, *NiBabies* can allow the use of certain pre-computed files during processing.
Initial support is limited to the following files:
- Anatomical mask in T1w space
- Antomical segmentation (aseg) in T1w space

To use pre-computed results, one or more [BIDS Derivatives](https://bids-specification.readthedocs.io/en/stable/05-derivatives/01-introduction.html#bids-derivatives) directories must be passed in to *NiBabies* using the `--derivatives` flag.
Derivative directories must include a [`dataset_description.json` and the required fields](https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html#derived-dataset-and-pipeline-description).
Additionally, files must include the `space-orig` key-value pair in the name.

A sample layout of a derivatives directory can be found below:

```
my_precomputed/
├── dataset_description.json
└── sub-01
    └── anat
        ├── sub-01_space-orig_desc-aseg_dseg.nii.gz
        ├── sub-01_space-orig_desc-brain_mask.json
        └── sub-01_space-orig_desc-brain_mask.nii.gz
```

## Multi-atlas segmentation with joint label fusion

By default, *NiBabies* will run [FSL FAST](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FAST) for tissue segmentation, and Infant FreeSurfer for segmentation labels.
However, you can instead use ANTs Joint Label Fusion to generate both, granted you provide multiple atlases with anatomicals / segmentations via the `--segmentation-atlases-dir` flag.
When using this approach, there are a few assumptions being made:
1. The anatomicals are brain masked.
1. The labeled segmentations follow the [FreeSurfer lookup table](https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT).

Here is an example layout of what the `--segmentation-atlases-dir` flag expects:

```
$ tree JLF-templates
JLF-templates/
├── Template01
│   ├── Segmentation.nii.gz
│   ├── T1w_brain.nii.gz
│   └── T2w_brain.nii.gz
└── Template02
    ├── Segmentation.nii.gz
    ├── T1w_brain.nii.gz
    └── T2w_brain.nii.gz
```

## More context on releases
Like other *NiPreps*, *NiBabies* follows Calendar Versioning ([CalVer](https://calver.org/)), in format of `YY.MINOR.MICRO`.
In short, here is a quick heuristic on how new releases should be looked at:
1. If the `YY` or `MINOR` has changed, it is a feature release, with substantial changes to the workflow.
1. If the `YY.MINOR` match the version you used, but the `MICRO` has changed, it is a bug-fix release.
Check the [release notes](https://github.com/nipreps/nibabies/releases) - if the fixes do not pertain to your data, there is no need to upgrade.

For more in-depth information, refer to the [*NiPreps* release documentation](https://www.nipreps.org/devs/releases/#principles).