# Tips and FAQs

## Leveraging precomputed results

Whether manual intervention is required, or you want to break up processing, *NiBabies* can reuse previously-computed files (either from NiBabies directly or a third party application) to be injected into the workflow directly.

:::{versionchanged} 24.0.0

In addition to the brain mask and anatomical segmentation, support was added for additional precomputed derivatives. To see which derivatives are supported, view [](outputs.md#anatomical-derivatives).
:::

To use pre-computed results, one or more [BIDS Derivatives](https://bids-specification.readthedocs.io/en/stable/05-derivatives/01-introduction.html#bids-derivatives) directories must be passed in to *NiBabies* using the `--derivatives` flag.
Derivative directories must include a [`dataset_description.json` and the required fields](https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html#derived-dataset-and-pipeline-description).
Additionally, files must include the `space-T1w` or `space-T2w` key-value pair in the filenames, and a matching sidecar JSON file with the `SpatialReference` field defined.

A sample layout of a derivatives directory can be found below:

```bash
my_precomputed/
├── dataset_description.json
└── sub-01
    └── anat
        ├── sub-01_desc-preproc_T2w.nii.gz
        ├── sub-01_space-T2w_desc-aseg_dseg.json
        ├── sub-01_space-T2w_desc-aseg_dseg.nii.gz
        ├── sub-01_space-T2w_desc-brain_mask.json
        └── sub-01_space-T2w_desc-brain_mask.nii.gz
```

In this example, `sub-01_desc-preproc_T2w.nii.gz` will be used as the T2w reference. The other files (the brain mask and segmentation), will be in the same space.

:::{warning}
If no anatomical reference is provided, the outputs must be in the same space as the raw anatomical data.
:::

:::{note}
If an aseg is provided, it will be used for surface generation.
:::

## Multi-atlas segmentation with joint label fusion

By default, *NiBabies* will run [FSL FAST](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FAST) for tissue segmentation, and Infant FreeSurfer for segmentation labels.

Alternatively, ANTs {abbr}`JLF (Joint Label Fusion)` can be used by providing a directory with one or more template images composed of anatomicals and segmentations. To pass in this directory, use the `--segmentation-atlases-dir` flag.
When using this approach, there are a few assumptions being made:

1. The anatomicals are brain masked.
1. The segmentation labels adhere to the [FreeSurfer LUT](https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/AnatomicalROI/FreeSurferColorLUT).

Here is an example layout of what the `--segmentation-atlases-dir` flag expects:

```bash
$ tree JLF-atlases

JLF-atlases/
├── dataset_description.json
├── participants.tsv
├── sub-01
│   ├── sub-01_desc-aseg_dseg.nii.gz
│   ├── [sub-01_T1w.json]  * optional
│   ├── sub-01_T1w.nii.gz
│   ├── [sub-01_T2w.json]  * optional
│   └── sub-01_T2w.nii.gz
├── sub-02
...
```

## More context on releases

Like other *NiPreps*, *NiBabies* follows Calendar Versioning ([CalVer](https://calver.org/)), in format of `YY.MINOR.MICRO`.
In short, here is a quick heuristic on how new releases should be looked at:

1. If the `YY` or `MINOR` has changed, it is a feature release, with substantial changes to the workflow.
1. If the `YY.MINOR` match the version you used, but the `MICRO` has changed, it is a bug-fix release.
Check the [release notes](https://github.com/nipreps/nibabies/releases) - if the fixes do not pertain to your data, there is no need to upgrade.

For more in-depth information, refer to the [*NiPreps* release documentation](https://www.nipreps.org/devs/releases/#principles).
