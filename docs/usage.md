# Usage

## The BIDS format

The *NiBabies* workflow takes as principal input the path of the dataset
that is to be processed.
The input dataset is required to be in valid
{abbr}`BIDS (The Brain Imaging Data Structure)` format,
and it must include at least one T1-weighted and
one T2-weighted structural image and
a BOLD series (unless using the `--anat-only` flag).

We highly recommend that you validate your dataset with the free, online
[BIDS Validator](http://bids-standard.github.io/bids-validator/).

### Participant Ages
*NiBabies* will attempt to automatically extract participant ages (in months) from the BIDS layout.
Specifically, these two files will be checked:
- [Sessions file](https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html#sessions-file): `<bids-root>/<subject>/subject_sessions.tsv`
- [Participants file](https://bids-specification.readthedocs.io/en/stable/03-modality-agnostic-files.html#participants-file): `<bids-root>/participants.tsv`

Either file should include `age` (or if you wish to be more explicit: `age_months`) columns, and it is
recommended to have an accompanying JSON file to further describe these fields, and explicitly state the values are in months.

Age is used to choose the optimal surface reconstruction method when the flag `--surface-recon-method auto` is employed for surface reconstruction:

 - M-CRIB-S: 0-3 months old (only if T2w and segmentation in T2w space is present)
 - infant FreeSurfer: 4-24 months old (or 0-24 months old if no T2w/T2w segmentation is present)
 - FreeSurfer (standard FreeSurfer developed from adult data): >24 months old

## The FreeSurfer license

*NiBabies* uses FreeSurfer tools, which require a license to run.

To obtain a FreeSurfer license, simply register for free at https://surfer.nmr.mgh.harvard.edu/registration.html.

FreeSurfer will search for a license key file first using the `$FS_LICENSE` environment variable and then in the default path to the license key file (`$FREESURFER_HOME`/license.txt). If `$FS_LICENSE` is set, the [`nibabies-wrapper`](#using-the-nibabies-wrapper) will automatically handle setting the license within the container.
Otherwise, you will need to use the `--fs-license-file` flag to ensure the license is available.


## Example command

The exact command to run *NiBabies* depends on the [Installation](./installation.md) method.
The common parts of the command follow the
[BIDS-Apps](https://github.com/BIDS-Apps) definition.
Example:

```Shell
$ nibabies data/bids_root/ out/ participant -w work/ --participant-label 01
```

Further information about BIDS and BIDS-Apps can be found at the
[NiPreps portal](https://www.nipreps.org/apps/framework/).

## Command-Line Arguments
```{argparse}
:ref: nibabies.cli.parser._build_parser
:prog: nibabies
:nodefaultconst:
```

## More information on command-line arguments

At minimum, the following *positional* arguments are required.

- **`bids_dir`** - the root folder of a BIDS valid dataset.
- **`output_dir`** - folder to store outputs and reports.
- **`level`** - processing stage to be run, currently can only be `participant`.

However, as infant brains can vastly differ depending on age, providing the following arguments is highly recommended:

- **`--participant-label`** - participant ID

- **`--session-id`** - session ID

- **`--segmentation-atlases-dir`** - directory containing pre-labeled segmentations to use for Joint Label Fusion.

:::{admonition} Tip
:class: tip

The segmentation directory layout should consist of one or more template directories containing:
* A segmented and labeled NIfTI that includes `Segmentation` in the filename.
* A brainmasked T1w NIfTI that includes `T1w` in the filename.

:::

## Using the nibabies wrapper

The wrapper will generate a Docker or Singularity command line for you, print it out for reporting purposes, and then execute it without further action needed.
For installation instructions, please see [](installation.md#installing-the-nibabies-wrapper)

### Sample Docker usage

```
$ nibabies-wrapper docker /path/to/data /path/to/output participant --fs-license-file /usr/freesurfer/license.txt

RUNNING: docker run --rm -e DOCKER_VERSION_8395080871=20.10.6 -it -v /path/to/data:/data:ro \
-v /path/to/output:/out -v /usr/freesurfer/license.txt:/opt/freesurfer/license.txt:ro \
nipreps/nibabies:23.0.0 /data /out participant
...
```

:::{admonition} Docker usage warning
:class: warning

When using Docker, the wrapper will default to using the same version of `nibabies` as the wrapper.
This can be overridden by using the `-i` flag to specify a particular Docker image.
:::

### Sample Singularity usage

```
$ nibabies-wrapper singularity /path/to/data /path/to/output participant -i nibabies-23.0.0.sif --fs-license-file /usr/freesurfer/license.txt

RUNNING: singularity run --cleanenv -B /path/to/data:/data:ro \
-B /path/to/output:/out -B /usr/freesurfer/license.txt:/opt/freesurfer/license.txt:ro \
nibabies-23.0.0.sif /data /out participant
...
```

:::{admonition} Singularity usage warning
:class: warning

Note that the `-i` flag is required when using Singularity, and should be the path to the already built Singularity image file.
:::

The command-line interface of the nibabies wrapper
------------------------------------------------

```{argparse}
:ref: nibabies_wrapper.get_parser
:prog: nibabies-wrapper
:nodefault:
:nodefaultconst:
```
