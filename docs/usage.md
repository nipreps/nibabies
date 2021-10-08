# Usage

## The BIDS format

The *NiBabies* workflow takes as principal input the path of the dataset
that is to be processed.
The input dataset is required to be in valid
{abbr}`BIDS (The Brain Imaging Data Structure)` format,
and it must include at least one T1-weighted and 
one T2-weighted structural image and
(unless disabled with a flag) a BOLD series.
We highly recommend that you validate your dataset with the free, online
[BIDS Validator](http://bids-standard.github.io/bids-validator/).

The exact command to run *NiBabies* depends on the [Installation] method.
The common parts of the command follow the
[BIDS-Apps](https://github.com/BIDS-Apps) definition.
Example:

```Shell
fmriprep data/bids_root/ out/ participant -w work/
```

Further information about BIDS and BIDS-Apps can be found at the
[NiPreps portal](https://www.nipreps.org/apps/framework/).

## Command-Line Arguments
```{argparse}
:ref: nibabies.cli.parser._build_parser
:prog: nibabies
:nodefault:
:nodefaultconst:
```

The command-line interface of the docker wrapper
------------------------------------------------

```{argparse}
:ref: nibabies_wrapper.get_parser
:prog: nibabies-wrapper
:nodefault:
:nodefaultconst:
```



