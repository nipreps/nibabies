from nipype.interfaces.base import (
    CommandLine,
    Directory,
    File,
    Str,
    TraitedSpec,
    traits,
)


class ReconNeonatalCortexInputSpec(TraitedSpec):
    work_dir = Directory(
        argstr="--work-dir %s",
        exists=True,
        hash_files=False,
        desc="Root working directory",
    )
    config = File(
        argstr="--config %s",
        exists=True,
        desc="Optional custom configuration file",
    )
    section = Str(
        "recon-neonatal-cortex",
        argstr="--section %s",
        default=True,
        desc="Configuration section name",
    )
    sessions = traits.Either(
        File(exists=True),
        traits.List(Str),
        argstr="--sessions %s",
        mandatory=True,
        desc="Either list of '{SubjectId}[-{SessionId}]' strings or path of CSV file",
    )
    brain = traits.Bool(
        argstr="--brain",
        desc="Reconstruct surface of brain mask",
    )
    white = traits.Bool(
        argstr="--white",
        desc="Reconstruct white surface",
    )
    cerebrum = traits.Bool(
        argstr="--cerebrum",
        desc="Reconstruct/keep initial white surface",
    )
    pial = traits.Bool(
        argstr='--pial',
        desc='Reconstruct pial surface',
    )
    regions_mask = traits.Bool(
        argstr='--regions-mask',
        desc="Create regions label image, implies --keep-regions-mask",
    )
    hindbrain = traits.Bool(
        argstr='--brainstem-and-cerebellum',
        desc="Reconstruct combined brainstem and cerebellum surface",
    )
    pial_outside_white = traits.Bool(
        argstr='--ensure-pial-is-outside-white-surface',
        desc='Ensure that pial surface is strictly outside the white surface',
    )
    join_internal_mesh = traits.Bool(
        argstr='--join-with-internal-mesh',
        desc='Join final mesh with internal (hemispheres) dividier mesh',
    )
    join_with_hindbrain = traits.Bool(
        argstr='--join-with-brainstem-and-cerebellum',
        desc="Merge cerebrum surface mesh with brainstem and cerebellum surface mesh",
    )
    nocut = traits.Bool(
        argstr='--nocut',
        desc='Save individual (closed) genus-0 surfaces for each hemisphere',
    )
    nocheck = traits.Bool(
        argstr='--nocheck',
        desc='Disable consistency and self-intersection checks of (intermediate) surface meshes',
    )
    keep_t1w = traits.Bool(
        argstr='--keep-t1w-image',
        desc="Keep resampled T1-weighted image even when no -debug option given",
    )
    keep_t2w = traits.Bool(
        argstr='--keep-t2w-image',
        desc="Keep resampled T2-weighted image even when no -debug option given",
    )
    keep_regions_mask = traits.Bool(
        argstr='--keep-regions-mask',
        desc="Keep regions label image even when no -debug option given",
    )
    force = traits.Bool(
        argstr='--force',
        desc='Overwrite existing output files',
    )
    join_tol = traits.Float(
        1,
        argstr='--jointol %g',
        desc='Join tolerance',
    )
    use_fast_collision = traits.Bool(
        argstr='-use-fast-collision',
        desc='Use the fast collision test for white and pial surfaces',
    )
    num_threads = traits.Int(
        argstr='--threads %d',
        desc='No. of cores to use for multi-threading',
    )


class ReconNeonatalCortex(CommandLine):
    """
    Example
    -------

    >>> rnc = ReconNeonatalCortex()
    >>> rnc.inputs.sessions = ["sub1-ses1", "sub1-ses2"]
    >>> rnc.cmdline
    'recon-neonatal-cortex --sessions sub1-ses1 sub1-ses2'
    """

    input_spec = ReconNeonatalCortexInputSpec
    _cmd = "recon-neonatal-cortex"
