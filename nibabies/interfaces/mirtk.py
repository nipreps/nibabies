from pathlib import Path

from nipype.interfaces.base import (
    CommandLine,
    Directory,
    File,
    Str,
    TraitedSpec,
    traits,
)


class ReconNeonatalCortexInputSpec(TraitedSpec):
    t1w_file = File(exists=True, desc="Input T1w file")
    t2w_file = File(exists=True, desc="Input T2w file")
    mask_file = File(exists=True, desc="Input mask file")
    labels_file = File(exists=True, desc="Input labels file")
    tissues_file = File(exists=True, desc="Input tissues file")

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
        desc='Number of cores to use for multi-threading',
    )


class ReconNeonatalCortexOutputSpec(TraitedSpec):
    config_file = File(desc='Configuration file used')
    output_dir = Directory(desc='Output directory')


class ReconNeonatalCortex(CommandLine):
    """
    Example
    -------

    >>> rnc = ReconNeonatalCortex()
    >>> rnc.inputs.sessions = ["sub1-ses1", "sub1-ses2"]
    >>> rnc.cmdline
    'mirtk recon-neonatal-cortex --sessions sub1-ses1 sub1-ses2'
    """

    input_spec = ReconNeonatalCortexInputSpec
    output_spec = ReconNeonatalCortexOutputSpec
    _cmd = "mirtk recon-neonatal-cortex"

    def _create_config(self):
        """Create a configuration file required to run the command."""
        from configparser import ConfigParser

        # Example config
        # https://github.com/DevelopmentalImagingMCRI/MCRIBS/blob/master/lib/Deformable/recon-neonatal-cortex.cfg

        out_dir = Path('out').absolute()
        out_dir.mkdir(exist_ok=True)
        if not (work_dir := self.inputs.work_dir):
            work_dir = Path('work').absolute()
            work_dir.mkdir(exist_ok=True)

        inputs = {}
        if self.inputs.t1w_file:
            inputs['input_t1w_image'] = self.inputs.t1w_file
        if self.inputs.t2w_file:
            inputs['input_t2w_image'] = self.inputs.t2w_file
        if self.inputs.mask_file:
            inputs['input_brain_mask'] = self.inputs.mask_file
        if self.inputs.labels_file:
            inputs['input_labels_image'] = self.inputs.labels_file
        if self.inputs.tissues_file:
            inputs['input_tissues_image'] = self.inputs.tissues_file

        conf = {
            'recon-neonatal-cortex': {
                # directory setup
                'out_dir': str(out_dir),
                'mesh_dir': str(out_dir / 'meshes'),
                'temp_dir': str(work_dir / 'tmp'),
                'logs_dir': str(out_dir / 'logs'),
                'work_dir': str(work_dir),
                'tissueseg_dir': '%(work_dir)s/TissueSeg',
                'tissuesegmcribs_dir': '%(work_dir)s/TissueSegMCRIBS',
                **inputs,  # input images
                # recon options
                'fill_wm_holes': False,
                'white_matter_labels': '51..82',
                'gray_matter_labels': '5..16,20..39',
                'deep_gray_matter_labels': '1..4,40..47,85..87',
                'lateral_ventricles_labels': '49,50',
                'corpus_callosum_labels': '48',
                'inter_hemisphere_labels': '40..47,85..87',
                'brainstem_labels': '19',
                'cerebellum_labels': '17,18',
                'subcortex_closing': 10,
                'brainstem_closing': 10,
                'cerebellum_closing': 10,
                'regions_mask': '%(out_dir)s/recon/regions.nii.gz',
                'cortical_hull_dmap': '%(out_dir)s/recon/cortical-hull-dmap.nii.gz',
                't1w_image': '%(temp_dir)s/t1w-image.nii.gz',
                't2w_image': '%(temp_dir)s/t2w-image.nii.gz',
                'brain_mask': '%(temp_dir)s/brain-mask.nii.gz',
                'white_matter_mask': '%(temp_dir)s/white-matter-mask.nii.gz',
                'gray_matter_mask': '%(temp_dir)s/gray-matter-mask.nii.gz',
                'deep_gray_matter_mask': '%(temp_dir)s/deep-gray-matter-mask.nii.gz',
                'corpus_callosum_mask': '%(temp_dir)s/corpus-callosum-mask.nii.gz',
                'ventricles_mask': '%(temp_dir)s/ventricles-mask.nii.gz',
                'ventricles_dmap': '%(temp_dir)s/ventricles-dmap.nii.gz',
                'brain_mesh': '%(mesh_dir)s/brain.vtp',
                'bs_cb_mesh': '%(mesh_dir)s/brainstem+cerebellum.vtp',
                'internal_mesh': '%(mesh_dir)s/internal.vtp',
                'cerebrum_mesh': '%(temp_dir)s/cerebrum.vtp',
                'right_cerebrum_mesh': '%(temp_dir)s/cerebrum-rh.vtp',
                'left_cerebrum_mesh': '%(temp_dir)s/cerebrum-lh.vtp',
                'white_mesh': '%(mesh_dir)s/white.vtp',
                'right_white_mesh': '%(mesh_dir)s/white-rh.vtp',
                'left_white_mesh': '%(mesh_dir)s/white-lh.vtp',
                'pial_mesh': '%(mesh_dir)s/pial.vtp',
                'right_pial_mesh': '%(mesh_dir)s/pial-rh.vtp',
                'left_pial_mesh': '%(mesh_dir)s/pial-lh.vtp',
            },
            'recon-neonatal-cortex white_model': {'distance': 1},
            'recon-neonatal-cortex pial_model': {'distance': 1},
        }
        parser = ConfigParser()
        parser.read_dict(conf)
        out_file = Path('recon-neonatal-cortex.ini').absolute()
        with out_file.open('w') as fp:
            parser.write(fp)
        return str(out_file)

    def _run_interface(self, runtime):
        if not self.inputs.config:
            self.inputs.config = self._create_config()
        return runtime
