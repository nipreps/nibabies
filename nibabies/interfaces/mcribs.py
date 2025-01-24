import shutil
from pathlib import Path

from nipype.interfaces.base import (
    CommandLine,
    CommandLineInputSpec,
    Directory,
    File,
    TraitedSpec,
    traits,
)


class MCRIBReconAllInputSpec(CommandLineInputSpec):
    # Input structure massaging
    outdir = Directory(
        exists=True,
        hash_files=False,
        desc='Path to save output, or path of existing MCRIBS output',
    )
    subjects_dir = Directory(
        exists=True,
        hash_files=False,
        desc='Path to FreeSurfer subjects directory',
    )
    subject_id = traits.Str(
        required=True,
        argstr='%s',
        position=-1,
        desc='Subject ID',
    )
    t1w_file = File(
        exists=True,
        copyfile=True,
        desc='T1w to be used for deformable (must be in register with T2w image)',
    )
    t2w_file = File(
        exists=True,
        copyfile=True,
        desc='T2w (Isotropic + N4 corrected)',
    )
    segmentation_file = File(
        exists=True,
        desc='Segmentation file (skips tissue segmentation)',
    )
    mask_file = File(exists=True, desc='T2w mask')

    # MCRIBS options
    conform = traits.Bool(
        False,
        usedefault=True,
        argstr='--conform',
        desc='Reorients to radiological, axial slice orientation. Resamples to isotropic voxels',
    )
    tissueseg = traits.Bool(
        False,
        usedefault=True,
        argstr='--tissueseg',
        desc='Perform tissue type segmentation',
    )
    surfrecon = traits.Bool(
        False,
        usedefault=True,
        argstr='--surfrecon',
        desc='Reconstruct surfaces',
    )
    surfrecon_method = traits.Enum(
        'Deformable',
        argstr='--surfreconmethod %s',
        requires=['surfrecon'],
        desc='Surface reconstruction method',
    )
    join_thresh = traits.Float(
        argstr='--deformablejointhresh %f',
        requires=['surfrecon'],
        desc='Join threshold parameter for Deformable',
    )
    fast_collision = traits.Bool(
        argstr='--deformablefastcollision',
        requires=['surfrecon'],
        desc='Use Deformable fast collision test',
    )
    autorecon_after_surf = traits.Bool(
        False,
        usedefault=True,
        argstr='--autoreconaftersurf',
        desc='Do all steps after surface reconstruction',
    )
    segstats = traits.Bool(
        argstr='--segstats',
        desc='Compute statistics on segmented volumes',
    )
    nthreads = traits.Int(
        argstr='-nthreads %d',
        desc='Number of threads for multithreading applications',
    )


class MCRIBReconAllOutputSpec(TraitedSpec):
    mcribs_dir = Directory(desc='MCRIBS output directory')
    subjects_dir = Directory(desc='FreeSurfer output directory')


class MCRIBReconAll(CommandLine):
    _cmd = 'MCRIBReconAll'
    input_spec = MCRIBReconAllInputSpec
    output_spec = MCRIBReconAllOutputSpec
    _no_run = False

    _expected_files = {
        'surfrecon': {
            'meshes': (
                'pial-lh-reordered.vtp',
                'pial-rh-reordered.vtp',
                'white-rh.vtp',
                'white-lh.vtp',
            )
        },
        'autorecon': {
            'mri': ('T2.mgz', 'aseg.presurf.mgz', 'ribbon.mgz', 'brain.mgz'),
            'label': ('lh.cortex.label', 'rh.cortex.label'),
            'stats': ('aseg.stats', 'brainvol.stats', 'lh.aparc.stats', 'rh.curv.stats'),
            'surf': (
                'lh.pial', 'rh.pial',
                'lh.white', 'rh.white',
                'lh.curv', 'rh.curv',
                'lh.thickness', 'rh.thickness'),
        }
    }  # fmt:skip

    @property
    def cmdline(self):
        cmd = super().cmdline
        if self.inputs.surfrecon and not self.inputs.t2w_file:
            raise AttributeError('Missing `t2w_file` input.')

        # If an output directory is provided, check if we can skip the run
        if self.inputs.outdir:
            sid = self.inputs.subject_id
            # Check MIRTK surface recon deformable
            if self.inputs.surfrecon:
                surfrecon_dir = Path(self.inputs.outdir) / sid / 'SurfReconDeformable' / sid
                if self._verify_outputs('surfrecon', surfrecon_dir):
                    self._no_run = True
            # Check FS directory population
            if self.inputs.autorecon_after_surf:
                fs_dir = Path(self.inputs.outdir) / sid / 'freesurfer' / sid
                if self._verify_outputs('autorecon', fs_dir):
                    self._no_run = True
                else:
                    self._no_run = False

            if self._no_run:
                return 'echo MCRIBReconAll: nothing to do'
        return cmd

    def _setup_directory_structure(self, mcribs_dir: Path) -> None:
        """
        Create the required structure for skipping steps.

        The directory tree
        ------------------

        <subject_id>/
        ├── RawT2
        │   └── <subject_id>.nii.gz
        ├── SurfReconDeformable
        │   └── <subject_id>
        │       └── temp
        │           └── t2w-image.nii.gz
        ├── TissueSeg
        │   ├── <subject_id>_all_labels.nii.gz
        │   └── <subject_id>_all_labels_manedit.nii.gz
        └── TissueSegDrawEM
            └── <subject_id>
                └── N4
                    └── <subject_id>.nii.gz
        """
        sid = self.inputs.subject_id
        mkdir_kw = {'parents': True, 'exist_ok': True}
        root = mcribs_dir / sid
        root.mkdir(**mkdir_kw)

        # T2w operations
        if self.inputs.t2w_file:
            t2w = root / 'RawT2' / f'{sid}.nii.gz'
            t2w.parent.mkdir(**mkdir_kw)
            if not t2w.exists():
                shutil.copy(self.inputs.t2w_file, str(t2w))

            if not self.inputs.conform:
                t2wiso = root / 'RawT2RadiologicalIsotropic' / f'{sid}.nii.gz'
                t2wiso.parent.mkdir(**mkdir_kw)
                if not t2wiso.exists():
                    t2wiso.symlink_to(f'../RawT2/{sid}.nii.gz')

                n4 = root / 'TissueSegDrawEM' / sid / 'N4' / f'{sid}.nii.gz'
                n4.parent.mkdir(**mkdir_kw)
                if not n4.exists():
                    n4.symlink_to(f'../../../RawT2/{sid}.nii.gz')

        # Segmentation
        if self.inputs.segmentation_file:
            # TissueSeg directive disabled
            tisseg = root / 'TissueSeg' / f'{sid}_all_labels.nii.gz'
            tisseg.parent.mkdir(**mkdir_kw)
            if not tisseg.exists():
                shutil.copy(self.inputs.segmentation_file, str(tisseg))
            manedit = tisseg.parent / f'{sid}_all_labels_manedit.nii.gz'
            if not manedit.exists():
                manedit.symlink_to(tisseg.name)

            if self.inputs.surfrecon:
                t2wseg = root / 'TissueSeg' / f'{sid}_t2w_restore.nii.gz'
                if not t2wseg.exists():
                    t2wseg.symlink_to(f'../RawT2/{sid}.nii.gz')

                surfrec = root / 'SurfReconDeformable' / sid / 'temp' / 't2w-image.nii.gz'
                surfrec.parent.mkdir(**mkdir_kw)
                if not surfrec.exists():
                    surfrec.symlink_to(f'../../../RawT2/{sid}.nii.gz')

                if self.inputs.mask_file:
                    surfrec_mask = surfrec.parent / 'brain-mask.nii.gz'
                    if not surfrec_mask.exists():
                        shutil.copy(self.inputs.mask_file, str(surfrec_mask))

        if self.inputs.surfrecon:
            # Create FreeSurfer layout to safeguard against cd-ing into missing directories
            for d in ('surf', 'mri', 'label', 'scripts', 'stats'):
                (root / 'freesurfer' / sid / d).mkdir(**mkdir_kw)

        # TODO?: T1w -> <subject_id>/RawT1RadiologicalIsotropic/<subjectid>.nii.gz
        return

    def _run_interface(self, runtime):
        # if users wish to preserve their runs
        mcribs_dir = self.inputs.outdir or Path(runtime.cwd) / 'mcribs'
        self._mcribs_dir = Path(mcribs_dir)
        if self.inputs.surfrecon:
            self._setup_directory_structure(self._mcribs_dir)
        # overwrite CWD to be in MCRIB subject's directory
        runtime.cwd = str(self._mcribs_dir / self.inputs.subject_id)
        return super()._run_interface(runtime)

    def _list_outputs(self):
        outputs = self._outputs().get()
        sid = self.inputs.subject_id
        if self.inputs.surfrecon:
            # verify surface reconstruction was successful
            surfrecon_dir = self._mcribs_dir / sid / 'SurfReconDeformable' / sid
            self._verify_outputs('surfrecon', surfrecon_dir, error=True)

        mcribs_fs = self._mcribs_dir / sid / 'freesurfer' / sid
        if self.inputs.autorecon_after_surf:
            self._verify_outputs('autorecon', mcribs_fs, error=True)

        outputs['mcribs_dir'] = str(self._mcribs_dir)
        if self.inputs.autorecon_after_surf and self.inputs.subjects_dir:
            dst = Path(self.inputs.subjects_dir) / self.inputs.subject_id
            if not dst.exists():
                shutil.copytree(mcribs_fs, dst)
                # Create a file to denote this SUBJECTS_DIR was derived from MCRIBS
                logfile = self._mcribs_dir / sid / 'logs' / f'{sid}.log'
                shutil.copyfile(logfile, (dst / 'scripts' / 'mcribs.log'))
            # Copy registration sphere to better match recon-all output
            for hemi in 'lr':
                orig = dst / 'surf' / f'{hemi}h.sphere.reg2'
                symbolic = Path(str(orig)[:-1])
                if orig.exists() and not symbolic.exists():
                    shutil.copyfile(orig, symbolic)
            outputs['subjects_dir'] = self.inputs.subjects_dir

        return outputs

    def _verify_outputs(self, step: str, root: Path, error: bool = False) -> bool:
        """
        Method to check to ensure the expected files are present successful.

        MCRIBReconAll does not return a failing exit code if a step failed, which leads
        this interface to be marked as completed without error in such cases.
        """

        expected = self._expected_files[step]
        for d, files in expected.items():
            for fl in files:
                if not (root / d / fl).exists():
                    if error:
                        raise FileNotFoundError(f'{step.capitalize()} missing: {fl}')
                    return False
        return True
