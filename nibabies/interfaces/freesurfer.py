import os
from pathlib import Path

from nipype.interfaces.freesurfer.base import FSCommand
from nipype.interfaces.freesurfer.preprocess import ReconAllOutputSpec
from nipype.interfaces.base import (
    traits, File, Directory,
    isdefined, CommandLineInputSpec, TraitedSpec
)

class InfantReconAllInputSpec(CommandLineInputSpec):
    subjects_dir = Directory(
        exists=True,
        hash_files=False,
        desc="path to subjects directory",
        genfile=True,
    )
    subject_id = traits.Str(
        "recon_all", argstr="-subjid %s", desc="subject name", usedefault=True
    )
    t1_file = File(
        exists=True,
        desc="path to T1w file",
    )
    age = traits.Range(
        low=0, high=24, argstr='--age %d', required=True, desc="Subject age in months"
    )
    outdir = Directory(
        argstr='--outdir %s',
        desc="Output directory where the reconall results are written."
             "The default location is <subjects_dir>/<subject_id>",
    )
    mask_file = traits.File(
        argstr='--masked %s',
        desc="Skull-stripped and INU-corrected T1 (skips skullstripping step)",
    )
    force = traits.Bool(
        argstr='--force',
        desc="Force all the processing to be (re)done",
    )

class InfantReconAllOutputSpec(TraitedSpec):
    outdir = Directory(exists=True, desc="Output directory")

class InfantReconAll(FSCommand):
    """
    Runs the infant recon all pipeline
    """

    _cmd = 'infant_recon_all'
    input_spec = InfantReconAllInputSpec
    output_spec = ReconAllOutputSpec

    def _parse_inputs(self):
        if not isdefined(self.inputs.subjects_dir):
            self.inputs.subjects_dir = _get_subjects_dir()
        super()._parse_inputs()

    def _run_interface(self, runtime):
        # make sure directory structure is intact
        subjdir = Path(self.inputs.subjects_dir) / self.inputs.subject_id
        subjdir.mkdir(parents=True, exist_ok=True)
        # T1 image is expected to be in a specific location
        if not (subjdir / 'mprage.nii.gz').exists() and not (subjdir / 'mprage.mgz').exists():
            if not isdefined(self.inputs.t1_file):
                raise Exception("T1 is required!")
            Path(self.inputs.t1_file).symlink_to(subjdir / 'mprage.nii.gz')
        return super()._run_interface(runtime)

    def _list_outputs(self):
        outputs = self._outputs().get()
        if isdefined(self.inputs.outdir):
            outputs["outdir"] = self.inputs.outdir
        else:
            outputs["outdir"] = str(Path(self.inputs.subjects_dir) / self.inputs.subject_id)
        return outputs


def _get_subjects_dir():
    return os.getenv('SUBJECTS_DIR') or os.getcwd()
