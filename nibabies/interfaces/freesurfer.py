import os
from pathlib import Path

from nipype.interfaces.base import (
    traits, File, Directory, CommandLine,
    isdefined, CommandLineInputSpec, TraitedSpec
)


class InfantReconAllInputSpec(CommandLineInputSpec):
    subjects_dir = Directory(
        exists=True,
        hash_files=False,
        desc="path to subjects directory",
    )
    subject_id = traits.Str(
        "recon_all", argstr="--subject %s", desc="subject name", required=True,
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
    aseg_file = File(
        argstr='--aseg %s',
        desc="Pre-computed segmentation file",
    )


class InfantReconAllOutputSpec(TraitedSpec):
    outdir = Directory(exists=True, desc="Output directory.")
    subject_id = traits.Str(desc="Subject name for whom to retrieve data")


class InfantReconAll(CommandLine):
    """
    Runs the infant recon all pipeline
    """

    _cmd = 'infant_recon_all'
    input_spec = InfantReconAllInputSpec
    output_spec = InfantReconAllOutputSpec

    def _run_interface(self, runtime):
        # make sure directory structure is intact
        if not isdefined(self.inputs.subjects_dir):
            self.inputs.subjects_dir = _set_subjects_dir()
        subjdir = Path(self.inputs.subjects_dir) / self.inputs.subject_id
        subjdir.mkdir(parents=True, exist_ok=True)
        # T1 image is expected to be in a specific location if no mask is present
        if not (subjdir / 'mprage.nii.gz').exists() and not (subjdir / 'mprage.mgz').exists():
            if isdefined(self.inputs.t1_file):
                Path(self.inputs.t1_file).symlink_to(subjdir / 'mprage.nii.gz')
            elif not isdefined(self.inputs.mask_file):
                raise RuntimeError("Neither T1 or mask present!")

        return super()._run_interface(runtime)

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['subject_id'] = self.inputs.subject_id
        if isdefined(self.inputs.outdir):
            outputs["outdir"] = self.inputs.outdir
        else:
            outputs["outdir"] = str(Path(self.inputs.subjects_dir) / self.inputs.subject_id)
        return outputs


def _set_subjects_dir():
    subjdir = os.getenv('SUBJECTS_DIR')
    if not subjdir:
        subjdir = os.getcwd()
        os.environ['SUBJECTS_DIR'] = subjdir
    return subjdir
