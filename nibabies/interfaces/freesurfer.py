import os
import logging
from pathlib import Path

from nipype.interfaces.base import (
    traits, File, Directory, CommandLine,
    isdefined, CommandLineInputSpec, TraitedSpec
)

from ..utils.misc import check_total_memory

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
        low=0,
        high=24,
        argstr='--age %d',
        desc="Subject age in months",
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
    newborn = traits.Bool(
        xor=['age'],
        argstr='--newborn',
        help="Use newborns from set",
    )
    aseg_file = File(
        argstr='--segfile',
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
        if not isdefined(self.inputs.outdir):
            subjdir = Path(self.inputs.subjects_dir) / self.inputs.subject_id
            try:
                subjdir.mkdir(parents=True, exist_ok=True)
            except OSError:
                raise OSError(
                    f"Current SUBJECTS_DIR <{subjdir}> cannot be written to. To fix this,"
                    "either define the input or unset the environmental variable."
                )
            # T1 image is expected to be in a specific location if no mask is present
            if not (subjdir / 'mprage.nii.gz').exists() and not (subjdir / 'mprage.mgz').exists():
                if isdefined(self.inputs.t1_file):
                    Path(self.inputs.t1_file).symlink_to(subjdir / 'mprage.nii.gz')
                elif not isdefined(self.inputs.mask_file):
                    raise RuntimeError("Neither T1 or mask present!")
        if isdefined(self.inputs.aseg_file):
            pass  # To be added in a future infant-FS release.

        # warn users that this might fail...
        if not check_total_memory(recommended_gb=20):
            import logging
            logging.getLogger('nipype.interface').warning(
                f"For best results, run {self._cmd} with at least 20GB available RAM."
            )

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
