import nipype.interfaces.freesurfer as fs
from nipype.interfaces.base import File, traits


class _MRICoregInputSpec(fs.registration.MRICoregInputSpec):
    reference_file = File(
        argstr='--ref %s',
        desc='reference (target) file',
        copyfile=False,
    )
    subject_id = traits.Str(
        argstr='--s %s',
        position=1,
        requires=['subjects_dir'],
        desc='freesurfer subject ID (implies ``reference_mask == '
        'aparc+aseg.mgz`` unless otherwise specified)',
    )


class MRICoreg(fs.MRICoreg):
    """
    Patched that allows setting both a reference file and the subjects dir.
    """

    input_spec = _MRICoregInputSpec
