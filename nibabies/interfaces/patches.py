from pathlib import Path

from nipype.interfaces import (
    freesurfer as fs,
)
from nipype.interfaces.ants.base import ANTSCommand, ANTSCommandInputSpec
from nipype.interfaces.ants.registration import (
    CompositeTransformUtil as _CompositeTransformUtil,
)
from nipype.interfaces.ants.registration import (
    CompositeTransformUtilOutputSpec as _CompositeTransformUtilOutputSpec,
)
from nipype.interfaces.base import File, InputMultiObject, TraitedSpec, traits


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


class ConcatXFMInputSpec(ANTSCommandInputSpec):
    transforms = InputMultiObject(
        traits.Either(File(exists=True), 'identity'),
        argstr='%s',
        mandatory=True,
        desc='transform files: will be applied in reverse order. For '
        'example, the last specified transform will be applied first.',
    )
    out_xfm = traits.File(
        'concat_xfm.h5',
        usedefault=True,
        argstr='--output [ %s, 1 ]',
        desc='output file name',
    )
    reference_image = File(
        argstr='--reference-image %s',
        mandatory=True,
        desc='reference image space that you wish to warp INTO',
        exists=True,
    )
    invert_transform_flags = InputMultiObject(traits.Bool())


class ConcatXFMOutputSpec(TraitedSpec):
    out_xfm = File(desc='Combined transform')


class ConcatXFM(ANTSCommand):
    """
    Streamed use of antsApplyTransforms to combine multiple xfms into a single file

    Examples
    --------

    >>> from nibabies.interfaces.patches import ConcatXFM
    >>> cxfm = ConcatXFM()
    >>> cxfm.inputs.transforms = [testdir / 'xfm0.h5', testdir / 'xfm1.h5']
    >>> cxfm.inputs.reference_image = testdir / 'anatomical.nii'
    >>> cxfm.cmdline  # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    'antsApplyTransforms --output [ concat_xfm.h5, 1 ] --reference-image .../anatomical.nii \
--transform .../xfm0.h5 --transform .../xfm1.h5'

    """

    _cmd = 'antsApplyTransforms'
    input_spec = ConcatXFMInputSpec
    output_spec = ConcatXFMOutputSpec

    def _get_transform_filenames(self):
        retval = []
        invert_flags = self.inputs.invert_transform_flags
        if not invert_flags:
            invert_flags = [False] * len(self.inputs.transforms)
        elif len(self.inputs.transforms) != len(invert_flags):
            raise ValueError(
                'ERROR: The invert_transform_flags list must have the same number '
                'of entries as the transforms list.'
            )

        for transform, invert in zip(self.inputs.transforms, invert_flags, strict=False):
            if invert:
                retval.append(f'--transform [ {transform}, 1 ]')
            else:
                retval.append(f'--transform {transform}')
        return ' '.join(retval)

    def _format_arg(self, opt, spec, val):
        if opt == 'transforms':
            return self._get_transform_filenames()
        return super()._format_arg(opt, spec, val)

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_xfm'] = Path(self.inputs.out_xfm).absolute()
        return outputs


class CompositeTransformUtilOutputSpec(_CompositeTransformUtilOutputSpec):
    out_transforms = traits.List(desc='list of ordered transform components')


class CompositeTransformUtil(_CompositeTransformUtil):
    """Outputs have changed in newer versions of ANTs."""

    output_spec = CompositeTransformUtilOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()

        # Ordering may change depending on forward/inverse transform
        # Forward: <prefix>_00_AffineTransform.mat, <prefix>_01_DisplacementFieldTransform.nii.gz
        # Inverse: <prefix>_01_AffineTransform.mat, <prefix>_00_DisplacementFieldTransform.nii.gz
        if self.inputs.process == 'disassemble':
            transforms = [
                str(Path(x).absolute())
                for x in sorted(Path().glob(f'{self.inputs.output_prefix}_*'))
            ]
            outputs['out_transforms'] = transforms

            # Potentially could be more than one affine / displacement per composite transform...
            outputs['affine_transform'] = [
                x for x in transforms if 'AffineTransform' in Path(x).name
            ][0]
            outputs['displacement_field'] = [
                x for x in transforms if 'DisplacementFieldTransform' in Path(x).name
            ][0]
        elif self.inputs.process == 'assemble':
            outputs['out_file'] = Path(self.inputs.out_file).absolute()
        return outputs
