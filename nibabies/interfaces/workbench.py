from nipype.interfaces.base import CommandLineInputSpec, File, traits, TraitedSpec
from nipype.interfaces.workbench.base import WBCommand


class CiftiDilateInputSpec(CommandLineInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        argstr="%s",
        position=0,
        desc="The input CIFTI file",
    )
    direction = traits.Enum(
        "ROW",
        "COLUMN",
        mandatory=True,
        argstr="%s",
        position=1,
        desc="Which dimension to dilate along, ROW or COLUMN",
    )
    surface_distance = traits.Int(
        mandatory=True,
        argstr="%d",
        position=2,
        desc="The distance to dilate on surfaces, in mm",
    )
    volume_distance = traits.Int(
        mandatory=True,
        argstr="%d",
        position=3,
        desc="The distance to dilate in the volume, in mm",
    )
    out_file = File(
        name_source=["in_file"],
        name_template="dilated_%s.nii",
        keep_extension=True,
        argstr="%s",
        position=4,
        desc="The dilated CIFTI file",
    )
    left_surface = File(
        exists=True,
        position=5,
        argstr="-left-surface %s",
        desc="Specify the left surface to use",
    )
    left_corrected_areas = File(
        exists=True,
        position=6,
        requires=["left_surface"],
        argstr="-left-corrected-areas %s",
        desc="vertex areas (as a metric) to use instead of computing them from the left surface.",
    )
    right_surface = File(
        exists=True,
        position=7,
        argstr="-right-surface %s",
        desc="Specify the right surface to use",
    )
    right_corrected_areas = File(
        exists=True,
        position=8,
        requires=["right_surface"],
        argstr="-right-corrected-areas %s",
        desc="vertex areas (as a metric) to use instead of computing them from the right surface",
    )
    cerebellum_surface = File(
        exists=True,
        position=9,
        argstr="-cerebellum-surface %s",
        desc="specify the cerebellum surface to use",
    )
    cerebellum_corrected_areas = File(
        exists=True,
        position=10,
        requires=["cerebellum_surface"],
        argstr="-cerebellum-corrected-areas %s",
        desc="vertex areas (as a metric) to use instead of computing them from the cerebellum "
        "surface",
    )
    bad_brainordinate_roi = File(
        exists=True,
        position=11,
        argstr="-bad-brainordinate-roi %s",
        desc="CIFTI dscalar or dtseries file, positive values denote brainordinates to have their "
        "values replaced",
    )
    nearest = traits.Bool(
        position=12,
        argstr="-nearest",
        desc="Use nearest good value instead of a weighted average",
    )
    merged_volume = traits.Bool(
        position=13,
        argstr="-merged-volume",
        desc="treat volume components as if they were a single component",
    )
    legacy_mode = traits.Bool(
        position=14,
        argstr="-legacy-mode",
        desc="Use the math from v1.3.2 and earlier for weighted dilation",
    )


class CiftiDilateOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="Dilated CIFTI file")


class CiftiDilate(WBCommand):
    """
    Dilate a CIFTI file.

    For all data values designated as bad, if they neighbor a good value or
    are within the specified distance of a good value in the same kind of
    model, replace the value with a distance weighted average of nearby good
    values, otherwise set the value to zero.  If -nearest is specified, it
    will use the value from the closest good value within range instead of a
    weighted average.  When the input file contains label data, nearest
    dilation is used on the surface, and weighted popularity is used in the
    volume.

    The -*-corrected-areas options are intended for dilating on group average
    surfaces, but it is only an approximate correction for the reduction of
    structure in a group average surface.

    If -bad-brainordinate-roi is specified, all values, including those with
    value zero, are good, except for locations with a positive value in the
    ROI.  If it is not specified, only values equal to zero are bad.

    """

    input_spec = CiftiDilateInputSpec
    output_spec = CiftiDilateOutputSpec
    _cmd = "wb_command -cifti-dilate"


class VolumeAffineResampleInputSpec(CommandLineInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        argstr="%s",
        position=0,
        desc="volume to resample",
    )
    volume_space = File(
        exists=True,
        mandatory=True,
        argstr="%s",
        position=1,
        desc="a volume file in the volume space you want for the output",
    )
    method = traits.Enum(
        "CUBIC", "ENCLOSING_VOXEL", "TRILINEAR",
        mandatory=True,
        argstr="%s",
        position=2,
        desc="The resampling method. The recommended methods are CUBIC "
             "(cubic spline) for most data, and ENCLOSING_VOXEL for label data.",
    )
    out_file = File(
        name_source=["in_file"],
        name_template="resampled_%s.nii.gz",
        keep_extension=True,
        argstr="%s",
        position=3,
        desc="the output volume",
    )
    affine = File(
        exists=True,
        mandatory=True,
        argstr="-affine %s",
        position=4,
        desc="the affine file to apply",
    )
    flirt = traits.Bool(
        argstr="-flirt %s %s",
        position=5,
        desc="Set ``True`` if ``affine`` is a FLIRT affine.",
    )
    flirt_source_volume = File(
        exists=True,
        desc="the source volume used when generating the affine; defaults to in_file",
        requires=['flirt'],
    )
    flirt_target_volume = File(
        exists=True,
        desc="the target volume used when generating the affine; defaults to volume_space",
        requires=['flirt'],
    )


class VolumeAffineResampleOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="the output volume")


class VolumeAffineResample(WBCommand):
    """
    Resample a volume file with an affine transformation.

    >>> from nibabies.interfaces.workbench import VolumeAffineResample
    >>> resample = VolumeAffineResample()
    >>> resample.inputs.in_file = data_dir /'functional.nii'
    >>> resample.inputs.volume_space = data_dir /'anatomical.nii'
    >>> resample.inputs.method = 'CUBIC'
    >>> resample.inputs.affine = data_dir / 'func_to_struct.mat'
    >>> resample.cmdline  #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    'wb_command -volume-resample .../functional.nii .../anatomical.nii CUBIC \
    resampled_functional.nii.gz -affine .../func_to_struct.mat'

    If the affine was generated with FLIRT, this should be indicated.
    By default, the interface will use the ``in_file`` and ``volume_space``
    for references.

    >>> resample.inputs.flirt = True
    >>> resample.cmdline  #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    'wb_command -volume-resample .../functional.nii .../anatomical.nii CUBIC \
    resampled_functional.nii.gz -affine .../func_to_struct.mat \
    -flirt .../functional.nii .../anatomical.nii'

    However, if other volumes were used to calculate the affine, they can
    be provided:

    >>> resample.inputs.flirt_source_volume = data_dir / 'epi.nii'
    >>> resample.inputs.flirt_target_volume = data_dir /'T1w.nii'
    >>> resample.cmdline  #doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
    'wb_command -volume-resample .../functional.nii .../anatomical.nii CUBIC \
    resampled_functional.nii.gz -affine .../func_to_struct.mat \
    -flirt .../epi.nii .../T1w.nii'
    """

    input_spec = VolumeAffineResampleInputSpec
    output_spec = VolumeAffineResampleOutputSpec
    _cmd = "wb_command -volume-resample"

    def _format_arg(self, opt, spec, val):
        if opt == "flirt" and val:
            val = (
                self.inputs.flirt_source_volume or self.inputs.in_file,
                self.inputs.flirt_target_volume or self.inputs.volume_space,
            )
        return super()._format_arg(opt, spec, val)
