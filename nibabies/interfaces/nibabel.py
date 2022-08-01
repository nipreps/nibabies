from pathlib import Path

from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)


class ReorientImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="Moving file")
    target_file = File(
        exists=True, xor=["target_orientation"], desc="Reference file to reorient to"
    )
    target_orientation = traits.Str(
        xor=["target_file"], desc="Axis codes of coordinate system to reorient to"
    )


class ReorientImageOutputSpec(TraitedSpec):
    out_file = File(desc="Reoriented file")


class ReorientImage(SimpleInterface):
    input_spec = ReorientImageInputSpec
    output_spec = ReorientImageOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = reorient_image(
            self.inputs.in_file,
            target_file=self.inputs.target_file,
            target_ornt=self.inputs.target_orientation,
        )
        return runtime


def reorient_image(
    in_file: str, *, target_file: str = None, target_ornt: str = None, newpath: str = None
) -> str:
    """
    Reorient an image.

    New orientation targets can be either another image, or a string representation of the
    orientation axis.

    Parameters
    ----------
    in_file : Image to be reoriented
    target_file : Reference image of desired orientation
    target_ornt : Orientation denoted by the first letter of each axis (i.e., "RAS", "LPI")
    """
    import nibabel as nb

    img = nb.load(in_file)
    img_axcodes = nb.aff2axcodes(img.affine)
    in_ornt = nb.orientations.axcodes2ornt(img_axcodes)

    if target_file:
        target_img = nb.load(target_file)
        target_ornt = nb.aff2axcodes(target_img.affine)

    out_ornt = nb.orientations.axcodes2ornt(target_ornt)
    ornt_xfm = nb.orientations.ornt_transform(in_ornt, out_ornt)
    reoriented = img.as_reoriented(ornt_xfm)

    if newpath is None:
        newpath = Path()
    out_file = str((Path(newpath) / "reoriented.nii.gz").absolute())
    reoriented.to_filename(out_file)
    return out_file
