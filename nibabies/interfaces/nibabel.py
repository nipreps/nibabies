from pathlib import Path

from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    InputMultiObject,
    SimpleInterface,
    TraitedSpec,
    traits,
)


class ReorientImageInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='Moving file')
    target_file = File(
        exists=True, xor=['target_orientation'], desc='Reference file to reorient to'
    )
    target_orientation = traits.Str(
        xor=['target_file'], desc='Axis codes of coordinate system to reorient to'
    )


class ReorientImageOutputSpec(TraitedSpec):
    out_file = File(desc='Reoriented file')


class ReorientImage(SimpleInterface):
    input_spec = ReorientImageInputSpec
    output_spec = ReorientImageOutputSpec

    def _run_interface(self, runtime):
        self._results['out_file'] = reorient_image(
            self.inputs.in_file,
            target_file=self.inputs.target_file,
            target_ornt=self.inputs.target_orientation,
        )
        return runtime


def reorient_image(
    in_file: str,
    *,
    target_file: str | None = None,
    target_ornt: str | None = None,
    newpath: str | None = None,
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
    out_file = str((Path(newpath) / 'reoriented.nii.gz').absolute())
    reoriented.to_filename(out_file)
    return out_file


class _CropToROIInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='label/mask volume to crop')
    margin_mm = traits.Float(
        10.0, usedefault=True, desc='margin to add around the bounding box, in mm'
    )


class _CropToROIOutputSpec(TraitedSpec):
    out_file = File(desc='cropped volume')


class CropToROI(SimpleInterface):
    """Crop a volume to the bounding box of its nonzero voxels, plus a margin."""

    input_spec = _CropToROIInputSpec
    output_spec = _CropToROIOutputSpec

    def _run_interface(self, runtime):
        self._results['out_file'] = crop_to_roi(
            self.inputs.in_file,
            margin_mm=self.inputs.margin_mm,
            newpath=runtime.cwd,
        )
        return runtime


def crop_to_roi(in_file: str, margin_mm: float = 10.0, newpath: str | None = None) -> str:
    """
    Crop a volume to the bounding box of its nonzero voxels, plus a margin.

    The output affine is adjusted so world coordinates are preserved.

    Parameters
    ----------
    in_file : A 3D label/mask volume.
    margin_mm : Margin to add on each side of the bounding box, in millimeters.
    """
    import nibabel as nb
    import numpy as np

    img = nb.load(in_file)
    data = np.asanyarray(img.dataobj)
    nonzero = np.argwhere(data > 0)
    if nonzero.size == 0:
        raise ValueError(f'No nonzero voxels found in {in_file}')

    zooms = np.asarray(img.header.get_zooms()[:3], dtype=float)
    margin_vox = np.ceil(margin_mm / zooms).astype(int)
    lo = np.maximum(nonzero.min(axis=0) - margin_vox, 0)
    hi = np.minimum(nonzero.max(axis=0) + 1 + margin_vox, np.asarray(data.shape[:3]))

    roi = data[lo[0] : hi[0], lo[1] : hi[1], lo[2] : hi[2]]
    affine = img.affine.copy()
    # Preserve world coordinates by shifting the origin to the new corner of the cropped volume
    affine[:3, 3] = nb.affines.apply_affine(img.affine, lo)

    out_file = str((Path(newpath or '.') / 'roi_cropped.nii.gz').absolute())
    img.__class__(roi, affine, img.header).to_filename(out_file)
    return out_file


class _MergeLabelROIsInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiObject(
        File(exists=True), mandatory=True, desc='per-structure label volumes'
    )
    template = File(exists=True, mandatory=True, desc='volume defining the output grid')
    overlap_method = traits.Enum(
        'first',
        'last',
        'fail',
        usedefault=True,
        desc='how to resolve voxels claimed by more than one structure: keep the first '
        'claim, keep the last claim, or raise an error',
    )


class _MergeLabelROIsOutputSpec(TraitedSpec):
    out_file = File(desc='merged label segmentation on the template grid')
    overlap_file = File(desc='per-voxel count of structures claiming each voxel')


class MergeLabelROIs(SimpleInterface):
    """Merge per-structure label volumes onto a common template grid."""

    input_spec = _MergeLabelROIsInputSpec
    output_spec = _MergeLabelROIsOutputSpec

    def _run_interface(self, runtime):
        out_file, overlap_file = merge_label_rois(
            self.inputs.in_files,
            self.inputs.template,
            overlap_method=self.inputs.overlap_method,
            newpath=runtime.cwd,
        )
        self._results['out_file'] = out_file
        self._results['overlap_file'] = overlap_file
        return runtime


def merge_label_rois(
    in_files: list,
    template: str,
    overlap_method: str = 'first',
    newpath: str | None = None,
):
    """
    Merge per-structure label volumes onto a common template grid.

    Each input is expected to share the template's voxel grid (same orientation and
    spacing), differing only by a cropped field of view. Voxels are pasted into the
    template grid.

    ``overlap_method`` decides on how to handle overlapping voxels.

    Parameters
    ----------
    overlap_method : ``'first'`` keeps the earliest claim, ``'last'`` keeps the latest,
        ``'fail'`` raises on the first collision.

    Returns
    -------
    out_file : The merged label segmentation on the template grid.
    overlap_file : Per-voxel count of structures claiming each voxel (collision map).
    """
    import nibabel as nb
    import numpy as np

    tpl = nb.load(template)
    out = np.zeros(tpl.shape[:3], dtype='int16')
    counts = np.zeros(tpl.shape[:3], dtype='int16')
    inv = np.linalg.inv(tpl.affine)

    for fname in in_files:
        img = nb.load(fname)
        data = np.asanyarray(img.dataobj)
        xfm = inv @ img.affine
        if not np.allclose(xfm[:3, :3], np.eye(3), atol=1e-3):
            raise ValueError(f'{fname} grid is not aligned with the template grid')
        offset = np.rint(xfm[:3, 3]).astype(int)
        for ijk in np.argwhere(data != 0):
            tgt = ijk + offset
            if np.any(tgt < 0) or np.any(tgt >= np.asarray(out.shape)):
                continue
            i, j, k = tgt
            counts[i, j, k] += 1
            if counts[i, j, k] > 1:  # overlap found
                match overlap_method:
                    case 'fail':
                        raise ValueError(f'Overlapping ROIs at voxel {(i, j, k)}')
                    case 'first':
                        continue  # keep the existing label
                    case _:
                        pass
            out[i, j, k] = data[tuple(ijk)]

    newpath = Path(newpath or '.')
    out_file = str((newpath / 'seg.nii.gz').absolute())
    overlap_file = str((newpath / 'overlap.nii.gz').absolute())
    nb.Nifti1Image(out, tpl.affine, tpl.header).to_filename(out_file)
    nb.Nifti1Image(counts, tpl.affine, tpl.header).to_filename(overlap_file)
    return out_file, overlap_file
