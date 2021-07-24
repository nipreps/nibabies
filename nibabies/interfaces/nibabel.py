"""Nibabel-based interfaces to eventually upstream to NiWorkflows."""
from nipype.interfaces.base import (
    traits,
    TraitedSpec,
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    InputMultiObject,
)


class _BinaryDilationInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="binary file to dilate")
    radius = traits.Float(3, usedefault=True, desc="structure element (ball) radius")
    iterations = traits.Range(low=0, value=1, usedefault=True, desc="repeat dilation")


class _BinaryDilationOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="the input file, after binary dilation")


class BinaryDilation(SimpleInterface):
    """Morphological binary dilation using Scipy."""

    input_spec = _BinaryDilationInputSpec
    output_spec = _BinaryDilationOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = _dilate(
            self.inputs.in_file,
            radius=self.inputs.radius,
            iterations=self.inputs.iterations,
            newpath=runtime.cwd,
        )
        return runtime


class MergeROIsInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiObject(File(exists=True), desc="ROI files to be merged")


class MergeROIsOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="NIfTI containing all ROIs")


class MergeROIs(SimpleInterface):
    """Combine multiple region of interest files (3D or 4D) into a single file"""
    input_spec = MergeROIsInputSpec
    output_spec = MergeROIsOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = _merge_rois(self.inputs.in_files, newpath=runtime.cwd)
        return runtime


def _dilate(in_file, radius=3, iterations=1, newpath=None):
    """Dilate (binary) input mask."""
    from pathlib import Path
    import numpy as np
    import nibabel as nb
    from scipy import ndimage
    from skimage.morphology import ball
    from nipype.utils.filemanip import fname_presuffix

    mask = nb.load(in_file)
    newdata = ndimage.binary_dilation(
        np.asanyarray(mask.dataobj) > 0,
        structure=ball(radius),
        iterations=iterations,
    )

    hdr = mask.header.copy()
    hdr.set_data_dtype("uint8")
    out_file = fname_presuffix(in_file, suffix="_dil", newpath=newpath or Path.cwd())
    mask.__class__(newdata.astype("uint8"), mask.affine, hdr).to_filename(out_file)
    return out_file


def _merge_rois(in_files, newpath=None):
    """
    Aggregate individual 4D ROI files together into a single subcortical NIfTI.
    All ROI images are sanity checked with regards to:
    1) Shape
    2) Affine
    3) Overlap

    If any of these checks fail, an ``AssertionError`` will be raised.
    """
    from pathlib import Path
    import nibabel as nb
    import numpy as np

    img = nb.load(in_files[0])
    data = np.array(img.dataobj)
    affine = img.affine
    header = img.header

    nonzero = np.any(data, axis=3)
    for roi in in_files[1:]:
        img = nb.load(roi)
        assert img.shape == data.shape, "Mismatch in image shape"
        assert np.allclose(img.affine, affine), "Mismatch in affine"
        roi_data = np.asanyarray(img.dataobj)
        roi_nonzero = np.any(roi_data, axis=3)
        assert not np.any(roi_nonzero & nonzero), "Overlapping ROIs"
        nonzero |= roi_nonzero
        data += roi_data
        del roi_data

    if newpath is None:
        newpath = Path()
    out_file = str((Path(newpath) / "combined.nii.gz").absolute())
    img.__class__(data, affine, header).to_filename(out_file)
    return out_file
