from nipype.interfaces.base import (
    traits,
    TraitedSpec,
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    OutputMultiObject,
    InputMultiObject,
)


class _IntensityClipInputSpec(BaseInterfaceInputSpec):
    in_file = File(
        exists=True, mandatory=True, desc="file which intensity will be clipped"
    )
    p_min = traits.Float(35.0, usedefault=True, desc="percentile for the lower bound")
    p_max = traits.Float(99.98, usedefault=True, desc="percentile for the upper bound")
    nonnegative = traits.Bool(
        True, usedefault=True, desc="whether input intensities must be positive"
    )
    dtype = traits.Enum(
        "int16", "float32", "uint8", usedefault=True, desc="output datatype"
    )
    invert = traits.Bool(False, usedefault=True, desc="finalize by inverting contrast")


class _IntensityClipOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="file after clipping")


class IntensityClip(SimpleInterface):
    """Clip the intensity range as prescribed by the percentiles."""

    input_spec = _IntensityClipInputSpec
    output_spec = _IntensityClipOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = _advanced_clip(
            self.inputs.in_file,
            p_min=self.inputs.p_min,
            p_max=self.inputs.p_max,
            nonnegative=self.inputs.nonnegative,
            dtype=self.inputs.dtype,
            invert=self.inputs.invert,
            newpath=runtime.cwd,
        )
        return runtime


def _advanced_clip(
    in_file, p_min=35, p_max=99.98, nonnegative=True, dtype="int16", invert=False, newpath=None,
):
    """
    Remove outliers at both ends of the intensity distribution and fit into a given dtype.
    This interface tries to emulate ANTs workflows' massaging that truncate images into
    the 0-255 range, and applies percentiles for clipping images.
    For image registration, normalizing the intensity into a compact range (e.g., uint8)
    is generally advised.
    To more robustly determine the clipping thresholds, data are removed of spikes
    with a median filter.
    Once the thresholds are calculated, the denoised data are thrown away and the thresholds
    are applied on the original image.
    """
    from pathlib import Path
    import nibabel as nb
    import numpy as np
    from scipy import ndimage
    from skimage.morphology import ball

    out_file = (Path(newpath or "") / "clipped.nii.gz").absolute()

    # Load data
    img = nb.squeeze_image(nb.load(in_file))
    assert len(img.shape) == 3, "Not a 3D image"
    data = img.get_fdata(dtype="float32")

    # Calculate stats on denoised version, to preempt outliers from biasing
    denoised = ndimage.median_filter(data, footprint=ball(3))

    a_min = np.percentile(
        denoised[denoised > 0] if nonnegative else denoised,
        p_min
    )
    a_max = np.percentile(
        denoised[denoised > 0] if nonnegative else denoised,
        p_max
    )

    # Clip and cast
    data = np.clip(data, a_min=a_min, a_max=a_max)
    data -= data.min()
    data /= data.max()

    if invert:
        data = 1.0 - data

    if dtype in ("uint8", "int16"):
        data = np.round(255 * data).astype(dtype)

    hdr = img.header.copy()
    hdr.set_data_dtype(dtype)
    img.__class__(data, img.affine, hdr).to_filename(out_file)

    return str(out_file)
