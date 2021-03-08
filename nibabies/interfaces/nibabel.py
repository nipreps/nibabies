"""Nibabel-based interfaces to eventually upstream to NiWorkflows."""
from nipype.interfaces.base import (
    traits,
    TraitedSpec,
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
)


class _BinaryDilationInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="binary file to dilate")
    radius = traits.Float(3, usedefault=True, desc="structure element (ball) radius")
    iterations = traits.Range(low=0, value=1, usedefault=True, desc="repeat dilation")


class _BinaryDilationOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="the input file, after binary dilation")


class BinaryDilation(SimpleInterface):
    """Brain extraction for EPI and GRE data."""

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
