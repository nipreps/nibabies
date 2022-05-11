from pathlib import Path

import nibabel as nb
import nibabel.gifti as ngi
import numpy as np
from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)

from .. import __version__


class _MaskGiftiInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="Input GIFTI (n-darrays)")
    mask_file = File(exists=True, mandatory=True, desc="Input mask (single binary darray)")
    threshold = traits.Float(
        desc="If mask is probabilistic, inclusion limit",
    )
    metadata = traits.Dict(
        desc="Metadata to insert into GIFTI",
    )


class _MaskGiftiOutputSpec(TraitedSpec):
    out_file = File(desc="Masked file")


class MaskGifti(SimpleInterface):
    """Mask file across GIFTI darrays"""

    input_spec = _MaskGiftiInputSpec
    output_spec = _MaskGiftiOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = _mask_gifti(
            self.inputs.in_file,
            self.inputs.mask_file,
            threshold=self.inputs.threshold or None,
            metadata=self.inputs.metadata,
            newpath=runtime.cwd,
        )
        return runtime


def _mask_gifti(in_file, mask_file, *, threshold=None, metadata=None, newpath=None):
    """
    Mask and create a GIFTI image.
    """
    metadata = metadata or {}

    img = nb.load(in_file)
    mask = nb.load(mask_file).agg_data()

    indices = np.nonzero(mask)[0]
    if threshold is not None:
        indices = np.where(mask > threshold)[0]

    data = img.agg_data()
    if isinstance(data, tuple):
        try:
            data = np.vstack(data)
        except Exception:
            raise NotImplementedError(f"Tricky GIFTI: {in_file} not supported.")
    else:
        data = data.T
    masked = data[:, indices]

    # rather than creating new GiftiDataArrays, just modify the data directly
    # and preserve the existing attributes
    for i, darr in enumerate(img.darrays):
        darr.data = masked[i]
        darr.dims = list(masked[i].shape)

    # Finalize by adding additional metadata to file
    metad = {
        **{"CreatedBy": f"MaskGifti (NiBabies-{__version__})"},
        **metadata,
    }
    if int(nb.__version__[0]) >= 4:  # API will change in 4.0.0
        existing_meta = img.meta or {}
        img.meta = ngi.GiftiMetaData({**metad, **existing_meta})
    else:
        meta = img.meta.data or []
        for k, v in metad.items():
            meta.append(ngi.GiftiNVPairs(k, v))
        img.meta.data = meta

    if newpath is None:
        newpath = Path()
    out_file = str((Path(newpath) / f"masked_{Path(in_file).name}").absolute())
    nb.save(img, out_file)
    return out_file
