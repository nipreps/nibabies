from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)
from nipype.utils.filemanip import fname_presuffix


class _DetectReferenceFrameInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='BOLD timeseries')
    ref_frame_start = traits.Int(
        mandatory=True, desc='Frame to start looking for a low-motion reference frame'
    )
    dummy_scans = traits.Either(
        None, traits.Int, usedefault=True, desc='Number of non-steady-state scans at the start.'
    )


class _DetectReferenceFrameOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='Single reference frame')
    frame_idx = traits.Int(desc='Frame index used')


class DetectReferenceFrame(SimpleInterface):
    """Select one reference frame for HMC and SDC correction"""

    input_spec = _DetectReferenceFrameInputSpec
    output_spec = _DetectReferenceFrameOutputSpec

    def _run_interface(self, runtime):
        out_path = fname_presuffix(self.inputs.in_file, suffix='_refframe', newpath=runtime.cwd)
        out_file, frame_idx = _detect_reference_frame(
            in_file=self.inputs.in_file,
            ref_frame_start=self.inputs.ref_frame_start,
            out_file=out_path,
            dummy_scans=self.inputs.dummy_scans,
        )
        self._results['out_file'] = out_file
        self._results['frame_idx'] = frame_idx
        return runtime


def _detect_reference_frame(
    in_file: str, ref_frame_start: int, out_file: str, dummy_scans: int | None = None
) -> tuple[str, int]:
    import warnings

    import nibabel as nb
    import numpy as np

    start_frame = max(ref_frame_start, dummy_scans) if dummy_scans else ref_frame_start

    img = nb.load(in_file)
    ts = img.get_fdata(dtype='float32')
    img_len = ts.shape[3]
    if start_frame >= img_len:
        warnings.warn(
            f'Caculating the BOLD reference starting on frame {start_frame} but only {img_len} '
            'volumes in BOLD file, so using last volume.',
            stacklevel=1,
        )
        start_frame = img_len - 1

    ts = ts[..., start_frame:]
    ts /= np.max(ts)
    ts_mean = np.nanmean(ts, axis=3)
    chosen_frame = np.argmin(np.sum((ts - ts_mean[..., np.newaxis]) ** 2, axis=(0, 1, 2)))
    chosen_frame_img = nb.Nifti1Image(np.squeeze(ts[..., chosen_frame]), affine=img.affine)
    nb.save(chosen_frame_img, out_file)
    return out_file, chosen_frame
