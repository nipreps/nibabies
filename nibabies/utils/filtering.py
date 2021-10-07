"""Signal processing filters."""


def truncation(
    in_file,
    clip_max=99.9,
    dtype="int16",
    out_file=None,
    out_max=1000,
    out_min=0,
    percentiles=(0.1, 95),
):
    """Truncate and clip the input image intensities."""
    from pathlib import Path
    import numpy as np
    import nibabel as nb
    from nipype.utils.filemanip import fname_presuffix

    try:
        info = np.iinfo(dtype)
    except ValueError:
        info = np.finfo(dtype)

    img = nb.load(in_file)
    hdr = img.header.copy()
    hdr.set_data_dtype(dtype)

    data = img.get_fdata()

    out_min = max(out_min, info.min)
    out_max = min(out_max, info.max)

    a_min = np.percentile(data.reshape(-1), percentiles[0])
    data -= a_min
    a_max = np.percentile(data.reshape(-1), percentiles[1])
    data *= out_max / a_max
    data = np.clip(data, info.min, info.max)

    if clip_max is not None:
        data = np.clip(data, 0, np.percentile(data.reshape(-1), clip_max))

    if out_file is None:
        out_file = fname_presuffix(Path(in_file).name, suffix="_trunc")

    out_file = str(Path(out_file).absolute())
    img.__class__(data.astype(dtype), img.affine, hdr).to_filename(out_file)
    return out_file


def gaussian_filter(in_file, sigma=None, out_file=None):
    """Filter input image by convolving with a Gaussian kernel."""
    from pathlib import Path
    import numpy as np
    import nibabel as nb
    from scipy.ndimage import gaussian_filter
    from nipype.utils.filemanip import fname_presuffix

    if out_file is None:
        out_file = fname_presuffix(Path(in_file).name, suffix="_gauss")
    out_file = str(Path(out_file).absolute())

    img = nb.load(in_file)
    if sigma is None:
        sigma = tuple(np.array(img.header.get_zooms()[:3]) * 2.0)
    img.__class__(gaussian_filter(img.dataobj, sigma), img.affine, img.header).to_filename(
        out_file
    )
    return out_file
