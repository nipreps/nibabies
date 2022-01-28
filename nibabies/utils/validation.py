import sys
import nibabel as nb
import numpy as np


def validate_t1w_derivatives(t1w_template, *, anat_mask=None, anat_aseg=None, atol=1e-5):
    """
    Validate anatomical derivatives.
    This function compares the input T1w's orientation and physical space to each derivative.

    Parameters
    ----------
    t1w_template : str
        T1w template
    anat_mask : str or None
        Precomputed anatomical brain mask
    anat_aseg : str or None
        Precomputed anatomical segmentations
    atol : float
        Absolute error tolerance between image origins

    Returns
    -------
    validated : dict
        A dictionary composed of derivative keys and validated filename values.
        Derivatives that failed validation will not be included.
    """

    validated = {}
    # T1w information
    t1w = nb.load(t1w_template)
    expected_ort = nb.aff2axcodes(t1w.affine)

    # Ensure orientation
    for name, deriv_fl in zip(('anat_mask', 'anat_aseg'), (anat_mask, anat_aseg)):
        if deriv_fl is None:
            continue
        img = nb.load(deriv_fl)
        if nb.aff2axcodes(img.affine) != expected_ort:
            print(
                f"Orientation mismatch between {name} <{deriv_fl}> and T1w <{t1w_template}>",
                file=sys.stderr,
            )
            continue
        if img.shape != t1w.shape or not np.allclose(t1w.affine, img.affine, atol=atol):
            print(
                f"Physical space mismatch between {name} <{deriv_fl}> and T1w <{t1w_template}>",
                file=sys.stderr,
            )
            continue
        validated[name] = deriv_fl
    return validated
