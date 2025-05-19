import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
from niworkflows.engine import Workflow, tag


@tag('anat.preproc')
def init_anat_preproc_wf(
    *,
    bspline_fitting_distance: int = 200,
    name: str = 'anat_preproc_wf',
) -> Workflow:
    """Polish up raw anatomical data.

    This workflow accepts T1w/T2w images as inputs (either raw or a merged template) and performs:
    - Intensity clipping
    - N4 Bias Field Correction

    The outputs of this workflow will be a structural reference used for subsequent processing.

    Inputs
    ------
    in_anat : :obj:`str`
        A single volume T1w/T2w image

    Outputs
    -------
    anat_preproc: :obj:`str`
        Preprocessed anatomical image (Denoising/INU/Clipping)
    """
    from nipype.interfaces.ants import N4BiasFieldCorrection
    from niworkflows.interfaces.header import ValidateImage
    from niworkflows.interfaces.nibabel import IntensityClip

    wf = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['in_anat']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['anat_preproc']),
        name='outputnode',
    )

    # validate image
    validate = pe.Node(ValidateImage(), name='anat_validate', run_without_submitting=True)
    clip = pe.Node(IntensityClip(p_min=10.0, p_max=99.5), name='clip')
    n4_correct = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            bspline_fitting_distance=bspline_fitting_distance,
            save_bias=True,
            copy_header=True,
            n_iterations=[50] * 5,
            convergence_threshold=1e-7,
            rescale_intensities=True,
            shrink_factor=4,
        ),
        name='n4_correct',
    )
    final_clip = pe.Node(IntensityClip(p_min=5.0, p_max=99.5), name='final_clip')

    wf.connect([
        (inputnode, validate, [('in_anat', 'in_file')]),
        (validate, clip, [('out_file', 'in_file')]),
        (clip, n4_correct, [('out_file', 'input_image')]),
        (n4_correct, final_clip, [('output_image', 'in_file')]),
        (final_clip, outputnode, [('out_file', 'anat_preproc')]),
    ])  # fmt:skip
    return wf


@tag('anat.csf_norm')
def init_csf_norm_wf(name: str = 'csf_norm_wf') -> Workflow:
    """Replace low intensity voxels within the CSF mask with the median value."""

    workflow = Workflow(name=name)
    workflow.__desc__ = (
        'The CSF mask was used to normalize the anatomical template by the median of voxels '
        'within the mask.'
    )
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['anat_preproc', 'anat_tpms']),
        name='inputnode',
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=['anat_preproc']), name='outputnode')

    # select CSF from BIDS-ordered list (GM, WM, CSF)
    select_csf = pe.Node(niu.Select(index=2), name='select_csf')
    norm_csf = pe.Node(niu.Function(function=_normalize_roi), name='norm_csf')

    workflow.connect([
        (inputnode, select_csf, [('anat_tpms', 'inlist')]),
        (select_csf, norm_csf, [('out', 'mask_file')]),
        (inputnode, norm_csf, [('anat_preproc', 'in_file')]),
        (norm_csf, outputnode, [('out', 'anat_preproc')]),
    ])  # fmt:skip

    return workflow


def init_conform_derivative_wf(
    *, in_file: str = None, name: str = 'conform_derivative_wf'
) -> pe.Workflow:
    """
    Ensure derivatives share the same space as anatomical references.

    This workflow is used when a derivative is provided without a reference.
    """
    from niworkflows.interfaces.header import MatchHeader
    from niworkflows.interfaces.images import Conform, TemplateDimensions

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['in_file', 'ref_file']), name='inputnode')
    inputnode.inputs.in_file = in_file
    outputnode = pe.Node(niu.IdentityInterface(fields=['out_file']), name='outputnode')

    ref_dims = pe.Node(TemplateDimensions(), name='ref_dims')
    conform = pe.Node(Conform(), name='conform')
    # Avoid mismatch tolerance from input
    match_header = pe.Node(MatchHeader(), name='match_header')

    workflow.connect([
        (inputnode, ref_dims, [('ref_file', 'anat_list')]),
        (ref_dims, conform, [
            ('target_zooms', 'target_zooms'),
            ('target_shape', 'target_shape'),
        ]),
        (inputnode, conform, [('in_file', 'in_file')]),
        (conform, match_header, [('out_file', 'in_file')]),
        (inputnode, match_header, [('ref_file', 'reference')]),
        (match_header, outputnode, [('out_file', 'out_file')]),
    ])  # fmt:skip
    return workflow


def _normalize_roi(in_file, mask_file, threshold=0.2, out_file=None):
    """Normalize low intensity voxels that fall within a given mask."""
    import nibabel as nb
    import numpy as np

    img = nb.load(in_file)
    img_data = np.asanyarray(img.dataobj)
    mask_img = nb.load(mask_file)
    # binary mask
    bin_mask = np.asanyarray(mask_img.dataobj) > threshold
    mask_data = bin_mask * img_data
    masked_data = mask_data[mask_data > 0]

    median = np.median(masked_data).astype(masked_data.dtype)
    normed_data = np.maximum(img_data, bin_mask * median)

    oimg = img.__class__(normed_data, img.affine, img.header)
    if not out_file:
        from nipype.utils.filemanip import fname_presuffix

        out_file = fname_presuffix(in_file, suffix='normed')
    oimg.to_filename(out_file)
    return out_file
