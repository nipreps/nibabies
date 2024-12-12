import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe
from niworkflows.engine.workflows import LiterateWorkflow


def init_anat_preproc_wf(
    *,
    bspline_fitting_distance: int = 200,
    name: str = 'anat_preproc_wf',
) -> LiterateWorkflow:
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

    wf = LiterateWorkflow(name=name)
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


def init_anat_csf_norm_wf(name='anat_csf_norm_wf') -> LiterateWorkflow:
    """Replace low intensity voxels within the CSF mask with the median value."""

    workflow = LiterateWorkflow(name=name)
    inputnode = niu.IdentityInterface(fields=['anat_preproc', 'anat_dseg'], name='inputnode')
    outputnode = niu.IdentityInterface(fields=['anat_preproc'], name='outputnode')

    applymask = pe.Node(ApplyMask(), name='applymask')

    norm2median = pe.Node(niu.Function(function=_normalize_csf), name='norm2median')
    # 1. mask brain with CSF mask
    # fslmaths input.nii.gz -mas aseg_label-CSF_mask.nii.gz input_CSF.nii.gz
    # 2. get median intensity of nonzero voxels in mask
    # fslstats input_CSF.nii.gz -P 50
    # 3. normalize CSF-masked T2w to the median
    # fslmaths input_CSF.nii.gz -bin -mul <median intensity from (> median_CSF.nii.gz
    # 4. make the modified T2w, setting voxel intensity to be the max between the original T2w's,
    # and the normalized mask from (3)'s:
    # fslmaths input.nii.gz -max median_CSF.nii.gz input_floorCSF.nii.gz

    return workflow


def _normalize_csf(in_file): ...
