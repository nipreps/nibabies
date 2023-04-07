import nipype.interfaces.utility as niu
import nipype.pipeline.engine as pe


def init_infant_epi_reference_wf(
    omp_nthreads: int,
    is_sbref: bool = False,
    start_frame: int = 17,
    name: str = 'infant_epi_reference_wf',
) -> pe.Workflow:
    """
    Workflow to generate a reference map from one or more infant EPI images.

    If any single-band references are provided, the reference map will be calculated from those.

    If no single-band references are provided, the BOLD files are used.
    To account for potential increased motion on the start of image acquisition, this
    workflow discards a bigger chunk of the initial frames.

    Parameters
    ----------
    omp_nthreads
        Maximum number of threads an individual process may use
    is_sbref
        A single-band reference is provided.
    start_frame
        BOLD frame to start creating the reference map from. Any earlier frames are discarded.

    Inputs
    ------
    bold_file
        BOLD EPI file
    sbref_file
        single-band reference EPI

    Outputs
    -------
    boldref_file
        The generated reference map
    boldref_mask
        Binary brain mask of the ``boldref_file``
    boldref_xfm
        Rigid-body transforms in LTA format

    """
    from niworkflows.workflows.epi.refmap import init_epi_reference_wf
    from sdcflows.interfaces.brainmask import BrainExtraction

    wf = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['epi_file']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['boldref_file', 'boldref_mask']),
        name='outputnode',
    )

    epi_reference_wf = init_epi_reference_wf(omp_nthreads, auto_bold_nss=False)

    boldref_mask = pe.Node(BrainExtraction(), name='boldref_mask')

    # fmt:off
    wf.connect([
        (inputnode, epi_reference_wf, [('epi_file', 'inputnode.in_files')]),
        (epi_reference_wf, boldref_mask, [('outputnode.epi_ref_file', 'in_file')]),
        (epi_reference_wf, outputnode, [('outputnode.epi_ref_file', 'boldref_file')]),
        (boldref_mask, outputnode, [('out_mask', 'boldref_mask')]),
    ])
    # fmt:on
    if not is_sbref:
        select_frames = pe.Node(
            niu.Function(function=_select_frames, output_names=['t_masks']),
            name='select_frames',
        )
        select_frames.inputs.start_frame = start_frame
        # fmt:off
        wf.connect([
            (inputnode, select_frames, [('epi_file', 'in_file')]),
            (select_frames, epi_reference_wf, [('t_masks', 'inputnode.t_masks')]),
        ])
        # fmt:on
    else:
        # Won't be used but needed to placate iternode
        # To consider: Add a check to ensure this is a 3D file
        epi_reference_wf.inputs.inputnode.t_masks = [True]
    return wf


def _select_frames(in_file: str, start_frame: int) -> list:
    import nibabel as nb
    import numpy as np

    img = nb.load(in_file)
    img_len = img.shape[3]
    if start_frame >= img_len:
        start_frame = img_len - 1
    t_mask = np.array([False] * img_len, dtype=bool)
    t_mask[start_frame:] = True
    return list(t_mask)
