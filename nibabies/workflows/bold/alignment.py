"""
Subcortical alignment into MNI space
"""

from nibabies.interfaces.nibabel import MergeROIs


def init_subcortical_mni_alignment_wf(*, vol_sigma=0.8, name='subcortical_mni_alignment_wf'):
    """
    Align individual subcortical structures into MNI space.

    This is a nipype workflow port of the DCAN infant pipeline:
    https://github.com/DCAN-Labs/dcan-infant-pipeline/blob\
    /master/fMRISurface/scripts/SubcorticalAlign_ROIs.sh


    Parameters
    ----------
    name : :obj:`str`
        Name of the workflow
    vol_sigma : :obj:`float`
        The sigma for the gaussian volume smoothing kernel, in mm

    Inputs
    ------
    bold_file : :obj:`str`
        BOLD file
    bold_roi : :obj:`str`
        File containing ROIs in BOLD space
    atlas_roi : :obj:`str`
        File containing ROIs in atlas space
    std_xfm : :obj:`str`
        File containing transform to the standard (MNI) space

    Outputs
    -------
    subcortical_file : :obj:`str`
        Volume file containing all ROIs individually aligned to standard
    """
    from nipype.pipeline import engine as pe
    from nipype.interfaces import utility as niu, fsl
    from ...interfaces.workbench import (
        CiftiCreateDenseTimeseries,
        CiftiCreateLabel,
        CiftiDilate,
        CiftiResample,
        CiftiSeparate,
        CiftiSmooth,
        VolumeAffineResample,
        VolumeAllLabelsToROIs,
        VolumeLabelExportTable,
        VolumeLabelImport,
    )
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["bold_file", "bold_roi", "atlas_roi", "atlas_xfm"]),
        name="inputnode",
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=["subcortical_file"]), name='outputnode')

    applyxfm_atlas = pe.Node(fsl.ApplyXFM(), name="applyxfm_atlas")
    vol_resample = pe.Node(
        VolumeAffineResample(method="ENCLOSING_VOXEL", flirt=True),
        name="vol_resample"
    )
    subj_rois = pe.Node(VolumeAllLabelsToROIs(label_map=1), name="subj_rois")
    split_rois = pe.Node(fsl.Split(dimension="t"), name="split_rois")
    atlas_rois = pe.Node(VolumeAllLabelsToROIs(label_map=1), name="atlas_rois")
    split_atlas_rois = pe.Node(fsl.Split(dimension="t"), name="split_atlas_rois")
    atlas_labels = pe.Node(VolumeLabelExportTable(label_map=1), name="atlas_labels")
    parse_labels = pe.Node(
        niu.Function(function=parse_roi_labels, output_names=["structures", "label_ids"]),
        name="parse_labels",
    )

    # The following is wrapped in a for-loop, iterating across each roi
    # Instead, we will use MapNodes and iter across the varying inputs
    roi2atlas = pe.MapNode(
        fsl.FLIRT(
            searchr_x=[-20, 20],
            searchr_y=[-20, 20],
            searchr_z=[-20, 20],
            interp="nearestneighbour",
        ),
        name="roi2atlas",
        iterfield=["in_file", "reference"],
    )
    applyxfm_roi = pe.MapNode(
        fsl.ApplyXFM(interp="spline"),
        iterfield=["reference", "in_matrix_file"],
        name='applyxfm_roi',
        mem_gb=4,
    )
    bold_mask_roi = pe.MapNode(
        fsl.ApplyMask(),
        iterfield=["in_file", "mask_file"],
        name='bold_mask_roi',
    )
    mul_roi = pe.MapNode(
        fsl.BinaryMaths(operation="mul"),
        iterfield=["in_file", "operand_value"],
        name='mul_roi',
    )
    mul_atlas_roi = pe.MapNode(
        fsl.BinaryMaths(operation="mul"),
        iterfield=["in_file", "operand_value"],
        name='mul_atlas_roi',
    )
    vol_label = pe.MapNode(
        VolumeLabelImport(drop_unused_labels=True),
        iterfield=["in_file"],
        name='vol_label',
    )
    vol_atlas_label = pe.MapNode(
        VolumeLabelImport(drop_unused_labels=True),
        iterfield=["in_file"],
        name='vol_atlas_label',
    )
    create_dtseries = pe.MapNode(
        CiftiCreateDenseTimeseries(),
        iterfield=["volume_data", "volume_structure_labels"],
        name='create_dtseries',
    )
    create_label = pe.MapNode(
        CiftiCreateLabel(),
        iterfield=["volume_label", "structure_label_volume"],
        name='create_label',
    )
    dilate = pe.MapNode(
        CiftiDilate(direction="COLUMN", surface_distance=0, volume_distance=10),
        iterfield=["in_file"],
        name="dilate"
    )
    resample = pe.MapNode(
        CiftiResample(
            direction="COLUMN",
            template_direction="COLUMN",
            surface_method="ADAP_BARY_AREA",
            volume_method="CUBIC",
            volume_predilate=10,
        ),
        iterfield=["in_file", "template"],
        name='resample',
    )
    smooth = pe.MapNode(
        CiftiSmooth(direction="COLUMN", fix_zeros_vol=True, sigma_surf=0, sigma_vol=vol_sigma),
        iterfield=["in_file"],
        name="smooth"
    )
    separate = pe.MapNode(
        CiftiSeparate(direction="COLUMN", volume_all_file='volume_all.nii.gz'),
        iterfield=["in_file"],
        name="separate"
    )
    fmt_agg_rois = pe.Node(
        niu.Function(
            function=format_agg_rois,
            output_names=["first_image", "op_files", "op_string"],
        ),
        name='fmt_agg_rois',
    )
    agg_rois = pe.Node(fsl.MultiImageMaths(), name='agg_rois')
    merge_rois = pe.Node(MergeROIs(), name="merge_rois")

    workflow = Workflow(name=name)
    # fmt: off
    workflow.connect([
        (inputnode, applyxfm_atlas, [
            ("bold_file", "in_file"),
            ("atlas_roi", "reference"),
            ("atlas_xfm", "in_matrix_file")]),
        (inputnode, vol_resample, [
            ("bold_roi", "in_file"),
            ("atlas_xfm", "affine"),
            ("bold_roi", "flirt_source_volume")]),
        (applyxfm_atlas, vol_resample, [
            ("out_file", "volume_space"),
            ("out_file", "flirt_target_volume")]),
        (inputnode, subj_rois, [("bold_roi", "in_file")]),
        (inputnode, atlas_rois, [("atlas_roi", "in_file")]),
        (subj_rois, split_rois, [("out_file", "in_file")]),
        (atlas_rois, split_atlas_rois, [("out_file", "in_file")]),
        (inputnode, atlas_labels, [("atlas_roi", "in_file")]),
        (atlas_labels, parse_labels, [("out_file", "label_file")]),
        # for loop across ROIs
        (split_rois, roi2atlas, [("out_files", "in_file")]),
        (split_atlas_rois, roi2atlas, [("out_files", "reference")]),
        (inputnode, applyxfm_roi, [("bold_file", "in_file")]),
        (split_atlas_rois, applyxfm_roi, [("out_files", "reference")]),
        (roi2atlas, applyxfm_roi, [("out_matrix_file", "in_matrix_file")]),
        (applyxfm_roi, bold_mask_roi, [("out_file", "in_file")]),
        (roi2atlas, bold_mask_roi, [("out_file", "mask_file")]),
        (roi2atlas, mul_roi, [("out_file", "in_file")]),
        (parse_labels, mul_roi, [("label_ids", "operand_value")]),
        (split_atlas_rois, mul_atlas_roi, [("out_files", "in_file")]),
        (parse_labels, mul_atlas_roi, [("label_ids", "operand_value")]),
        (mul_roi, vol_label, [("out_file", "in_file")]),
        (atlas_labels, vol_label, [("out_file", "label_list_file")]),
        (mul_atlas_roi, vol_atlas_label, [("out_file", "in_file")]),
        (atlas_labels, vol_atlas_label, [("out_file", "label_list_file")]),
        (bold_mask_roi, create_dtseries, [("out_file", "volume_data")]),
        (vol_label, create_dtseries, [("out_file", "volume_structure_labels")]),
        (vol_atlas_label, create_label, [
            ("out_file", "volume_label"),
            ("out_file", "structure_label_volume")]),
        (create_dtseries, dilate, [("out_file", "in_file")]),
        (dilate, resample, [("out_file", "in_file")]),
        (create_label, resample, [("out_file", "template")]),
        (resample, smooth, [("out_file", "in_file")]),
        (smooth, separate, [("out_file", "in_file")]),
        # end loop
        (mul_roi, fmt_agg_rois, [("out_file", "rois")]),
        (fmt_agg_rois, agg_rois, [
            ("first_image", "in_file"),
            ("op_files", "operand_files"),
            ("op_string", "op_string")]),
        (separate, merge_rois, [("volume_all_file", "in_files")]),
        (merge_rois, outputnode, [("out_file", "subcortical_file")]),
    ])
    # fmt: on
    return workflow


def parse_roi_labels(label_file):
    """
    Parse a label file composed of one or more sets of:
    <labelname>
    <key> <red> <green> <blue> <alpha>

    Return a list of structure names and label keys.

    >>> structs, ids = parse_roi_labels(test_data / "labelfile.txt")
    >>> structs
    ['CEREBELLUM_LEFT', 'THALAMUS_LEFT', 'CAUDATE_LEFT']
    >>> ids
    [8, 10, 11]
    """

    with open(label_file) as fp:
        lines = fp.readlines()
    if len(lines) % 2 == 1:
        raise RuntimeError("Label file is incomplete or invalid")
    structs, label_ids = [], []
    for idx, line in enumerate(lines):
        if idx % 2 == 0:
            structs.append(line.strip())
        else:
            label_ids.append(int(line.split(' ', 1)[0]))
    return structs, label_ids


def format_agg_rois(rois):
    """
    Helper function to format MultiImageMaths command.

    Parameters
    ----------
    rois : `list` of `str`s
        List of files

    Returns
    -------
    first_image
    op_files
    op_string

    """
    return rois[0], rois[1:], ("-add %s " * (len(rois) - 1)).strip()
