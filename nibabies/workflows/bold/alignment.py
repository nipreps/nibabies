"""
Subcortical alignment into MNI space
"""

from nipype.interfaces import fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from nibabies.data import load as load_data
from nibabies.interfaces.nibabel import CropToROI, MergeLabelROIs
from nibabies.interfaces.workbench import VolumeLabelImport


def init_subcortical_rois_wf(*, name: str = 'subcortical_rois_wf'):
    """
    Refine segmentations into volumes of expected CIFTI subcortical structures.


    Parameters
    ----------
    name : :obj:`str`
        Name of the workflow

    Inputs
    ------
    MNIInfant_aseg : :obj:`str`
        FreeSurfer's aseg in MNIInfant space

    Outputs
    -------
    MNIInfant_rois : :obj:`str`
        Subcortical ROIs in MNIInfant space
    MNI152_rois : :obj:`str`
        Subcortical ROIs in `MNI152NLin6Asym` space
    """
    from niworkflows.interfaces.nibabel import MapLabels
    from templateflow.api import get as get_template

    # TODO: Implement BOLD refinement once InfantFS outputs subj/mri/wmparc.mgz
    # The code is found at
    # https://github.com/DCAN-Labs/dcan-infant-pipeline/blob/
    # 0e9c2fe32fb4a5032d0a2a3e0905ad97fa52b398/PostFreeSurfer/scripts/
    # FreeSurfer2CaretConvertAndRegisterNonlinear.sh
    # Lines 70-78 & 116-127
    # #
    # For now, just use the aseg

    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['MNIInfant_aseg']), name='inputnode')
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['MNIInfant_rois', 'MNI152_rois']),
        name='outputnode',
    )
    # Fetch the HCP volumetric template
    tpl_rois = get_template(
        'MNI152NLin6Asym', resolution=2, atlas='HCP', suffix='dseg', raise_empty=True
    )
    outputnode.inputs.MNI152_rois = tpl_rois

    # This will only used for the wmparc in subject space
    # For now, define it and don't run it
    # TODO: Move to TemplateFlow

    # tpl_avgwmparc = load_data("tpl-MNI152NLin6Asym_res-01_desc-avgwmparc_dseg.nii.gz")
    # applywarp_tpl = pe.Node(
    #     fsl.ApplyWarp(in_file=tpl_avgwmparc, ref_file=tpl_rois, interp="nn"),
    #     name="applywarp_std"
    # )

    map_labels = pe.Node(
        MapLabels(mappings_file=load_data('FreeSurferLabelRemappings.json')),
        name='map_labels',
    )

    refine_bold_rois = pe.Node(
        VolumeLabelImport(
            label_list_file=load_data('FreeSurferSubcorticalLabelTableLut.txt'),
            discard_others=True,
        ),
        name='refine_bold_rois',
    )

    # fmt: off
    workflow.connect([
        (inputnode, map_labels, [('MNIInfant_aseg', 'in_file')]),
        (map_labels, refine_bold_rois, [('out_file', 'in_file')]),
        # (applywarp_tpl, refine_std_rois, [("out_file", "in_file")]),
        (refine_bold_rois, outputnode, [('out_file', 'MNIInfant_rois')]),
    ])
    # fmt: on
    return workflow


def init_subcortical_mni_alignment_wf(
    *,
    vol_sigma: float = 0.8,
    debug: bool = False,
    name: str = 'subcortical_mni_alignment_wf',
):
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
    debug : :obj:`bool`
        Generate the additional collision-overlay QC reportlet

    Inputs
    ------
    MNIInfant_bold : :obj:`str`
        BOLD file in MNI Infant space
    MNIInfant_rois : :obj:`str`
        File containing ROIs in MNI Infant space
    MNI152_rois : :obj:`str`
        File containing ROIs in MNI152NLin6Asym space

    Outputs
    -------
    subcortical_volume : :obj:`str`
        Volume file containing all ROIs individually aligned to standard
    subcortical_labels : :obj:`str`
        Volume file containing all labels
    subcortical_seg : :obj:`str`
        Label volume of the aligned subcortical structures on the standard grid,
        for QC against the reference atlas ROIs
    out_report : :obj:`str`
        QC reportlet (SVG) overlaying the aligned segmentation against the reference ROIs
    out_report_overlap : :obj:`str`
        QC reportlet (SVG) overlaying the segmentation with the inter-structure collision
        layer; only produced when ``debug`` is set
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nibabel import Binarize, MergeROIs
    from niworkflows.interfaces.reportlets.masks import ROIsPlot
    from templateflow.api import get as get_template

    from ...interfaces.reports import SubcorticalAlignmentReport
    from ...interfaces.workbench import (
        CiftiCreateDenseTimeseries,
        CiftiCreateLabel,
        CiftiDilate,
        CiftiResample,
        CiftiSeparate,
        CiftiSmooth,
        VolumeAllLabelsToROIs,
        VolumeLabelExportTable,
    )

    # reuse saved atlas to atlas transform
    atlas_xfm = load_data('MNIInfant_to_MNI1526NLinAsym.mat')
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['MNIInfant_bold', 'MNIInfant_rois', 'MNI152_rois']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'subcortical_volume',
                'subcortical_labels',
                'subcortical_seg',
                'out_report',
                'out_report_overlap',  # Undefined unless debug
            ]
        ),
        name='outputnode',
    )

    subj_rois = pe.Node(VolumeAllLabelsToROIs(label_map=1), name='subj_rois')
    split_rois = pe.Node(fsl.Split(dimension='t'), name='split_rois')
    atlas_rois = pe.Node(VolumeAllLabelsToROIs(label_map=1), name='atlas_rois')
    split_atlas_rois = pe.Node(fsl.Split(dimension='t'), name='split_atlas_rois')
    atlas_labels = pe.Node(VolumeLabelExportTable(label_map=1), name='atlas_labels')
    parse_labels = pe.Node(
        niu.Function(function=parse_roi_labels, output_names=['structures', 'label_ids']),
        name='parse_labels',
    )

    roi_bbox = pe.MapNode(
        CropToROI(),
        iterfield=['in_file'],
        name='roi_bbox',
    )

    # The following is wrapped in a for-loop, iterating across each roi
    # Instead, we will use MapNodes and iter across the varying inputs
    roi2atlas = pe.MapNode(
        fsl.FLIRT(
            in_matrix_file=atlas_xfm,
            searchr_x=[-20, 20],
            searchr_y=[-20, 20],
            searchr_z=[-20, 20],
            interp='nearestneighbour',
        ),
        name='roi2atlas',
        iterfield=['in_file', 'reference'],
    )
    applyxfm_roi = pe.MapNode(
        fsl.ApplyXFM(interp='spline'),
        iterfield=['reference', 'in_matrix_file'],
        name='applyxfm_roi',
        mem_gb=1,
    )
    bold_mask_roi = pe.MapNode(
        fsl.ApplyMask(),
        iterfield=['in_file', 'mask_file'],
        name='bold_mask_roi',
    )
    mul_roi = pe.MapNode(
        fsl.BinaryMaths(operation='mul'),
        iterfield=['in_file', 'operand_value'],
        name='mul_roi',
    )
    mul_atlas_roi = pe.MapNode(
        fsl.BinaryMaths(operation='mul'),
        iterfield=['in_file', 'operand_value'],
        name='mul_atlas_roi',
    )
    vol_label = pe.MapNode(
        VolumeLabelImport(drop_unused_labels=True),
        iterfield=['in_file'],
        name='vol_label',
    )
    vol_atlas_label = pe.MapNode(
        VolumeLabelImport(drop_unused_labels=True),
        iterfield=['in_file'],
        name='vol_atlas_label',
    )
    create_dtseries = pe.MapNode(
        CiftiCreateDenseTimeseries(),
        iterfield=['volume_data', 'volume_structure_labels'],
        name='create_dtseries',
    )
    create_label = pe.MapNode(
        CiftiCreateLabel(),
        iterfield=['volume_label', 'structure_label_volume'],
        name='create_label',
    )
    dilate = pe.MapNode(
        CiftiDilate(direction='COLUMN', surface_distance=0, volume_distance=10),
        iterfield=['in_file'],
        name='dilate',
    )
    resample = pe.MapNode(
        CiftiResample(
            direction='COLUMN',
            template_direction='COLUMN',
            surface_method='ADAP_BARY_AREA',
            volume_method='CUBIC',
            volume_predilate=10,
        ),
        iterfield=['in_file', 'template'],
        name='resample',
    )
    smooth = pe.MapNode(
        CiftiSmooth(direction='COLUMN', fix_zeros_vol=True, sigma_surf=0, sigma_vol=vol_sigma),
        iterfield=['in_file'],
        name='smooth',
    )
    separate = pe.MapNode(
        CiftiSeparate(direction='COLUMN', volume_all_file='volume_all.nii.gz'),
        iterfield=['in_file'],
        name='separate',
    )

    merge_segs = pe.Node(MergeLabelROIs(), name='merge_segs')

    bg_img = str(
        get_template('MNI152NLin6Asym', resolution=2, desc=None, suffix='T1w', raise_empty=True)
    )
    subcortical_rpt = pe.Node(
        SubcorticalAlignmentReport(anat=bg_img),
        name='subcortical_rpt',
    )
    merge_rois = pe.Node(MergeROIs(), name='merge_rois')

    workflow = Workflow(name=name)
    workflow.connect([
        (inputnode, subj_rois, [('MNIInfant_rois', 'in_file')]),
        (inputnode, atlas_rois, [('MNI152_rois', 'in_file')]),
        (subj_rois, split_rois, [('out_file', 'in_file')]),
        (atlas_rois, split_atlas_rois, [('out_file', 'in_file')]),
        (inputnode, atlas_labels, [('MNI152_rois', 'in_file')]),
        (atlas_labels, parse_labels, [('out_file', 'label_file')]),
        (split_atlas_rois, roi_bbox, [('out_files', 'in_file')]),
        # for loop across ROIs
        (split_rois, roi2atlas, [('out_files', 'in_file')]),
        (roi_bbox, roi2atlas, [('out_file', 'reference')]),
        (inputnode, applyxfm_roi, [('MNIInfant_bold', 'in_file')]),
        (roi_bbox, applyxfm_roi, [('out_file', 'reference')]),
        (roi2atlas, applyxfm_roi, [('out_matrix_file', 'in_matrix_file')]),
        (applyxfm_roi, bold_mask_roi, [('out_file', 'in_file')]),
        (roi2atlas, bold_mask_roi, [('out_file', 'mask_file')]),
        (roi2atlas, mul_roi, [('out_file', 'in_file')]),
        (parse_labels, mul_roi, [('label_ids', 'operand_value')]),
        (split_atlas_rois, mul_atlas_roi, [('out_files', 'in_file')]),
        (parse_labels, mul_atlas_roi, [('label_ids', 'operand_value')]),
        (mul_roi, vol_label, [('out_file', 'in_file')]),
        (atlas_labels, vol_label, [('out_file', 'label_list_file')]),
        (mul_atlas_roi, vol_atlas_label, [('out_file', 'in_file')]),
        (atlas_labels, vol_atlas_label, [('out_file', 'label_list_file')]),
        (bold_mask_roi, create_dtseries, [('out_file', 'volume_data')]),
        (vol_label, create_dtseries, [('out_file', 'volume_structure_labels')]),
        (vol_atlas_label, create_label, [
            ('out_file', 'volume_label'),
            ('out_file', 'structure_label_volume')]),
        (create_dtseries, dilate, [('out_file', 'in_file')]),
        (dilate, resample, [('out_file', 'in_file')]),
        (create_label, resample, [('out_file', 'template')]),
        (resample, smooth, [('out_file', 'in_file')]),
        (smooth, separate, [('out_file', 'in_file')]),

        (mul_roi, merge_segs, [('out_file', 'in_files')]),
        (inputnode, merge_segs, [('MNI152_rois', 'template')]),
        (merge_segs, outputnode, [('out_file', 'subcortical_seg')]),
        (separate, merge_rois, [('volume_all_file', 'in_files')]),
        (merge_rois, outputnode, [('out_file', 'subcortical_volume')]),
        (inputnode, outputnode, [('MNI152_rois', 'subcortical_labels')]),
        (inputnode, subcortical_rpt, [('MNI152_rois', 'reference')]),
        (merge_segs, subcortical_rpt, [('out_file', 'moving')]),
        (subcortical_rpt, outputnode, [('out_report', 'out_report')]),
    ])  # fmt:skip

    if debug:
        binarize_seg = pe.Node(Binarize(thresh_low=0), name='binarize_seg')
        binarize_collision = pe.Node(Binarize(thresh_low=1), name='binarize_collision')
        mrg_seg_collision = pe.Node(niu.Merge(2), name='mrg_seg_collision')
        overlap_rpt = pe.Node(
            ROIsPlot(colors=['b', 'r'], generate_report=True, in_file=bg_img),
            name='overlap_rpt',
        )
        workflow.connect([
            (merge_segs, binarize_seg, [('out_file', 'in_file')]),
            (merge_segs, binarize_collision, [('overlap_file', 'in_file')]),
            (binarize_seg, mrg_seg_collision, [('out_mask', 'in1')]),
            (binarize_collision, mrg_seg_collision, [('out_mask', 'in2')]),
            (mrg_seg_collision, overlap_rpt, [('out', 'in_rois')]),
            (overlap_rpt, outputnode, [('out_report', 'out_report_overlap')]),
        ])  # fmt:skip

    return workflow


def parse_roi_labels(label_file: str):
    """
    Parse a label file composed of one or more sets of:
    <labelname>
    <key> <red> <green> <blue> <alpha>

    Return a list of structure names and label keys.

    >>> structs, ids = parse_roi_labels(datadir / "FreeSurferSubcorticalLabelTableLut.txt")
    >>> structs
    ['ACCUMBENS_LEFT', 'ACCUMBENS_RIGHT', 'AMYGDALA_LEFT', ...]
    >>> ids
    [26, 58, 18, ...]
    """

    with open(label_file) as fp:
        lines = fp.readlines()
    if len(lines) % 2 == 1:
        raise RuntimeError('Label file is incomplete or invalid')
    structs, label_ids = [], []
    for idx, line in enumerate(lines):
        if idx % 2 == 0:
            structs.append(line.strip())
        else:
            label_ids.append(int(line.split(' ', 1)[0]))
    return structs, label_ids
