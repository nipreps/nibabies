"""
Subcortical alignment into MNI space
"""

from nipype.interfaces.fsl.maths import MultiImageMaths
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from nipype.interfaces.workbench.cifti import CiftiSmooth
from ...interfaces.workbench import (
    CiftiCreateDenseFromTemplate,
    CiftiCreateDenseTimeseries,
    CiftiCreateLabel,
    CiftiDilate,
    CiftiResample,
    CiftiSeparate,
    VolumeAffineResample,
    VolumeAllLabelsToROIs,
    VolumeLabelExportTable,
    VolumeLabelImport,
)


def gen_subcortical_alignment_wf(repetition_time, rois, name='subcortical_alignment_wf'):
    """
    Align individual subcortical structures into MNI space.

    This is a nipype workflow port of the DCAN infant pipeline.
    https://github.com/DCAN-Labs/dcan-infant-pipeline/blob/247e19e5441cc814cea2f23720caeeb6c6aeadf8/fMRISurface/scripts/SubcorticalAlign_ROIs.sh

    Parameters
    ----------
    repetition_time : :obj:`int`
        BOLD file's TR
    rois : :obj:`list`
        ROIs to align
    name : :obj:`str`
        Name of the workflow

    Inputs
    ------
    label_file : :obj:`str`
        File containing ROI labels
    std_xfm : :obj:`str`
        File containing transform to the standard (MNI) space

    Outputs
    -------

    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow


    # inputs
    # ${ROIFolder}/Atlas_ROIs.${GrayordinatesResolution}.nii.gz
    # ${ROIFolder}/ROIs.${GrayordinatesResolution}.nii.gz
    # $HCPPIPEDIR_Templates/InfMNI_2AdultMNI_Step2.mat
    inputnode = pe.Node(
        niu.IdentityInterface(fields=["bold_roi", "atlas_roi", "atlas_xfm"]),
        name="inputnode",
    )
    # inputnode = pe.Node(
    #     niu.IdentityInterface(fields=["bold_file", "rois", "label_file", "std_xfm"]),
    #     name='inputnode',
    # )
    outputnode = pe.Node(niu.IdentityInterface(fields=[]), name='outputnode')


    #generate altas-roi space fMRI cifti for subcortical data
    # flirt -in ${VolumefMRI}.nii.gz -ref ${ROIFolder}/Atlas_ROIs.${GrayordinatesResolution}.nii.gz -applyxfm -init $HCPPIPEDIR_Templates/InfMNI_2AdultMNI_Step2.mat -out ${VolumefMRI}_2MNI.nii.gz
    applyxfm_atlas = pe.Node(fsl.FLIRT(), name="applyxfm_atlas")

    # ${CARET7DIR}/wb_command -volume-affine-resample ${ROIFolder}/ROIs.${GrayordinatesResolution}.nii.gz $HCPPIPEDIR_Templates/InfMNI_2AdultMNI_Step2.mat ${VolumefMRI}_2MNI.nii.gz ENCLOSING_VOXEL ${ResultsFolder}/ROIs.${GrayordinatesResolution}.nii.gz -flirt ${ROIFolder}/ROIs.${GrayordinatesResolution}.nii.gz ${VolumefMRI}_2MNI.nii.gz
    volume_resample = pe.Node(VolumeAffineResample(), name="volume_resample")
    # ${CARET7DIR}/wb_command -cifti-create-dense-timeseries ${WD}/${NameOffMRI}_temp_orig_atlas.dtseries.nii -volume ${VolumefMRI}_2MNI.nii.gz ${ROIFolder}/Atlas_ROIs.${GrayordinatesResolution}.nii.gz
    create_dense = pe.Node(CiftiCreateDenseTimeseries(), name="create_dense")

    # #splitting atlas and subject volume label files into individual ROI files for registration
    # ${CARET7DIR}/wb_command -volume-all-labels-to-rois ${ROIFolder}/ROIs.${GrayordinatesResolution}.nii.gz 1 ${WD}/sub_allroi.nii.gz
    subj_rois = pe.Node(VolumeAllLabelsToROIs(label_map=1), name="subj_rois")
    # fslsplit ${WD}/sub_allroi.nii.gz sub_ROI -t
    split_rois = pe.Node(fsl.Split(dimension="t"), name="split_rois")
    # ${CARET7DIR}/wb_command -volume-all-labels-to-rois ${ROIFolder}/Atlas_ROIs.${GrayordinatesResolution}.nii.gz 1 ${WD}/atl_allroi.nii.gz
    atlas_rois = pe.Node(VolumeAllLabelsToROIs(label_map=1), name="atlas_rois")
    # fslsplit ${WD}/atl_allroi.nii.gz atl_ROI -t
    split_atlas_rois = pe.Node(fsl.Split(dimension="t"), name="split_atlas_rois")

    # #exporting table to generate independent label files
    # ${CARET7DIR}/wb_command -volume-label-export-table ${ROIFolder}/Atlas_ROIs.${GrayordinatesResolution}.nii.gz 1 ${WD}/labelfile.txt
    atlas_labels = pe.Node(VolumeLabelExportTable(label_map=1), name="atlas_labels")

    # #initialize workbench command for creating dense time series.
    # DTSCommand="${CARET7DIR}/wb_command -cifti-create-dense-from-template ${WD}/${NameOffMRI}_temp_orig_atlas.dtseries.nii ${WD}/${NameOffMRI}_temp_atlas.dtseries.nii -series ${TR} 0.0 "
    # #initialize command to make sub2atl_label_ROI.2.nii.gz
    # Sub2AtlCmd=""

    ### The following is wrapped in a for-loop, iterating across each loop

    #perform linear mapping from subject to atlas ROI
    # flirt -in ${WD}/sub_${ROIname} -ref ${WD}/atl_${ROIname} -searchrx -20 20 -searchry -20 20 -searchrz -20 20 -o sub2atl_${ROIname} -interp nearestneighbour -omat ${WD}/sub2atl_${ROInum}.mat
    roi2atlas = pe.MapNode(
        fsl.FLIRT(
            searchr_x=[-20, 20],
            searchr_y=[-20, 20],
            searchr_z=[-20, 20],
            interp="nearestneighbour",
        ),
        name="roi2atlas",
    )

    # flirt -in ${VolumefMRI} -ref ${WD}/atl_${ROIname} -applyxfm -init ${WD}/sub2atl_${ROInum}.mat -o sub2atl_vol_${ROIname} -interp spline
    applyxfm = pe.MapNode(fsl.ApplyXFM(interp="spline"), name='applyxfm')

    # #masking BOLD volumetric data to ROI only
    # fslmaths sub2atl_vol_${ROIname} -mas sub2atl_${ROIname} sub2atl_vol_masked_${ROIname}
    bold_mask_roi = pe.MapNode(fsl.ApplyMask(), name='bold_mask_roi')

    # #extract information from label file (values on even numbered line preceded by label on odd)
    parse_labels = pe.Node(
        niu.Function(function=parse_roi_labels, output_names=["structures", "label_id"]), name="parse_labels")

    # #multiply ROI volume by value -- needed for creating a volume label file
    # fslmaths ${WD}/sub2atl_${ROIname} -mul $roi_value ${WD}/sub2atl_label_${ROIname}
    # fslmaths ${WD}/atl_${ROIname} -mul $roi_value ${WD}/atl_label_${ROIname}
    mul_roi = pe.MapNode(fsl.BinaryMaths(operation="mul"), name='mul_roi')
    mul_atlas_roi = pe.MapNode(fsl.BinaryMaths(operation="mul"), name='mul_atlas_roi')

    # #use wb_command to create the volume label file, needed for creating the individual dtseries
    # ${CARET7DIR}/wb_command -volume-label-import ${WD}/sub2atl_label_${ROIname} ${WD}/labelfile.txt sub2atl_vol_label_${ROIname} -drop-unused-labels
    # ${CARET7DIR}/wb_command -volume-label-import ${WD}/atl_label_${ROIname} ${WD}/labelfile.txt atl_vol_label_${ROIname} -drop-unused-labels
    vol_label = pe.MapNode(VolumeLabelImport(drop_unused_labels=True), name='vol_label')
    vol_atlas_label = pe.MapNode(VolumeLabelImport(drop_unused_labels=True), name='vol_atlas_label')

    # #create the individual dtseries
    # ${CARET7DIR}/wb_command -cifti-create-dense-timeseries ${WD}/${NameOffMRI}_temp_subject_${ROInum}.dtseries.nii -volume sub2atl_vol_masked_${ROIname} ${WD}/sub2atl_vol_label_${ROIname}
    # # Maybe here, too.
    create_dtseries = pe.MapNode(CiftiCreateDenseTimeseries(), name='create_dtseries')

    # #create the cifti label file from the volume label file (why?????)
    # #${CARET7DIR}/wb_command -cifti-create-label ${WD}/sub_${NameOffMRI}_temp_template_${ROInum}.dlabel.nii -volume ${WD}/sub2atl_vol_label_${ROIname} ${WD}/sub2atl_vol_label_${ROIname}
    # ${CARET7DIR}/wb_command -cifti-create-label ${WD}/atl_${NameOffMRI}_temp_template_${ROInum}.dlabel.nii -volume ${WD}/atl_vol_label_${ROIname} ${WD}/atl_vol_label_${ROIname}
    create_label = pe.MapNode(CiftiCreateLabel(), name='create_label')

    # #dilate the timeseries
    # ${CARET7DIR}/wb_command -cifti-dilate ${WD}/${NameOffMRI}_temp_subject_${ROInum}.dtseries.nii COLUMN 0 10 ${WD}/${NameOffMRI}_temp_subject_${ROInum}_dilate.dtseries.nii
    dilate = pe.MapNode(
        CiftiDilate(direction="COLUMN", surface_distance=0, volume_distance=10),
        name="dilate"
    )

    # #perform resampling - resample into Atlas space, not subject space.
    # ${CARET7DIR}/wb_command -cifti-resample ${WD}/${NameOffMRI}_temp_subject_${ROInum}_dilate.dtseries.nii COLUMN ${WD}/atl_${NameOffMRI}_temp_template_${ROInum}.dlabel.nii COLUMN ADAP_BARY_AREA CUBIC ${WD}/${NameOffMRI}_temp_atlas_${ROInum}.dtseries.nii -volume-predilate 10
    resample = pe.MapNode(
        CiftiResample(
            direction="COLUMN",
            template_direction="COLUMN",
            surface_method="ADAP_BARY_AREA",
            volume_method="CUBIC",
            volume_predilate=10,
        ),
        name='resample',
    )

    # #perform smoothing
    # ${CARET7DIR}/wb_command -cifti-smoothing ${WD}/${NameOffMRI}_temp_atlas_${ROInum}.dtseries.nii 0 ${Sigma} COLUMN ${WD}/${NameOffMRI}_temp_subject_dilate_resample_smooth_${ROInum}.dtseries.nii -fix-zeros-volume
    smooth = pe.MapNode(CiftiSmooth(direction="COLUMN", fix_zeros_vol=True), name="smooth")

    # #split back into a volumetric timeseries file
    # ${CARET7DIR}/wb_command -cifti-separate ${WD}/${NameOffMRI}_temp_subject_dilate_resample_smooth_${ROInum}.dtseries.nii COLUMN -volume-all ${ResultsFolder}/${NameOffMRI}_${ROInum}.nii.gz
    separate = pe.MapNode(CiftiSeparate(direction="COLUMN", volume_all=True), name="separate")

    # #add input to wb_command to grow iteratively
    # DTSCommand="${CARET7DIR}/wb_command -cifti-create-dense-from-template ${WD}/${NameOffMRI}_temp_orig_atlas.dtseries.nii ${WD}/${NameOffMRI}_temp_atlas.dtseries.nii -series ${TR} 0.0"
    # DTSCommand="${DTSCommand} -volume ${roi_name} ${ResultsFolder}/${NameOffMRI}_${ROInum}.nii.gz"
    create_dtseries = pe.MapNode(
        CiftiCreateDenseFromTemplate(series=True, series_step=repetition_time, series_start=0),
        name='create_dtseries',
    )

    # fslmaths <file> -add <sub2atl_label_
    # Sub2AtlCmd="${Sub2AtlCmd}-add ${WD}/sub2atl_label_${ROIname} "
    operations = "-add %s " * len(rois) - 1
    agg_rois = pe.MapNode(fsl.MultiImageMaths(op_string=operations.strip()), name='agg_rois')

    wf = Workflow(name=name)
    wf.connect([
        # flirt output: sub2atl_${ROIname}
        # applyxfm output: sub2atl_vol_${ROIname}
        # bold_mask_roi output: sub2atl_vol_masked_${ROIname}
        (inputnode, roi2atlas, [("rois", "in_file")]),
        (inputnode, applyxfm, [("bold_file", "in_file")]),
        (inputnode, parse_labels, [("label_file", "label_file")]),
        (applyxfm, bold_mask_roi, [("out_file", "in_file")]),
        (bold_mask_roi, create_dtseries, [("out_file", "volume_data")]),
        (roi2atlas, mul_roi, [("out_file", "in_file")]),
    ])

def parse_roi_labels(label_file):
    """
    Parse a label file composed of one or more sets of:
    <labelname>
    <key> <red> <green> <blue> <alpha>

    Return a list of structure names and label keys.

    Example
    -------
    CEREBELLUM_LEFT
    8 230 148 34 255
    THALAMUS_LEFT
    10 0 118 14 255

    TODO: Add unit test
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
            label_ids.append(line.split(' ', 1)[0])
    return structs, label_ids


def _format_volume_rois(rois, structs):
    """
    Format volume arguments for CiftiCreateDenseFromTemplate.
    """

    return [(struct, roi) for struct, roi in zip(structs, rois)]
