"""
Subcortical alignment into MNI space
"""

from nipype.interfaces import fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from nibabies.data import load as load_data
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


SUBCORTICAL_GROUPS = {
    'basal_ganglia_left': ('ACCUMBENS_LEFT', 'CAUDATE_LEFT', 'PALLIDUM_LEFT', 'PUTAMEN_LEFT'),
    'basal_ganglia_right': ('ACCUMBENS_RIGHT', 'CAUDATE_RIGHT', 'PALLIDUM_RIGHT', 'PUTAMEN_RIGHT'),
    'mtl_left': ('HIPPOCAMPUS_LEFT', 'AMYGDALA_LEFT'),
    'mtl_right': ('HIPPOCAMPUS_RIGHT', 'AMYGDALA_RIGHT'),
    'cerebellum_left': ('CEREBELLUM_LEFT',),
    'cerebellum_right': ('CEREBELLUM_RIGHT',),
    'thalamus_left': ('THALAMUS_LEFT',),
    'thalamus_right': ('THALAMUS_RIGHT',),
    'diencephalon_brainstem': (
        'DIENCEPHALON_VENTRAL_LEFT',
        'DIENCEPHALON_VENTRAL_RIGHT',
        'BRAIN_STEM',
    ),
}
"""Neighboring CIFTI subcortical structures that share a single affine alignment."""


def init_subcortical_mni_alignment_wf(
    *,
    vol_sigma: float = 0.8,
    groups: dict[str, tuple[str, ...]] | None = None,
    name: str = 'subcortical_mni_alignment_wf',
):
    """
    Align groups of neighboring subcortical structures into MNI space.

    Adapted from the DCAN infant pipeline:
    https://github.com/DCAN-Labs/dcan-infant-pipeline/blob\
    /master/fMRISurface/scripts/SubcorticalAlign_ROIs.sh

    Differences from the DCAN implementation:

    * Neighboring structures are aligned together, as groups rather than one structure at a time.
    * The aligned groups are combined by summing their volumes, rather than through a dense
      CIFTI template, to reduce memory usage.

    Parameters
    ----------
    name : :obj:`str`
        Name of the workflow
    vol_sigma : :obj:`float`
        The sigma for the gaussian volume smoothing kernel, in mm
    groups : :obj:`dict` or None
        Mapping of group name to CIFTI structure names, covering every structure exactly
        once. Defaults to :data:`SUBCORTICAL_GROUPS`.

    Inputs
    ------
    MNIInfant_bold : :obj:`str`
        BOLD file in MNI Infant space
    MNIInfant_rois : :obj:`str`
        Subject subcortical label volume in MNI Infant space
    MNI152_rois : :obj:`str`
        HCP subcortical atlas label volume in MNI152NLin6Asym space

    Outputs
    -------
    subcortical_volume : :obj:`str`
        BOLD volume on the atlas grid, with each structure group aligned separately
    subcortical_labels : :obj:`str`
        Atlas label volume defining the subcortical grayordinates (``MNI152_rois``)
    """
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from niworkflows.interfaces.nibabel import MapLabels, MergeROIs
    from niworkflows.interfaces.utility import KeySelect

    from ...interfaces.workbench import (
        CiftiCreateDenseTimeseries,
        CiftiCreateLabel,
        CiftiDilate,
        CiftiResample,
        CiftiSeparate,
        CiftiSmooth,
        VolumeLabelExportTable,
    )

    if groups is None:
        groups = SUBCORTICAL_GROUPS
    mappings = group_label_mappings(groups, load_data('FreeSurferSubcorticalLabelTableLut.txt'))

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['MNIInfant_bold', 'MNIInfant_rois', 'MNI152_rois']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['subcortical_volume', 'subcortical_labels']),
        name='outputnode',
    )

    select_subcortical_group = pe.Node(
        KeySelect(fields=['mappings'], keys=list(groups), mappings=mappings),
        name='select_subcortical_group',
        run_without_submitting=True,
    )
    select_subcortical_group.iterables = ('key', list(groups))
    subj_group = pe.Node(MapLabels(), name='subj_group')
    atlas_group = pe.Node(MapLabels(), name='atlas_group')
    atlas_labels = pe.Node(VolumeLabelExportTable(label_map=1), name='atlas_labels')

    roi2atlas = pe.Node(
        fsl.FLIRT(
            searchr_x=[-20, 20],
            searchr_y=[-20, 20],
            searchr_z=[-20, 20],
            interp='nearestneighbour',
        ),
        name='roi2atlas',
    )
    applyxfm_roi = pe.Node(
        fsl.ApplyXFM(interp='spline'),
        name='applyxfm_roi',
        mem_gb=4,
    )
    vol_label = pe.Node(VolumeLabelImport(drop_unused_labels=True), name='vol_label')
    vol_atlas_label = pe.Node(VolumeLabelImport(drop_unused_labels=True), name='vol_atlas_label')
    create_dtseries = pe.Node(CiftiCreateDenseTimeseries(), name='create_dtseries')
    create_label = pe.Node(CiftiCreateLabel(), name='create_label')
    dilate = pe.Node(
        CiftiDilate(direction='COLUMN', surface_distance=0, volume_distance=10), name='dilate'
    )
    resample = pe.Node(
        CiftiResample(
            direction='COLUMN',
            template_direction='COLUMN',
            surface_method='ADAP_BARY_AREA',
            volume_method='CUBIC',
            volume_predilate=10,
        ),
        name='resample',
    )
    smooth = pe.Node(
        CiftiSmooth(direction='COLUMN', fix_zeros_vol=True, sigma_surf=0, sigma_vol=vol_sigma),
        name='smooth',
    )
    separate = pe.Node(
        CiftiSeparate(direction='COLUMN', volume_all_file='volume_all.nii.gz'), name='separate'
    )
    join_subcortical_groups = pe.JoinNode(
        niu.IdentityInterface(fields=['volumes']),
        joinsource='select_subcortical_group',
        joinfield=['volumes'],
        name='join_subcortical_groups',
    )
    merge_rois = pe.Node(MergeROIs(), name='merge_rois')

    workflow = Workflow(name=name)
    workflow.connect([
        (inputnode, subj_group, [('MNIInfant_rois', 'in_file')]),
        (inputnode, atlas_group, [('MNI152_rois', 'in_file')]),
        (select_subcortical_group, subj_group, [('mappings', 'mappings')]),
        (select_subcortical_group, atlas_group, [('mappings', 'mappings')]),
        (inputnode, atlas_labels, [('MNI152_rois', 'in_file')]),
        # per structure group (select_subcortical_group iterables)
        (subj_group, roi2atlas, [('out_file', 'in_file')]),
        (atlas_group, roi2atlas, [('out_file', 'reference')]),
        (inputnode, applyxfm_roi, [('MNIInfant_bold', 'in_file')]),
        (atlas_group, applyxfm_roi, [('out_file', 'reference')]),
        (roi2atlas, applyxfm_roi, [('out_matrix_file', 'in_matrix_file')]),
        (roi2atlas, vol_label, [('out_file', 'in_file')]),
        (atlas_labels, vol_label, [('out_file', 'label_list_file')]),
        (atlas_group, vol_atlas_label, [('out_file', 'in_file')]),
        (atlas_labels, vol_atlas_label, [('out_file', 'label_list_file')]),
        (applyxfm_roi, create_dtseries, [('out_file', 'volume_data')]),
        (vol_label, create_dtseries, [('out_file', 'volume_structure_labels')]),
        (vol_atlas_label, create_label, [
            ('out_file', 'volume_label'),
            ('out_file', 'structure_label_volume')]),
        (create_dtseries, dilate, [('out_file', 'in_file')]),
        (dilate, resample, [('out_file', 'in_file')]),
        (create_label, resample, [('out_file', 'template')]),
        (resample, smooth, [('out_file', 'in_file')]),
        (smooth, separate, [('out_file', 'in_file')]),
        (separate, join_subcortical_groups, [('volume_all_file', 'volumes')]),
        (join_subcortical_groups, merge_rois, [('volumes', 'in_files')]),
        (merge_rois, outputnode, [('out_file', 'subcortical_volume')]),
        (inputnode, outputnode, [('MNI152_rois', 'subcortical_labels')]),
    ])  # fmt:skip
    return workflow


def parse_roi_labels(label_file: str) -> dict[str, int]:
    """
    Parse a label file composed of one or more sets of:
    <labelname>
    <key> <red> <green> <blue> <alpha>

    Return a mapping of structure names to label keys, in file order.

    >>> labels = parse_roi_labels(datadir / "FreeSurferSubcorticalLabelTableLut.txt")
    >>> list(labels.keys())
    ['ACCUMBENS_LEFT', 'ACCUMBENS_RIGHT', 'AMYGDALA_LEFT', ...]
    >>> list(labels.values())
    [26, 58, 18, ...]
    """

    with open(label_file) as fp:
        lines = fp.readlines()
    if len(lines) % 2 == 1:
        raise RuntimeError('Label file is incomplete or invalid')
    return {
        name.strip(): int(entry.split(' ', 1)[0])
        for name, entry in zip(lines[::2], lines[1::2], strict=True)
    }


def group_label_mappings(groups: dict[str, tuple[str, ...]], label_file: str) -> list[dict]:
    """
    Build one label remapping per group, zeroing every label outside that group.

    Raises :obj:`ValueError` unless ``groups`` covers every structure in ``label_file``
    exactly once.

    >>> mappings = group_label_mappings(
    ...     SUBCORTICAL_GROUPS, datadir / "FreeSurferSubcorticalLabelTableLut.txt"
    ... )
    >>> len(mappings)
    9
    >>> sorted(mappings[0])  # basal_ganglia_left keeps 11, 12, 13, 26
    [8, 10, 16, 17, 18, 28, 47, 49, 50, 51, 52, 53, 54, 58, 60]
    """
    labels = parse_roi_labels(label_file)
    grouped = [val for values in groups.values() for val in values]

    unknown = sorted(set(grouped) - set(labels))
    duped = sorted({s for s in grouped if grouped.count(s) > 1})
    missing = sorted(set(labels) - set(grouped))
    if unknown or duped or missing:
        raise ValueError(
            'Subcortical groups must cover each structure exactly once '
            f'(unknown: {unknown}, duplicated: {duped}, missing: {missing})'
        )

    return [
        {key: 0 for struct, key in labels.items() if struct not in members}
        for members in groups.values()
    ]
