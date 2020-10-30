from pkg_resources import resource_filename as pkgr_fn
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype.interfaces.ants.segmentation import JointFusion
from niworkflows.interfaces.fixes import (
    FixHeaderRegistration as Registration,
    FixHeaderApplyTransforms as ApplyTransforms,
)
from smriprep.utils.misc import apply_lut as _apply_bids_lut
from smriprep.workflows.anatomical import _aseg_to_three, _split_segments

from ...config import DEFAULT_MEMORY_MIN_GB


def init_anat_seg_wf(
    age_months=None,
    anat_modality="T1w",
    template_dir=None,
    sloppy=False,
    omp_nthreads=None,
    name="anat_seg_wf",
):
    """Calculate segmentation from collection of OHSU atlases"""

    if template_dir is None:  # TODO: set default
        raise NotImplementedError(
            "No default has been set yet. Please pass segmentations."
        )
    if anat_modality != "T1w":
        raise NotImplementedError(
            "Only T1w images are currently accepted for ANTs LabelFusion."
        )

    wf = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=["anat_brain"]), name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["anat_aseg", "anat_dseg", "anat_tpms"]),
        name="outputnode",
    )

    tmpl_anats, tmpl_segs = _parse_segmentation_atlases(anat_modality, template_dir)

    # register to templates
    ants_params = "testing" if sloppy else "precise"
    # Register to each subject space
    norm = pe.MapNode(
        Registration(
            from_file=pkgr_fn(
                "niworkflows.data", f"antsBrainExtraction_{ants_params}.json"
            )
        ),
        name="norm",
        iterfield=["fixed_image", "moving_image"],
        n_procs=omp_nthreads,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )
    norm.inputs.moving_image = tmpl_anats
    norm.inputs.float = True

    apply_atlas = pe.MapNode(
        ApplyTransforms(
            dimension=3,
            interpolation="NearestNeighbor",
            float=True,
        ),
        iterfield=["transforms", "input_image"],
        name="apply_atlas",
    )
    apply_atlas.inputs.input_image = tmpl_anats

    apply_seg = pe.MapNode(
        ApplyTransforms(dimension=3, interpolation="MultiLabel"),  # NearestNeighbor?
        name="apply_seg",
        iterfield=["transforms", "input_image"],
    )
    apply_seg.inputs.input_image = tmpl_segs

    jointfusion = pe.Node(
        JointFusion(
            dimension=3,
            out_label_fusion="fusion_labels.nii.gz",
        ),
        name="jointfusion",
    )

    # Convert FreeSurfer aseg to three-tissue segmentation
    lut_anat_dseg = pe.Node(
        niu.Function(function=_apply_bids_lut), name="lut_anat_dseg"
    )
    lut_anat_dseg.inputs.lut = _aseg_to_three

    # split each tissue into individual masks
    split_seg = pe.Node(niu.Function(function=_split_segments), name="split_seg")

    # fmt: off
    wf.connect([
        (inputnode, norm, [('anat_brain', 'fixed_image')]),
        (norm, apply_atlas, [('forward_transforms', 'transforms')]),
        (inputnode, apply_atlas, [('anat_brain', 'reference_image')]),
        (norm, apply_seg, [('forward_transforms', 'transforms')]),
        (inputnode, apply_seg, [('anat_brain', 'reference_image')]),
        (inputnode, jointfusion, [(('anat_brain', _to_list), 'target_image')]),
        (apply_atlas, jointfusion, [('output_image', 'atlas_image')]),
        (apply_seg, jointfusion, [('output_image', 'atlas_segmentation_image')]),
        (jointfusion, outputnode, [('out_label_fusion', 'anat_aseg')]),
        (jointfusion, lut_anat_dseg, [('out_label_fusion', 'in_dseg')]),
        (lut_anat_dseg, outputnode, [('out', 'anat_dseg')]),
        (lut_anat_dseg, split_seg, [('out', 'in_file')]),
        (split_seg, outputnode, [('out', 'anat_tpms')]),
    ])
    # fmt: on

    return wf


def _parse_segmentation_atlases(anat_modality, template_dir):
    """
    Parse segmentation templates directory for anatomical and segmentation files.

    This is currently hardcoded to match DCAN lab templates.
    Will need to rethink standardization for more general cases.
    """
    from pathlib import Path

    anats, segs = [], []

    for f in Path(template_dir).glob("**/*.nii*"):
        if "Segmentation" in f.name:
            segs.append(str(f.absolute()))
        elif anat_modality in f.name:
            anats.append(str(f.absolute()))

    assert anats
    assert segs
    # there should matching files per template
    assert len(anats) == len(segs)

    return sorted(anats), sorted(segs)


def _to_list(x):
    return [x]
