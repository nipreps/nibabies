from pkg_resources import resource_filename as pkgr_fn
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from nipype.interfaces.ants.segmentation import JointFusion
from niworkflows.interfaces.fixes import (
    FixHeaderRegistration as Registration,
    FixHeaderApplyTransforms as ApplyTransforms,
)
from smriprep.utils.misc import apply_lut as _apply_bids_lut
from smriprep.workflows.anatomical import _aseg_to_three, _split_segments, _probseg_fast2bids

from ...config import DEFAULT_MEMORY_MIN_GB


def init_anat_seg_wf(
    age_months=None,
    anat_modality="T1w",
    template_dir=None,
    sloppy=False,
    omp_nthreads=1,
    name="anat_seg_wf",
):
    """
    Calculate brain tissue segmentations from either:
     A) a collection of manually segmented templates
     B) FSL's FAST
    """

    if anat_modality != "T1w":
        raise NotImplementedError(
            "Only T1w images are currently accepted for the segmentation workflow."
        )

    wf = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=["anat_brain"]), name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["anat_aseg", "anat_dseg", "anat_tpms"]),
        name="outputnode",
    )

    # Coerce segmentation labels to BIDS
    lut_anat_dseg = pe.Node(niu.Function(function=_apply_bids_lut), name="lut_anat_dseg")

    if template_dir is None:
        # Use FSL FAST for segmentations
        anat_dseg = pe.Node(
            fsl.FAST(segments=True, no_bias=True, probability_maps=True),
            name="anat_dseg",
            mem_gb=3,
        )
        lut_anat_dseg.inputs.lut = (0, 3, 1, 2)  # Maps: 0 -> 0, 3 -> 1, 1 -> 2, 2 -> 3.
        fast2bids = pe.Node(
            niu.Function(function=_probseg_fast2bids),
            name="fast2bids",
            run_without_submitting=True,
        )

        wf.connect(
            [
                (
                    inputnode,
                    anat_dseg,
                    [
                        ("anat_brain", "in_files"),
                    ],
                ),
                (
                    anat_dseg,
                    lut_anat_dseg,
                    [
                        ("partial_volume_map", "in_dseg"),
                    ],
                ),
                (
                    lut_anat_dseg,
                    outputnode,
                    [
                        ("out", "anat_dseg"),
                    ],
                ),
                (
                    anat_dseg,
                    fast2bids,
                    [
                        ("partial_volume_files", "inlist"),
                    ],
                ),
                (
                    fast2bids,
                    outputnode,
                    [
                        ("out", "anat_tpms"),
                    ],
                ),
            ]
        )
        return wf

    # Otherwise, register to templates and run ANTs JointFusion
    lut_anat_dseg.inputs.lut = _aseg_to_three()
    tmpl_anats, tmpl_segs = _parse_segmentation_atlases(anat_modality, template_dir)

    # register to templates
    ants_params = "testing" if sloppy else "precise"
    # Register to each subject space
    norm = pe.MapNode(
        Registration(
            from_file=pkgr_fn("niworkflows.data", f"antsBrainExtraction_{ants_params}.json")
        ),
        name="norm",
        iterfield=["moving_image"],
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
            num_threads=omp_nthreads,
        ),
        name="jointfusion",
    )

    jf_label = pe.Node(niu.Function(function=_to_dtype), name="jf_label")

    # split each tissue into individual masks
    split_seg = pe.Node(niu.Function(function=_split_segments), name="split_seg")
    to_list = pe.Node(niu.Function(function=_to_list), name="to_list")

    # fmt: off
    wf.connect([
        (inputnode, norm, [('anat_brain', 'fixed_image')]),
        (norm, apply_atlas, [('forward_transforms', 'transforms')]),
        (inputnode, apply_atlas, [('anat_brain', 'reference_image')]),
        (norm, apply_seg, [('forward_transforms', 'transforms')]),
        (inputnode, apply_seg, [('anat_brain', 'reference_image')]),
        (inputnode, to_list, [('anat_brain', 'in_file')]),
        (to_list, jointfusion, [('out', 'target_image')]),
        (apply_atlas, jointfusion, [('output_image', 'atlas_image')]),
        (apply_seg, jointfusion, [('output_image', 'atlas_segmentation_image')]),
        (jointfusion, jf_label, [('out_label_fusion', 'in_file')]),
        (jf_label, outputnode, [('out', 'anat_aseg')]),
        (jf_label, lut_anat_dseg, [('out', 'in_dseg')]),
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


def _to_list(in_file):
    return [in_file]


def _to_dtype(in_file, dtype="uint8"):
    """
    Freesurfer's ``mri_convert`` complains about unsigned 32-bit integers.
    Since we may plan using the JLF segmentation within ``infant_recon_all``,
    better to make this change now.
    """
    import nibabel as nb
    import numpy as np
    from pathlib import Path

    img = nb.load(in_file)
    out_file = Path(f"labels{''.join(Path(in_file).suffixes)}").absolute()

    new_data = np.asanyarray(img.get_fdata(), dtype=dtype)
    img.set_data_dtype(dtype)
    img.__class__(new_data, img.affine, img.header).to_filename(out_file)
    return str(out_file)
