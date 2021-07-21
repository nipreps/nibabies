"""Script for testing the subcortical MNI alignment"""
from pathlib import Path


def init_workflow(bold_file, bold_roi, bold_atlas_roi, atlas_xfm, TR=None):
    from nibabies import config
    from nibabies.workflows.bold.alignment import init_subcortical_mni_alignment_wf

    if TR is None:
        # guess TR from header
        import nibabel as nb
        img = nb.load(bold_file)
        assert len(img.shape) > 3, "Not a 4D file"
        TR = img.header['pixdim'][4]

    wf = init_subcortical_mni_alignment_wf(repetition_time=TR)
    wf.inputs.inputnode.bold_file = bold_file
    wf.inputs.inputnode.bold_roi = bold_roi
    wf.inputs.inputnode.atlas_roi = bold_atlas_roi
    wf.inputs.inputnode.atlas_xfm = atlas_xfm

    wf.base_dir = Path('workdir').absolute()
    return wf


if __name__ == "__main__":
    from argparse import ArgumentParser, RawTextHelpFormatter

    parser = ArgumentParser(
        description="DCAN subcortical MNI alignment",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "bold_file",
        type=Path,
        help="the input BOLD file",
    )
    parser.add_argument(
        "bold_roi",
        type=Path,
        help="segmentations in BOLD space",
    )
    parser.add_argument(
        "bold_atlas_roi",
        type=Path,
        help="segmentations in ROI space, unrefined",
    )
    parser.add_argument(
        "atlas_xfm",
        type=Path,
        help="transformation of input BOLD file to MNI space",
    )
    parser.add_argument(
        "--tr",
        type=float,
        help="BOLD repetition time. If not provided, NIfTI header information is used",
    )
    opts = parser.parse_args()
    init_workflow(
        opts.bold_file.absolute(),
        opts.bold_roi.absolute(),
        opts.bold_atlas_roi.absolute(),
        opts.atlas_xfm.absolute(),
        TR=opts.tr,
    ).run(plugin="MultiProc")
