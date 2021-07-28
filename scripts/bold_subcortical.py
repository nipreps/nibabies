"""Script for testing the subcortical MNI alignment"""
from pathlib import Path


def init_workflow(bold_file, bold_roi, bold_atlas_roi, vol_sigma):
    from nibabies.workflows.bold.alignment import init_subcortical_mni_alignment_wf

    wf = init_subcortical_mni_alignment_wf(vol_sigma=vol_sigma)
    wf.inputs.inputnode.bold_file = bold_file
    wf.inputs.inputnode.bold_roi = bold_roi
    wf.inputs.inputnode.atlas_roi = bold_atlas_roi

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
        "--vol-sigma",
        type=float,
        default=0.8,
        help="The sigma for the gaussian volume smoothing kernel, in mm",
    )
    parser.add_argument(
        "--nipype-plugin",
        default="MultiProc",
        help="Nipype plugin to run workflow with",
    )
    opts = parser.parse_args()
    wf = init_workflow(
        opts.bold_file.absolute(),
        opts.bold_roi.absolute(),
        opts.bold_atlas_roi.absolute(),
        vol_sigma=opts.vol_sigma,
    )

    wf.config['execution']['crashfile_format'] = 'txt'
    wf.run(plugin=opts.nipype_plugin)
