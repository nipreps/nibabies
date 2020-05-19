"""Command Line Interface."""
from pathlib import Path


def get_parser():
    """Build parser object."""
    from os import cpu_count
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    from .._version import get_versions

    parser = ArgumentParser(
        description="""\
nibabies-be -- Atlas-based brain extraction tool of the \
ANTs package.\
""",
        formatter_class=ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input_image",
        type=Path,
        help="The target image for brain extraction.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="nibabies-be v{}".format(get_versions()["version"]),
    )
    parser.add_argument(
        "--template",
        choices=("MNIInfant",),
        default="MNIInfant",
        help="The TemplateFlow ID of the reference template.",
    )
    parser.add_argument(
        "--cohort",
        type=int,
        choices=range(1,12),
        help="TemplateFlow cohort ID of the reference template"
    )
    parser.add_argument(
        "--omp-nthreads",
        type=int,
        default=cpu_count(),
        help="Number of CPUs available for multithreading processes.",
    )
    parser.add_argument(
        "--nprocs",
        type=int,
        default=cpu_count(),
        help="Number of processes that can be run in parallel.",
    )
    parser.add_argument(
        "-m",
        "--mri-scheme",
        default="T2w",
        choices=("T2w",),
        help="select a particular MRI scheme",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("results").absolute(),
        help="path where intermediate results should be stored",
    )
    parser.add_argument(
        "-w",
        "--work-dir",
        type=Path,
        default=Path("work").absolute(),
        help="path where intermediate results should be stored",
    )
    parser.add_argument(
        "--sloppy",
        dest="debug",
        action="store_true",
        default=False,
        help="Use low-quality tools for speed - TESTING ONLY",
    )
    return parser


def main(argv=None):
    """Entry point."""
    from ..workflows.brain_extraction import init_infant_brain_extraction_wf

    opts = get_parser().parse_args(argv)
    template_specs = {'resolution': 2 if opts.debug else 1}
    if opts.cohort:
        template_specs['cohort'] = opts.cohort
    be = init_infant_brain_extraction_wf(
        debug=opts.debug,
        in_template=opts.template,
        template_specs=template_specs,
        mri_scheme=opts.mri_scheme,
        omp_nthreads=opts.omp_nthreads,
        output_dir=opts.output_dir,
    )
    be.base_dir = opts.work_dir
    be.inputs.inputnode.in_files = opts.input_image

    nipype_plugin = {"plugin": "Linear"}
    if opts.nprocs > 1:
        nipype_plugin["plugin"] = "MultiProc"
        nipype_plugin["plugin_args"] = {
            "nproc": opts.nprocs,
            "raise_insufficient": False,
            "maxtasksperchild": 1,
        }
    be.run(**nipype_plugin)


if __name__ == "__main__":
    raise RuntimeError(
        "Please `pip install` nibabies and use the `nibabies-be` command."
    )