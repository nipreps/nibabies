"""Main runner"""
import logging
from pathlib import Path
import sys


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
        "command",
        choices=('bew', 'bew+surf'),
        help="Specific nibabies commandline workflow",
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
        choices=("MNIInfant", "UNCInfant"),
        default="UNCInfant",
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
        default="T1w",
        choices=("T1w", "T2w"),
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
    parser.add_argument(
        "--age-months",
        dest="age_months",
        type=int,
        help="Age (in months)",
    )
    parser.add_argument(
        "--subject",
        dest="subject_id",
        help="subject ID (if running infant recon-all)"
    )
    return parser


def main(argv=None):
    """Entry point."""
    from nipype import config

    opts = get_parser().parse_args(argv)
    template_specs = {}
    if opts.template == 'MNIInfant':
        template_specs = {'resolution': 2 if opts.debug else 1}

    # specify cohort
    if opts.cohort:
        template_specs['cohort'] = opts.cohort
    elif opts.age_months:
        if opts.age_months <= 2:
            cohort = 1
        elif opts.age_months < 12:
            cohort = 2
        else:
            cohort = 3
        template_specs['cohort'] = cohort

    if opts.command == 'bew':
        from ..workflows.brain_extraction import init_infant_brain_extraction_wf
        wf = init_infant_brain_extraction_wf(
            ants_affine_init=True,
            debug=opts.debug,
            in_template=opts.template,
            template_specs=template_specs,
            mri_scheme=opts.mri_scheme,
            omp_nthreads=opts.omp_nthreads,
            output_dir=opts.output_dir,
        )
    elif opts.command == 'bew+surf':
        from ..workflows.base import init_infant_anat_wf
        wf = init_infant_anat_wf(
            template_name=opts.template,
            template_specs=template_specs,
            age_months=opts.age_months,
            mri_scheme=opts.mri_scheme,
            omp_nthreads=opts.omp_nthreads,
            output_dir=opts.output_dir,
            subject_id=opts.subject_id,
        )
    else:
        print(f"No workflow for command: {opts.command}", file=sys.stderr)
        sys.exit(1)

    # Run the workflow
    wf.inputs.inputnode.in_files = opts.input_image
    wf.base_dir = opts.work_dir
    nipype_plugin = {"plugin": "Linear"}
    if opts.nprocs > 1:
        nipype_plugin["plugin"] = "MultiProc"
        nipype_plugin["plugin_args"] = {
            "nproc": opts.nprocs,
            "raise_insufficient": False,
            "maxtasksperchild": 1,
        }
    wf.base_dir = opts.work_dir

    logging.getLogger('nipype.interface').setLevel("DEBUG")
    wf.run(**nipype_plugin)

if __name__ == "__main__":
    raise RuntimeError(
        "Please `pip install` this and run via the commandline interfaces, `nibabies <command>`"
    )