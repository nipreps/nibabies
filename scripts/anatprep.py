"""Script for testing the initial preprocessing steps of T1w and T2w."""


def init_workflow(bids_path, output_path, participant_label, workdir=None):
    """Create the preprocessing workflow."""
    from nipype.pipeline import engine as pe
    from nibabies.workflows.anatomical.preproc import init_anat_average_wf
    from nibabies.workflows.anatomical.registration import init_coregistration_wf
    from nibabies.workflows.anatomical.brain_extraction import (
        init_infant_brain_extraction_wf,
    )
    from nibabies.workflows.anatomical.outputs import init_coreg_report_wf

    wf = pe.Workflow(name="nibabies_anat")
    for subid in participant_label:
        sub_wf = pe.Workflow(name=f"nibabies_anat_{subid}")
        t1w_files = list(
            (bids_path / f"sub-{subid}" / "anat").glob(f"sub-{subid}*_T1w.nii.gz")
        )

        t2w_files = list(
            (bids_path / f"sub-{subid}" / "anat").glob(f"sub-{subid}*_T2w.nii.gz")
        )

        t1w_ref = init_anat_average_wf(
            num_maps=len(t1w_files), name="t1w_ref", omp_nthreads=8
        )
        t2w_ref = init_anat_average_wf(
            num_maps=len(t2w_files), name="t2w_ref", omp_nthreads=8
        )

        t1w_ref.inputs.inputnode.in_files = [str(f) for f in t1w_files]
        t2w_ref.inputs.inputnode.in_files = [str(f) for f in t2w_files]

        be = init_infant_brain_extraction_wf(omp_nthreads=8, age_months=2)
        cr = init_coregistration_wf(omp_nthreads=8, sloppy=True)

        rpt = init_coreg_report_wf(output_dir=str(output_path.absolute()))
        rpt.inputs.inputnode.source_file = [str(f) for f in t1w_files]

        # fmt:off
        sub_wf.connect([
            (t2w_ref, be, [("outputnode.out_file", "inputnode.in_t2w")]),
            (t1w_ref, cr, [("outputnode.out_file", "inputnode.in_t1w")]),
            (be, cr, [
                ("outputnode.t2w_preproc", "inputnode.in_t2w_preproc"),
                ("outputnode.out_mask", "inputnode.in_mask"),
                ("outputnode.out_probmap", "inputnode.in_probmap"),
            ]),
            (cr, rpt, [
                ("outputnode.t1w_preproc", "inputnode.t1w_preproc"),
                ("outputnode.t2w_preproc", "inputnode.t2w_preproc"),
                ("outputnode.t1w_mask", "inputnode.in_mask"),
            ]),
        ])
        # fmt:on
        wf.add_nodes([sub_wf])

    if workdir:
        wf.base_dir = workdir
    return wf


if __name__ == "__main__":
    from pathlib import Path
    import re
    from argparse import ArgumentParser
    from argparse import RawTextHelpFormatter

    parser = ArgumentParser(
        description="sMRIPrep-infants: Structural MRI PREProcessing workflows",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "bids_dir",
        action="store",
        type=Path,
        help="the root folder of a BIDS valid dataset (sub-XXXXX folders should "
        "be found at the top level in this folder).",
    )
    parser.add_argument(
        "output_dir",
        action="store",
        type=Path,
        help="the output path for the outcomes of preprocessing and visual " "reports",
    )
    parser.add_argument(
        "--participant-label",
        "--participant_label",
        action="store",
        nargs="+",
        help="a space delimited list of participant identifiers or a single "
        "identifier (the sub- prefix can be removed)",
    )

    opts = parser.parse_args()

    participant_label = [
        re.sub(r"^sub-", "", p) for p in opts.participant_label
    ]
    init_workflow(
        opts.bids_dir,
        opts.output_dir,
        participant_label,
        workdir=Path.cwd().absolute() / "workdir"
    ).run()
