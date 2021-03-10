# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Writing outputs."""
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from ...interfaces import DerivativesDataSink


def init_coreg_report_wf(*, output_dir, name="coreg_report_wf"):
    """
    Generate and store a report in the right location.

    Parameters
    ----------
    output_dir : :obj:`str`
        Directory in which to save derivatives
    name : :obj:`str`
        Workflow name (default: coreg_report_wf)

    Inputs
    ------
    source_file
        Input reference T1w image
    t1w_preproc
        Preprocessed T1w image.
    t2w_preproc
        Preprocessed T2w image, aligned with the T1w image.
    in_mask
        Brain mask.

    """
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )

    workflow = Workflow(name=name)

    inputfields = [
        "source_file",
        "t1w_preproc",
        "t2w_preproc",
        "in_mask",
    ]
    inputnode = pe.Node(niu.IdentityInterface(fields=inputfields), name="inputnode")
    # Generate reportlets showing spatial normalization
    norm_rpt = pe.Node(
        SimpleBeforeAfter(before_label="T2w", after_label="T1w"),
        name="norm_rpt",
        mem_gb=0.1,
    )

    ds_t1w_t2w_report = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir, space="T2w", suffix="T1w", datatype="figures"
        ),
        name="ds_t1w_t2w_report",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, norm_rpt, [("t2w_preproc", "before"),
                               ("t1w_preproc", "after"),
                               ("in_mask", "wm_seg")]),
        (inputnode, ds_t1w_t2w_report, [("source_file", "source_file")]),
        (norm_rpt, ds_t1w_t2w_report, [("out_report", "in_file")]),
    ])
    # fmt:on

    return workflow
