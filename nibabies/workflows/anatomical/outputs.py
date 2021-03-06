# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Writing outputs."""
from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from niworkflows.engine.workflows import LiterateWorkflow as Workflow

from ..interfaces import DerivativesDataSink


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
        Preprocessed T2w image.
    t2w_mask
        Brain mask in T2w space.
    t1w2t2w_xfm
        Mapping to resample T1w-space data into T2w-space.

    """
    from niworkflows.interfaces.fixes import (
        FixHeaderApplyTransforms as ApplyTransforms,
    )
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )

    workflow = Workflow(name=name)

    inputfields = [
        "source_file",
        "t1w_preproc",
        "t2w_preproc",
        "t2w_mask",
        "t1w2t2w_xfm",
    ]
    inputnode = pe.Node(niu.IdentityInterface(fields=inputfields), name="inputnode")

    map_t1w = pe.Node(ApplyTransforms(interpolation="BSpline"), name="map_t1w")
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
        (inputnode, map_t1w, [("t1w_preproc", "input_file"),
                              ("t2w_preproc", "reference_file"),
                              ("t1w2t2w_xfm", "transforms")]),
        (inputnode, norm_rpt, [("t2w_preproc", "before"),
                               ("t2w_mask", "wm_seg")]),
        (inputnode, ds_t1w_t2w_report, [("source_file", "source_file")]),
        (norm_rpt, ds_t1w_t2w_report, [("out_report", "in_file")]),
    ])
    # fmt:on

    return workflow
