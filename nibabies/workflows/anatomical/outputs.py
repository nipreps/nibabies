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


def init_anat_reports_wf(*, freesurfer, output_dir, sloppy, name="anat_reports_wf"):
    """
    Patched workflow for reports to allow no resolution for templates
    Set up a battery of datasinks to store reports in the right location.
    Parameters
    ----------
    freesurfer : :obj:`bool`
        FreeSurfer was enabled
    output_dir : :obj:`str`
        Directory in which to save derivatives
    name : :obj:`str`
        Workflow name (default: anat_reports_wf)
    Inputs
    ------
    source_file
        Input T1w image
    std_t1w
        T1w image resampled to standard space
    std_mask
        Mask of skull-stripped template
    subject_dir
        FreeSurfer SUBJECTS_DIR
    subject_id
        FreeSurfer subject ID
    t1w_conform_report
        Conformation report
    t1w_preproc
        The T1w reference map, which is calculated as the average of bias-corrected
        and preprocessed T1w images, defining the anatomical space.
    t1w_dseg
        Segmentation in T1w space
    t1w_mask
        Brain (binary) mask estimated by brain extraction.
    template
        Template space and specifications
    """
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from niworkflows.interfaces.reportlets.masks import ROIsPlot
    from smriprep.interfaces.templateflow import TemplateFlowSelect
    from smriprep.workflows.outputs import (
        _fmt, _empty_report, _rpt_masks, _drop_cohort, _pick_cohort,
    )
    from ...utils.patches import set_tf_resolution

    workflow = Workflow(name=name)

    inputfields = [
        "source_file",
        "t1w_conform_report",
        "t1w_preproc",
        "t1w_dseg",
        "t1w_mask",
        "template",
        "std_t1w",
        "std_mask",
        "subject_id",
        "subjects_dir",
    ]
    inputnode = pe.Node(niu.IdentityInterface(fields=inputfields), name="inputnode")

    seg_rpt = pe.Node(
        ROIsPlot(colors=["b", "magenta"], levels=[1.5, 2.5]), name="seg_rpt"
    )

    t1w_conform_check = pe.Node(
        niu.Function(function=_empty_report),
        name="t1w_conform_check",
        run_without_submitting=True,
    )

    ds_t1w_conform_report = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir, desc="conform", datatype="figures"
        ),
        name="ds_t1w_conform_report",
        run_without_submitting=True,
    )

    ds_t1w_dseg_mask_report = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir, suffix="dseg", datatype="figures"
        ),
        name="ds_t1w_dseg_mask_report",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, t1w_conform_check, [('t1w_conform_report', 'in_file')]),
        (t1w_conform_check, ds_t1w_conform_report, [('out', 'in_file')]),
        (inputnode, ds_t1w_conform_report, [('source_file', 'source_file')]),
        (inputnode, ds_t1w_dseg_mask_report, [('source_file', 'source_file')]),
        (inputnode, seg_rpt, [('t1w_preproc', 'in_file'),
                              ('t1w_mask', 'in_mask'),
                              ('t1w_dseg', 'in_rois')]),
        (seg_rpt, ds_t1w_dseg_mask_report, [('out_report', 'in_file')]),
    ])
    # fmt:on

    # Generate reportlets showing spatial normalization
    tf_select = pe.Node(
        TemplateFlowSelect(), name="tf_select", run_without_submitting=True
    )

    set_tf_res = pe.Node(niu.Function(function=set_tf_resolution), name='set_tf_res')
    set_tf_res.inputs.sloppy = sloppy

    norm_msk = pe.Node(
        niu.Function(
            function=_rpt_masks,
            output_names=["before", "after"],
            input_names=["mask_file", "before", "after", "after_mask"],
        ),
        name="norm_msk",
    )
    norm_rpt = pe.Node(SimpleBeforeAfter(), name="norm_rpt", mem_gb=0.1)
    norm_rpt.inputs.after_label = "Participant"  # after

    ds_std_t1w_report = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir, suffix="T1w", datatype="figures"
        ),
        name="ds_std_t1w_report",
        run_without_submitting=True,
    )

    # fmt:off
    workflow.connect([
        (inputnode, set_tf_res, [(('template', _drop_cohort), 'template')]),
        (set_tf_res, tf_select, [('out', 'resolution')]),
        (inputnode, tf_select, [(('template', _drop_cohort), 'template'),
                                (('template', _pick_cohort), 'cohort')]),
        (inputnode, norm_rpt, [('template', 'before_label')]),
        (inputnode, norm_msk, [('std_t1w', 'after'),
                               ('std_mask', 'after_mask')]),
        (tf_select, norm_msk, [('t1w_file', 'before'),
                               ('brain_mask', 'mask_file')]),
        (norm_msk, norm_rpt, [('before', 'before'),
                              ('after', 'after')]),
        (inputnode, ds_std_t1w_report, [
            (('template', _fmt), 'space'),
            ('source_file', 'source_file')]),
        (norm_rpt, ds_std_t1w_report, [('out_report', 'in_file')]),
    ])
    # fmt:on

    if freesurfer:
        from smriprep.interfaces.reports import FSSurfaceReport

        recon_report = pe.Node(FSSurfaceReport(), name="recon_report")
        recon_report.interface._always_run = True

        ds_recon_report = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir, desc="reconall", datatype="figures"
            ),
            name="ds_recon_report",
            run_without_submitting=True,
        )
        # fmt:off
        workflow.connect([
            (inputnode, recon_report, [('subjects_dir', 'subjects_dir'),
                                       ('subject_id', 'subject_id')]),
            (recon_report, ds_recon_report, [('out_report', 'in_file')]),
            (inputnode, ds_recon_report, [('source_file', 'source_file')])
        ])
        # fmt:on

    return workflow
