# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Writing outputs."""

from __future__ import annotations

import typing as ty
from pathlib import Path

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine import Workflow, tag
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from smriprep.workflows.outputs import init_template_iterator_wf

from nibabies.interfaces import DerivativesDataSink

if ty.TYPE_CHECKING:
    from niworkflows.utils.spaces import SpatialReferences

BIDS_TISSUE_ORDER = ('GM', 'WM', 'CSF')


@tag('anat.coreg-report')
def init_coreg_report_wf(*, output_dir, name='coreg_report_wf'):
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
        'source_file',
        't1w_preproc',
        't2w_preproc',
        'in_mask',
    ]
    inputnode = pe.Node(niu.IdentityInterface(fields=inputfields), name='inputnode')
    # Generate reportlets showing spatial normalization
    norm_rpt = pe.Node(
        SimpleBeforeAfter(before_label='T2w', after_label='T1w'),
        name='norm_rpt',
        mem_gb=0.1,
    )

    ds_t1w_t2w_report = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir, space='T2w', suffix='T1w', datatype='figures'
        ),
        name='ds_t1w_t2w_report',
        run_without_submitting=True,
    )

    # fmt: off
    workflow.connect([
        (inputnode, norm_rpt, [('t2w_preproc', 'before'),
                               ('t1w_preproc', 'after'),
                               ('in_mask', 'wm_seg')]),
        (inputnode, ds_t1w_t2w_report, [('source_file', 'source_file')]),
        (norm_rpt, ds_t1w_t2w_report, [('out_report', 'in_file')]),
    ])
    # fmt: on

    return workflow


@tag('anat.reports')
def init_anat_reports_wf(
    *,
    spaces: SpatialReferences,
    surface_recon: ty.Literal['freesurfer', 'infantfs', 'mcribs'] | None,
    output_dir: str,
    sloppy: bool,
    name='anat_reports_wf',
) -> Workflow:
    """
    Patched workflow for reports to allow no resolution for templates
    Set up a battery of datasinks to store reports in the right location.
    Parameters
    ----------
    recon_method : :obj:`bool`
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
    from niworkflows.interfaces.reportlets.masks import ROIsPlot
    from niworkflows.interfaces.reportlets.registration import (
        SimpleBeforeAfterRPT as SimpleBeforeAfter,
    )
    from smriprep.workflows.outputs import (
        _empty_report,
        _rpt_masks,
    )

    workflow = Workflow(name=name)

    inputfields = [
        'source_file',
        'anat_preproc',
        'anat_dseg',
        'anat_mask',
        'template',
        'anat2std_xfm',
        # Optional
        'subject_id',
        'subjects_dir',
        'anat_conform_report',
    ]
    inputnode = pe.Node(niu.IdentityInterface(fields=inputfields), name='inputnode')

    seg_rpt = pe.Node(ROIsPlot(colors=['b', 'magenta'], levels=[1.5, 2.5]), name='seg_rpt')

    anat_conform_check = pe.Node(
        niu.Function(function=_empty_report),
        name='anat_conform_check',
        run_without_submitting=True,
    )

    ds_anat_conform_report = pe.Node(
        DerivativesDataSink(base_directory=output_dir, desc='conform', datatype='figures'),
        name='ds_anat_conform_report',
        run_without_submitting=True,
    )

    ds_anat_dseg_mask_report = pe.Node(
        DerivativesDataSink(base_directory=output_dir, suffix='dseg', datatype='figures'),
        name='ds_anat_dseg_mask_report',
        run_without_submitting=True,
    )

    workflow.connect([
        (inputnode, anat_conform_check, [('anat_conform_report', 'in_file')]),
        (anat_conform_check, ds_anat_conform_report, [('out', 'in_file')]),
        (inputnode, ds_anat_conform_report, [('source_file', 'source_file')]),
        (inputnode, ds_anat_dseg_mask_report, [('source_file', 'source_file')]),
        (inputnode, seg_rpt, [('anat_preproc', 'in_file'),
                              ('anat_mask', 'in_mask'),
                              ('anat_dseg', 'in_rois')]),
        (seg_rpt, ds_anat_dseg_mask_report, [('out_report', 'in_file')]),
    ])  # fmt:skip

    if spaces._cached is not None and spaces.cached.references:
        template_iterator_wf = init_template_iterator_wf(spaces=spaces, sloppy=sloppy)
        t1w_std = pe.Node(
            ApplyTransforms(
                dimension=3,
                default_value=0,
                float=True,
                interpolation='LanczosWindowedSinc',
            ),
            name='t1w_std',
        )
        mask_std = pe.Node(
            ApplyTransforms(
                dimension=3,
                default_value=0,
                float=True,
                interpolation='MultiLabel',
            ),
            name='mask_std',
        )

        norm_msk = pe.Node(
            niu.Function(
                function=_rpt_masks,
                output_names=['before', 'after'],
                input_names=['mask_file', 'before', 'after', 'after_mask'],
            ),
            name='norm_msk',
        )
        norm_rpt = pe.Node(SimpleBeforeAfter(), name='norm_rpt', mem_gb=0.1)
        norm_rpt.inputs.after_label = 'Participant'  # after

        ds_std_t1w_report = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='T1w', datatype='figures'),
            name='ds_std_t1w_report',
            run_without_submitting=True,
        )

        workflow.connect([
            (inputnode, template_iterator_wf, [
                ('template', 'inputnode.template'),
                ('anat2std_xfm', 'inputnode.anat2std_xfm'),
            ]),
            (inputnode, t1w_std, [('anat_preproc', 'input_image')]),
            (inputnode, mask_std, [('anat_mask', 'input_image')]),
            (template_iterator_wf, t1w_std, [
                ('outputnode.anat2std_xfm', 'transforms'),
                ('outputnode.std_t1w', 'reference_image'),
            ]),
            (template_iterator_wf, mask_std, [
                ('outputnode.anat2std_xfm', 'transforms'),
                ('outputnode.std_t1w', 'reference_image'),
            ]),
            (template_iterator_wf, norm_rpt, [('outputnode.space', 'before_label')]),
            (t1w_std, norm_msk, [('output_image', 'after')]),
            (mask_std, norm_msk, [('output_image', 'after_mask')]),
            (template_iterator_wf, norm_msk, [
                ('outputnode.std_t1w', 'before'),
                ('outputnode.std_mask', 'mask_file'),
            ]),
            (norm_msk, norm_rpt, [
                ('before', 'before'),
                ('after', 'after'),
            ]),
            (inputnode, ds_std_t1w_report, [('source_file', 'source_file')]),
            (template_iterator_wf, ds_std_t1w_report, [('outputnode.space', 'space')]),
            (norm_rpt, ds_std_t1w_report, [('out_report', 'in_file')]),
        ])  # fmt:skip

    if not surface_recon:
        return workflow

    from smriprep.interfaces.reports import FSSurfaceReport

    recon_report = pe.Node(FSSurfaceReport(), name='recon_report')
    recon_report.interface._always_run = True

    if surface_recon == 'freesurfer':
        recon_desc = 'reconall'
    elif surface_recon == 'infantfs':
        recon_desc = 'infantfs'
    elif surface_recon == 'mcribs':
        recon_desc = 'mcribs'

    ds_recon_report = pe.Node(
        DerivativesDataSink(base_directory=output_dir, desc=recon_desc, datatype='figures'),
        name='ds_recon_report',
        run_without_submitting=True,
    )
    workflow.connect([
        (inputnode, recon_report, [('subjects_dir', 'subjects_dir'),
                                   ('subject_id', 'subject_id')]),
        (recon_report, ds_recon_report, [('out_report', 'in_file')]),
        (inputnode, ds_recon_report, [('source_file', 'source_file')])
    ])  # fmt: skip

    return workflow


@tag('anat.derivatives')
def init_anat_derivatives_wf(
    *,
    bids_root: Path | str,
    output_dir: Path | str,
    spaces: SpatialReferences,
    cifti_output: bool,
    num_t1w: int | None,
    num_t2w: int | None,
    surface_recon: ty.Literal['freesurfer', 'infantfs', 'mcribs'] | None,
    tpm_labels: tuple[str, str, str] = BIDS_TISSUE_ORDER,
    name: str = 'anat_derivatives_wf',
):
    """
    Set up a battery of datasinks to store derivatives in the right location.
    Parameters
    ----------
    bids_root : :obj:`str`
        Root path of BIDS dataset
    freesurfer : :obj:`bool`
        FreeSurfer was enabled
    num_t1w : :obj:`int`
        Number of T1w images
    output_dir : :obj:`str`
        Directory in which to save derivatives
    name : :obj:`str`
        Workflow name (default: anat_derivatives_wf)
    tpm_labels : :obj:`tuple`
        Tissue probability maps in order
    Inputs
    ------
    template
        Template space and specifications
    source_files
        List of input T1w images
    t1w_ref_xfms
        List of affine transforms to realign input T1w images
    t1w_preproc
        The T1w reference map, which is calculated as the average of bias-corrected
        and preprocessed T1w images, defining the anatomical space.
    t1w_mask
        Mask of the ``t1w_preproc``
    t1w_dseg
        Segmentation in T1w space
    t1w_tpms
        Tissue probability maps in T1w space
    anat2std_xfm
        Nonlinear spatial transform to resample imaging data given in anatomical space
        into standard space.
    std2anat_xfm
        Inverse transform of ``anat2std_xfm``
    std_t1w
        T1w reference resampled in one or more standard spaces.
    std_mask
        Mask of skull-stripped template, in standard space
    std_dseg
        Segmentation, resampled into standard space
    std_tpms
        Tissue probability maps in standard space
    t1w2fsnative_xfm
        LTA-style affine matrix translating from T1w to
        FreeSurfer-conformed subject space
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed
        subject space to T1w
    surfaces
        GIFTI surfaces (gray/white boundary, midthickness, pial, inflated)
    t1w_fs_aseg
        FreeSurfer's aseg segmentation, in native T1w space
    t1w_fs_aparc
        FreeSurfer's aparc+aseg segmentation, in native T1w space
    t2w_source_files
        List of input T2w images
    t2w_preproc
        The T2w image in T1w space.
    cifti_morph
        Morphometric CIFTI-2 dscalar files
    cifti_density
        Grayordinate density
    cifti_metadata
        JSON files containing metadata dictionaries
    """
    from niworkflows.interfaces.nibabel import ApplyMask
    from niworkflows.interfaces.utility import KeySelect
    from smriprep.workflows.outputs import (
        _bids_relative,
        _combine_cohort,
        _drop_path,
        _fmt_cohort,
        _is_native,
    )

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'template',
                # T1w
                't1w_source_files',
                't1w_ref_xfms',
                't1w_preproc',
                # T2w
                't2w_source_files',
                't2w_ref_xfms',
                't2w_preproc',
                # Can be in either T1w/T2w space
                'anat_mask',
                'anat_dseg',
                'anat_tpms',
                'anat2std_xfm',
                'std2anat_xfm',
                # FS
                'anat2fsnative_xfm',
                'fsnative2anat_xfm',
                'anat_fs_aseg',
                'anat_fs_aparc',
                'anat_ribbon',
                'surfaces',
                'morphometrics',
                # CIFTI
                'cifti_metadata',
                'cifti_density',
                'cifti_morph',
                'sphere_reg',
                'sphere_reg_fsLR',
            ]
        ),
        name='inputnode',
    )
    # The preferred space to use for to/from entities
    source_files = 't1w_source_files' if num_t1w else 't2w_source_files'
    space = 'T1w' if num_t1w else 'T2w'

    if num_t1w:
        raw_sources = pe.Node(niu.Function(function=_bids_relative), name='t1w_raw_sources')
        raw_sources.inputs.bids_root = bids_root

        ds_t1w_preproc = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='preproc', compress=True),
            name='ds_t1w_preproc',
            run_without_submitting=True,
        )
        ds_t1w_preproc.inputs.SkullStripped = False

        if num_t1w > 1:
            # Please note the dictionary unpacking to provide the from argument.
            # It is necessary because from is a protected keyword (not allowed as argument name).
            ds_t1w_ref_xfms = pe.MapNode(
                DerivativesDataSink(
                    base_directory=output_dir,
                    to='T1w',
                    mode='image',
                    suffix='xfm',
                    extension='txt',
                    **{'from': 'orig'},
                ),
                iterfield=['source_file', 'in_file'],
                name='ds_t1w_ref_xfms',
                run_without_submitting=True,
            )
            # fmt:off
            workflow.connect([
                (inputnode, ds_t1w_ref_xfms, [('t1w_source_files', 'source_file'),
                                              ('t1w_ref_xfms', 'in_file')]),
            ])
            # fmt:on

    if num_t2w:
        if not num_t1w:
            raw_sources = pe.Node(niu.Function(function=_bids_relative), name='t2w_raw_sources')
            raw_sources.inputs.bids_root = bids_root

        ds_t2w_preproc = pe.Node(
            DerivativesDataSink(base_directory=output_dir, desc='preproc', compress=True),
            name='ds_t2w_preproc',
            run_without_submitting=True,
        )
        ds_t2w_preproc.inputs.SkullStripped = False

        if num_t1w:
            ds_t2w_preproc.inputs.space = 'T1w'

        if num_t2w > 1:
            # Please note the dictionary unpacking to provide the from argument.
            # It is necessary because from is a protected keyword (not allowed as argument name).
            ds_t2w_ref_xfms = pe.MapNode(
                DerivativesDataSink(
                    base_directory=output_dir,
                    to='T1w',
                    mode='image',
                    suffix='xfm',
                    extension='txt',
                    **{'from': 'orig'},
                ),
                iterfield=['source_file', 'in_file'],
                name='ds_t2w_ref_xfms',
                run_without_submitting=True,
            )
            # fmt:off
            workflow.connect([
                (inputnode, ds_t2w_ref_xfms, [('t2w_source_files', 'source_file'),
                                              ('t2w_ref_xfms', 'in_file')]),
            ])
            # fmt:on

    ds_anat_mask = pe.Node(
        DerivativesDataSink(base_directory=output_dir, desc='brain', suffix='mask', compress=True),
        name='ds_anat_mask',
        run_without_submitting=True,
    )
    ds_anat_mask.inputs.Type = 'Brain'

    ds_anat_dseg = pe.Node(
        DerivativesDataSink(base_directory=output_dir, suffix='dseg', compress=True),
        name='ds_anat_dseg',
        run_without_submitting=True,
    )

    ds_anat_tpms = pe.Node(
        DerivativesDataSink(base_directory=output_dir, suffix='probseg', compress=True),
        name='ds_anat_tpms',
        run_without_submitting=True,
    )
    ds_anat_tpms.inputs.label = tpm_labels

    ds_anat_ribbon = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc='ribbon',
            suffix='mask',
            extension='.nii.gz',
            compress=True,
        ),
        name='ds_anat_ribbon',
        run_without_submitting=True,
    )

    if num_t1w:
        workflow.connect([
            (inputnode, ds_t1w_preproc, [('t1w_preproc', 'in_file'),
                                         ('t1w_source_files', 'source_file')]),
        ])  # fmt:skip

    if num_t2w:
        workflow.connect([
            (inputnode, ds_t2w_preproc, [('t2w_preproc', 'in_file'),
                                         ('t2w_source_files', 'source_file')]),
        ])  # fmt:skip

    # fmt:off
    workflow.connect([
        (inputnode, raw_sources, [(source_files, 'in_files')]),
        (inputnode, ds_anat_mask, [('anat_mask', 'in_file'),
                                   (source_files, 'source_file')]),
        (inputnode, ds_anat_tpms, [('anat_tpms', 'in_file'),
                                   (source_files, 'source_file')]),
        (inputnode, ds_anat_dseg, [('anat_dseg', 'in_file'),
                                   (source_files, 'source_file')]),
        (inputnode, ds_anat_ribbon, [('anat_ribbon', 'in_file'),
                                     (source_files, 'source_file')]),
        (raw_sources, ds_anat_mask, [('out', 'RawSources')]),
    ])
    # fmt:on

    # Transforms
    if spaces.get_spaces(nonstandard=False, dim=(3,)):
        ds_std2anat_xfm = pe.MapNode(
            DerivativesDataSink(base_directory=output_dir, to=space, mode='image', suffix='xfm'),
            iterfield=('in_file', 'from'),
            name='ds_std2anat_xfm',
            run_without_submitting=True,
        )

        ds_anat2std_xfm = pe.MapNode(
            DerivativesDataSink(
                base_directory=output_dir, mode='image', suffix='xfm', **{'from': space}
            ),
            iterfield=('in_file', 'to'),
            name='ds_anat2std_xfm',
            run_without_submitting=True,
        )

        # fmt:off
        workflow.connect([
            (inputnode, ds_anat2std_xfm, [
                ('anat2std_xfm', 'in_file'),
                (('template', _combine_cohort), 'to'),
                (source_files, 'source_file')]),
            (inputnode, ds_std2anat_xfm, [
                ('std2anat_xfm', 'in_file'),
                (('template', _combine_cohort), 'from'),
                (source_files, 'source_file')]),
        ])
        # fmt:on

    # Write derivatives in standard spaces specified by --output-spaces
    if spaces._cached is not None and spaces.cached.references:
        from niworkflows.interfaces.fixes import (
            FixHeaderApplyTransforms as ApplyTransforms,
        )
        from niworkflows.interfaces.nibabel import GenerateSamplingReference
        from niworkflows.interfaces.space import SpaceDataSource
        from smriprep.interfaces.templateflow import TemplateFlowSelect

        spacesource = pe.Node(SpaceDataSource(), name='spacesource', run_without_submitting=True)
        spacesource.iterables = (
            'in_tuple',
            [(s.fullname, s.spec) for s in spaces.cached.get_standard(dim=(3,))],
        )

        gen_tplid = pe.Node(
            niu.Function(function=_fmt_cohort),
            name='gen_tplid',
            run_without_submitting=True,
        )

        select_xfm = pe.Node(
            KeySelect(fields=['anat2std_xfm']),
            name='select_xfm',
            run_without_submitting=True,
        )
        select_tpl = pe.Node(TemplateFlowSelect(), name='select_tpl', run_without_submitting=True)

        gen_ref = pe.Node(GenerateSamplingReference(), name='gen_ref', mem_gb=0.01)

        mask_anat = pe.Node(ApplyMask(), name='mask_anat')

        # Resample T1w-space inputs
        anat2std_t1w = pe.Node(
            ApplyTransforms(
                dimension=3,
                default_value=0,
                float=True,
                interpolation='LanczosWindowedSinc',
            ),
            name='anat2std_t1w',
        )

        anat2std_mask = pe.Node(ApplyTransforms(interpolation='MultiLabel'), name='anat2std_mask')
        anat2std_dseg = pe.Node(ApplyTransforms(interpolation='MultiLabel'), name='anat2std_dseg')
        anat2std_tpms = pe.MapNode(
            ApplyTransforms(dimension=3, default_value=0, float=True, interpolation='Gaussian'),
            iterfield=['input_image'],
            name='anat2std_tpms',
        )

        ds_std_t1w = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir,
                desc='preproc',
                compress=True,
            ),
            name='ds_std_t1w',
            run_without_submitting=True,
        )
        ds_std_t1w.inputs.SkullStripped = True

        ds_std_mask = pe.Node(
            DerivativesDataSink(
                base_directory=output_dir, desc='brain', suffix='mask', compress=True
            ),
            name='ds_std_mask',
            run_without_submitting=True,
        )
        ds_std_mask.inputs.Type = 'Brain'

        ds_std_dseg = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='dseg', compress=True),
            name='ds_std_dseg',
            run_without_submitting=True,
        )

        ds_std_tpms = pe.Node(
            DerivativesDataSink(base_directory=output_dir, suffix='probseg', compress=True),
            name='ds_std_tpms',
            run_without_submitting=True,
        )

        set_tpl_res = pe.Node(
            niu.Function(function=_set_tpl_res),
            name='set_tpl_res',
            run_without_submitting=True,
            mem_gb=0.1,
        )

        # CRITICAL: the sequence of labels here (CSF-GM-WM) is that of the output of FSL-FAST
        #           (intensity mean, per tissue). This order HAS to be matched also by the ``tpms``
        #           output in the data/io_spec.json file.
        ds_std_tpms.inputs.label = tpm_labels

        preproc_file = 't1w_preproc' if num_t1w else 't2w_preproc'

        # fmt: off
        workflow.connect([
            (inputnode, mask_anat, [(preproc_file, 'in_file'),
                                    ('anat_mask', 'in_mask')]),
            (mask_anat, anat2std_t1w, [('out_file', 'input_image')]),
            (inputnode, anat2std_mask, [('anat_mask', 'input_image')]),
            (inputnode, anat2std_dseg, [('anat_dseg', 'input_image')]),
            (inputnode, anat2std_tpms, [('anat_tpms', 'input_image')]),
            (inputnode, gen_ref, [(preproc_file, 'moving_image')]),
            (inputnode, select_xfm, [
                ('anat2std_xfm', 'anat2std_xfm'),
                ('template', 'keys')]),
            (spacesource, gen_tplid, [('space', 'template'),
                                      ('cohort', 'cohort')]),
            (gen_tplid, select_xfm, [('out', 'key')]),
            (spacesource, select_tpl, [('space', 'template'),
                                       ('cohort', 'cohort')]),
            (spacesource, set_tpl_res, [('space', 'space'),
                                        ('resolution', 'resolution')]),
            (set_tpl_res, select_tpl, [('out', 'resolution')]),
            (spacesource, gen_ref, [(('resolution', _is_native), 'keep_native')]),
            (select_tpl, gen_ref, [('t1w_file', 'fixed_image')]),
            (anat2std_t1w, ds_std_t1w, [('output_image', 'in_file')]),
            (anat2std_mask, ds_std_mask, [('output_image', 'in_file')]),
            (anat2std_dseg, ds_std_dseg, [('output_image', 'in_file')]),
            (anat2std_tpms, ds_std_tpms, [('output_image', 'in_file')]),
            (select_tpl, ds_std_mask, [(('brain_mask', _drop_path), 'RawSources')]),
        ])
        # fmt: on
        workflow.connect(
            # Connect apply transforms nodes
            [
                (gen_ref, n, [('out_file', 'reference_image')])
                for n in (anat2std_t1w, anat2std_mask, anat2std_dseg, anat2std_tpms)
            ]
            + [
                (select_xfm, n, [('anat2std_xfm', 'transforms')])
                for n in (anat2std_t1w, anat2std_mask, anat2std_dseg, anat2std_tpms)
            ]
            # Connect the source_file input of these datasinks
            + [
                (inputnode, n, [(source_files, 'source_file')])
                for n in (ds_std_t1w, ds_std_mask, ds_std_dseg, ds_std_tpms)
            ]
            # Connect the space input of these datasinks
            + [
                (
                    spacesource,
                    n,
                    [('space', 'space'), ('cohort', 'cohort'), ('resolution', 'resolution')],
                )
                for n in (ds_std_t1w, ds_std_mask, ds_std_dseg, ds_std_tpms)
            ]
        )

    if not surface_recon:
        return workflow

    from niworkflows.interfaces.nitransforms import ConcatenateXFMs
    from niworkflows.interfaces.surf import Path2BIDS

    # FS native space transforms
    lta2itk_fwd = pe.Node(ConcatenateXFMs(), name='lta2itk_fwd', run_without_submitting=True)
    lta2itk_inv = pe.Node(ConcatenateXFMs(), name='lta2itk_inv', run_without_submitting=True)
    ds_anat_fsnative = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            mode='image',
            to='fsnative',
            suffix='xfm',
            extension='txt',
            **{'from': space},
        ),
        name='ds_anat_fsnative',
        run_without_submitting=True,
    )
    ds_fsnative_anat = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            mode='image',
            to=space,
            suffix='xfm',
            extension='txt',
            **{'from': 'fsnative'},
        ),
        name='ds_fsnative_anat',
        run_without_submitting=True,
    )
    # Surfaces
    name_surfs = pe.MapNode(
        Path2BIDS(), iterfield='in_file', name='name_surfs', run_without_submitting=True
    )
    ds_surfs = pe.MapNode(
        DerivativesDataSink(base_directory=output_dir, extension='.surf.gii'),
        iterfield=['in_file', 'hemi', 'suffix'],
        name='ds_surfs',
        run_without_submitting=True,
    )
    name_regs = pe.MapNode(
        Path2BIDS(), iterfield='in_file', name='name_regs', run_without_submitting=True
    )
    ds_regs = pe.MapNode(
        DerivativesDataSink(
            base_directory=output_dir,
            space='fsaverage',
            desc='reg',
            suffix='sphere',
            extension='.surf.gii',
        ),
        iterfield=['in_file', 'hemi'],
        name='ds_regs',
        run_without_submitting=True,
    )
    name_reg_fsLR = pe.MapNode(
        Path2BIDS(), iterfield='in_file', name='name_reg_fsLR', run_without_submitting=True
    )
    ds_reg_fsLR = pe.MapNode(
        DerivativesDataSink(
            base_directory=output_dir,
            space='dHCP' if surface_recon == 'mcribs' else 'fsLR',
            desc='reg',
            suffix='sphere',
            extension='.surf.gii',
        ),
        iterfield=['in_file', 'hemi'],
        name='ds_reg_fsLR',
        run_without_submitting=True,
    )
    # Morphometrics
    name_morphs = pe.MapNode(
        Path2BIDS(),
        iterfield='in_file',
        name='name_morphs',
        run_without_submitting=True,
    )
    ds_morphs = pe.MapNode(
        DerivativesDataSink(base_directory=output_dir, extension='.shape.gii'),
        iterfield=['in_file', 'hemi', 'suffix'],
        name='ds_morphs',
        run_without_submitting=True,
    )
    # Parcellations
    ds_anat_fsaseg = pe.Node(
        DerivativesDataSink(base_directory=output_dir, desc='aseg', suffix='dseg', compress=True),
        name='ds_anat_fsaseg',
        run_without_submitting=True,
    )
    ds_anat_fsparc = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir, desc='aparcaseg', suffix='dseg', compress=True
        ),
        name='ds_anat_fsparc',
        run_without_submitting=True,
    )

    # fmt: off
    workflow.connect([
        (inputnode, lta2itk_fwd, [('anat2fsnative_xfm', 'in_xfms')]),
        (inputnode, lta2itk_inv, [('fsnative2anat_xfm', 'in_xfms')]),
        (inputnode, ds_anat_fsnative, [(source_files, 'source_file')]),
        (lta2itk_fwd, ds_anat_fsnative, [('out_xfm', 'in_file')]),
        (inputnode, ds_fsnative_anat, [(source_files, 'source_file')]),
        (lta2itk_inv, ds_fsnative_anat, [('out_xfm', 'in_file')]),
        (inputnode, name_surfs, [('surfaces', 'in_file')]),
        (inputnode, ds_surfs, [('surfaces', 'in_file'),
                               (source_files, 'source_file')]),
        (name_surfs, ds_surfs, [('hemi', 'hemi'),
                                ('suffix', 'suffix')]),
        (inputnode, name_regs, [('sphere_reg', 'in_file')]),
        (inputnode, ds_regs, [('sphere_reg', 'in_file'),
                              (source_files, 'source_file')]),
        (name_regs, ds_regs, [('hemi', 'hemi')]),
        (inputnode, name_reg_fsLR, [('sphere_reg_fsLR', 'in_file')]),
        (inputnode, ds_reg_fsLR, [('sphere_reg_fsLR', 'in_file'),
                                  (source_files, 'source_file')]),
        (name_reg_fsLR, ds_reg_fsLR, [('hemi', 'hemi')]),
        (inputnode, name_morphs, [('morphometrics', 'in_file')]),
        (inputnode, ds_morphs, [('morphometrics', 'in_file'),
                                (source_files, 'source_file')]),
        (name_morphs, ds_morphs, [('hemi', 'hemi'),
                                  ('suffix', 'suffix')]),
        (inputnode, ds_anat_fsaseg, [('anat_fs_aseg', 'in_file'),
                                     (source_files, 'source_file')]),
        (inputnode, ds_anat_fsparc, [('anat_fs_aparc', 'in_file'),
                                     (source_files, 'source_file')]),
    ])
    # fmt: on
    if cifti_output:
        ds_cifti_morph = pe.MapNode(
            DerivativesDataSink(
                base_directory=output_dir,
                suffix=['curv', 'sulc', 'thickness'],
                compress=False,
                space='fsLR',
            ),
            name='ds_cifti_morph',
            run_without_submitting=True,
            iterfield=['in_file', 'meta_dict', 'suffix'],
        )
        # fmt:off
        workflow.connect([
            (inputnode, ds_cifti_morph, [('cifti_morph', 'in_file'),
                                         (source_files, 'source_file'),
                                         ('cifti_density', 'density'),
                                         (('cifti_metadata', _read_jsons), 'meta_dict')])
        ])
        # fmt:on
    return workflow


def init_ds_seg_wf(
    *,
    output_dir: str,
    seg_type: str,
    extra_entities: dict | None = None,
):
    """
    Set up a battery of datasinks to store derivatives in the right location.

    Parameters
    ----------
    bids_root : :obj:`str`
        Root path of BIDS dataset
    output_dir : :obj:`str`
        Directory in which to save derivatives
    seg_type : :obj:`str`
        Type of segmentation (aseg, aparcaseg, etc)
    extra_entities : :obj:`dict` or None
        Additional entities to add to filename
    name : :obj:`str`
        Workflow name (default: ds_anat_segs_wf)

    Inputs
    ------
    in_seg
        Input segmentation, in native anatomical space
    source_files
        List of input anatomical images
    """
    workflow = Workflow(name=f'ds_{seg_type}_wf')

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['source_files', 'in_seg']),
        name='inputnode',
    )

    extra_entities = extra_entities or {}

    ds_seg = pe.Node(
        DerivativesDataSink(
            base_directory=output_dir,
            desc=seg_type,
            suffix='dseg',
            compress=True,
            **extra_entities,
        ),
        name='ds_anat_fsaseg',
        run_without_submitting=True,
    )

    workflow.connect([
        (inputnode, ds_seg, [
            ('in_seg', 'in_file'),
            ('source_files', 'source_file'),
        ]),
    ])  # fmt:skip

    return workflow


def _set_tpl_res(space, resolution):
    if space in ('UNCInfant', 'Fischer344'):
        from nipype.interfaces.base import Undefined

        return Undefined
    try:
        return int(resolution)
    except ValueError:
        return 1


def _read_jsons(in_file):
    from json import loads
    from pathlib import Path

    return [loads(Path(f).read_text()) for f in in_file]
