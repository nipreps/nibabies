from pathlib import Path

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl


def init_infant_anat_wf(
    *,
    age_months,
    t1w,
    t2w,
    anat_modality,
    bids_root,
    existing_derivatives,
    freesurfer,
    longitudinal,
    omp_nthreads,
    output_dir,
    segmentation_atlases,
    skull_strip_mode,
    skull_strip_template,
    sloppy,
    spaces,
    name="infant_anat_wf",
):
    """

      - T1w reference: realigning and then averaging anatomical images.
      - Brain extraction and INU (bias field) correction.
      - Brain tissue segmentation.
      - Spatial normalization to standard spaces.
      - Surface reconstruction with FreeSurfer_.

    Outputs
    -------

    anat_preproc
        The anatomical reference map, which is calculated as the average of bias-corrected
        and preprocessed anatomical images, defining the anatomical space.
    anat_brain
        Skull-stripped ``anat_preproc``
    anat_mask
        Brain (binary) mask estimated by brain extraction.
    anat_dseg
        Brain tissue segmentation of the preprocessed structural image, including
        gray-matter (GM), white-matter (WM) and cerebrospinal fluid (CSF).
    anat_tpms
        List of tissue probability maps corresponding to ``t1w_dseg``.
    std_preproc
        T1w reference resampled in one or more standard spaces.
    std_mask
        Mask of skull-stripped template, in MNI space
    std_dseg
        Segmentation, resampled into MNI space
    std_tpms
        List of tissue probability maps in MNI space
    subjects_dir
        FreeSurfer SUBJECTS_DIR
    anat2std_xfm
        Nonlinear spatial transform to resample imaging data given in anatomical space
        into standard space.
    std2anat_xfm
        Inverse transform of the above.
    subject_id
        FreeSurfer subject ID
    anat2fsnative_xfm
        LTA-style affine matrix translating from T1w to
        FreeSurfer-conformed subject space
    fsnative2anat_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed
        subject space to T1w
    surfaces
        GIFTI surfaces (gray/white boundary, midthickness, pial, inflated)
    """
    from nipype.interfaces.base import Undefined
    from nipype.interfaces.ants.base import Info as ANTsInfo
    from niworkflows.interfaces.images import ValidateImage
    from smriprep.workflows.anatomical import init_anat_template_wf, _probseg_fast2bids, _pop
    from smriprep.workflows.norm import init_anat_norm_wf
    from smriprep.workflows.outputs import (
        init_anat_reports_wf,
        init_anat_derivatives_wf,
    )

    from ...utils.misc import fix_multi_source_name
    from .brain_extraction import init_infant_brain_extraction_wf
    from .segmentation import init_anat_seg_wf
    from .surfaces import init_infant_surface_recon_wf

    # for now, T1w only
    num_t1w = len(t1w) if t1w else 0
    num_t2w = len(t2w) if t2w else 0

    wf = pe.Workflow(name=name)
    desc = """Anatomical data preprocessing

: """
    desc += f"""\
A total of {num_t1w} T1w and {num_t2w} T2w images were found within the input
BIDS dataset."""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=["t1w", "t2w", "subject_id", "subjects_dir"]
        ),  # FLAIR / ROI?
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "anat_preproc",
                "anat_brain",
                "anat_mask",
                "anat_dseg",
                "anat_tpms",
                "anat_ref_xfms",
                "std_preproc",
                "std_brain",
                "std_dseg",
                "std_tpms",
                "subjects_dir",
                "subject_id",
                "anat2std_xfm",
                "std2anat_xfm",
                "anat2fsnative_xfm",
                "fsnative2anat_xfm",
                "surfaces",
                "anat_aseg",
                "anat_aparc",
            ]
        ),
        name="outputnode",
    )

    # Connect reportlets workflows
    anat_reports_wf = init_anat_reports_wf(
        freesurfer=freesurfer,
        output_dir=output_dir,
    )

    if existing_derivatives:
        raise NotImplementedError("Reusing derivatives is not yet supported.")

    desc += """
All of the T1-weighted images were corrected for intensity non-uniformity (INU)
""" if num_t1w > 1 else """\
The T1-weighted (T1w) image was corrected for intensity non-uniformity (INU)
"""
    desc += """\
with `N4BiasFieldCorrection` [@n4], distributed with ANTs {ants_ver} \
[@ants, RRID:SCR_004757]"""
    desc += '.\n' if num_t1w > 1 else ", and used as T1w-reference throughout the workflow.\n"

    desc += """\
The T1w-reference was then skull-stripped with a modified implementation of
the `antsBrainExtraction.sh` workflow (from ANTs), using {skullstrip_tpl}
as target template.
Brain tissue segmentation of cerebrospinal fluid (CSF),
white-matter (WM) and gray-matter (GM) was performed on
the brain-extracted T1w using ANTs JointFusion, distributed with ANTs {ants_ver}.
"""

    wf.__desc__ = desc.format(
        ants_ver=ANTsInfo.version() or '(version unknown)',
        skullstrip_tpl=skull_strip_template.fullname,
    )
    # Define output workflows
    anat_reports_wf = init_anat_reports_wf(freesurfer=freesurfer, output_dir=output_dir)
    # HACK: remove resolution from TFSelect
    anat_reports_wf.get_node('tf_select').inputs.resolution = Undefined

    anat_derivatives_wf = init_anat_derivatives_wf(
        bids_root=bids_root,
        freesurfer=freesurfer,
        num_t1w=num_t1w,
        output_dir=output_dir,
        spaces=spaces,
    )
    # HACK: remove resolution from TFSelect
    anat_derivatives_wf.get_node('select_tpl').inputs.resolution = Undefined

    # Multiple T1w files -> generate average reference
    t1w_template_wf = init_anat_template_wf(
        longitudinal=False,
        omp_nthreads=omp_nthreads,
        num_t1w=num_t1w,
    )

    use_t2w = False
    if num_t2w:
        t2w_template_wf = init_t2w_template_wf(
            longitudinal=longitudinal,
            omp_nthreads=omp_nthreads,
            num_t2w=num_t2w,
        )
        wf.connect(inputnode, 't2w', t2w_template_wf, 'inputnode.t2w')
        # TODO: determine cutoff (< 8 months)
        use_t2w = True

    anat_validate = pe.Node(
        ValidateImage(),
        name='anat_validate',
        run_without_submitting=True,
    )

    # INU + Brain Extraction
    if skull_strip_mode != 'force':
        raise NotImplementedError("Skull stripping is currently required.")

    brain_extraction_wf = init_infant_brain_extraction_wf(
        age_months=age_months,
        mri_scheme=anat_modality.capitalize(),
        ants_affine_init=True,
        skull_strip_template=skull_strip_template.space,
        template_specs=skull_strip_template.spec,
        omp_nthreads=omp_nthreads,
        output_dir=Path(output_dir),
        sloppy=sloppy,
        use_t2w=use_t2w,
    )
    # Ensure single outputs
    be_buffer = pe.Node(
        niu.IdentityInterface(
            fields=["anat_preproc", "anat_brain"]
        ),
        name='be_buffer'
    )

    # Segmentation - initial implementation should be simple: JLF
    anat_seg_wf = init_anat_seg_wf(
        age_months=age_months,
        anat_modality=anat_modality.capitalize(),
        template_dir=segmentation_atlases,
        sloppy=sloppy,
        omp_nthreads=omp_nthreads,
    )

    # Spatial normalization (requires segmentation)
    anat_norm_wf = init_anat_norm_wf(
        debug=sloppy,
        omp_nthreads=omp_nthreads,
        templates=spaces.get_spaces(nonstandard=False, dim=(3,)),
    )
    # HACK: remove resolution from TFSelect
    anat_norm_wf.get_node('tf_select').inputs.resolution = Undefined
    # HACK: requires patched niworkflows to allow setting resolution to none
    anat_norm_wf.get_node('registration').inputs.template_resolution = None

    # fmt: off
    if use_t2w:
        wf.connect(t2w_template_wf, 'outputnode.t2w_ref', brain_extraction_wf, 'inputnode.t2w')

    wf.connect([
        (inputnode, t1w_template_wf, [
            ('t1w', 'inputnode.t1w'),
        ]),
        (t1w_template_wf, outputnode, [
            ('outputnode.t1w_realign_xfm', 'anat_ref_xfms'),
        ]),
        (t1w_template_wf, anat_validate, [
            ('outputnode.t1w_ref', 'in_file'),
        ]),
        (anat_validate, brain_extraction_wf, [
            ('out_file', 'inputnode.t1w'),
        ]),
        (brain_extraction_wf, be_buffer, [
            (('outputnode.t1w_corrected', _pop), 'anat_preproc'),
            (('outputnode.t1w_corrected_brain', _pop), 'anat_brain'),
            (('outputnode.t1w_mask', _pop), 'anat_mask'),
        ]),
        (be_buffer, outputnode, [
            ('anat_preproc', 'anat_preproc'),
            ('anat_brain', 'anat_brain'),
            ('anat_mask', 'anat_mask'),
        ]),
        (be_buffer, anat_seg_wf, [
            ('anat_brain', 'inputnode.anat_brain'),
        ]),
        (anat_seg_wf, outputnode, [
            ('outputnode.anat_dseg', 'anat_dseg'),
        ]),
        (anat_seg_wf, anat_norm_wf, [
            ('outputnode.anat_dseg', 'inputnode.moving_segmentation'),
            ('outputnode.anat_tpms', 'inputnode.moving_tpms'),
        ]),
        (be_buffer, anat_norm_wf, [
            ('anat_preproc', 'inputnode.moving_image'),
            ('anat_mask', 'inputnode.moving_mask'),
        ]),
        (anat_norm_wf, outputnode, [
            ('poutputnode.standardized', 'std_preproc'),
            ('poutputnode.std_mask', 'std_mask'),
            ('poutputnode.std_dseg', 'std_dseg'),
            ('poutputnode.std_tpms', 'std_tpms'),
            ('outputnode.template', 'template'),
            ('outputnode.anat2std_xfm', 'anat2std_xfm'),
            ('outputnode.std2anat_xfm', 'std2anat_xfm'),
        ]),
        (inputnode, anat_norm_wf, [
            (('t1w', fix_multi_source_name), 'inputnode.orig_t1w'),  # anat_validate? not used...
        ]),
    ])

    wf.connect([
        # reports
        (inputnode, anat_reports_wf, [
            ('t1w', 'inputnode.source_file'),
        ]),
        (outputnode, anat_reports_wf, [
            ('anat_preproc', 'inputnode.t1w_preproc'),
            ('anat_mask', 'inputnode.t1w_mask'),
            ('anat_dseg', 'inputnode.t1w_dseg'),
            ('std_preproc', 'inputnode.std_t1w'),
            ('std_mask', 'inputnode.std_mask'),
        ]),
        (t1w_template_wf, anat_reports_wf, [
            ('outputnode.out_report', 'inputnode.t1w_conform_report'),
        ]),
        (anat_norm_wf, anat_reports_wf, [
            ('poutputnode.template', 'inputnode.template'),
        ]),
        # derivatives
        (t1w_template_wf, anat_derivatives_wf, [
            ('outputnode.t1w_valid_list', 'inputnode.source_files'),
            ('outputnode.t1w_realign_xfm', 'inputnode.t1w_ref_xfms'),
        ]),
        (be_buffer, anat_derivatives_wf, [
            ('anat_mask', 'inputnode.t1w_mask'),
            ('anat_preproc', 'inputnode.t1w_preproc'),
        ]),
        (anat_norm_wf, anat_derivatives_wf, [
            ('outputnode.template', 'inputnode.template'),
            ('outputnode.anat2std_xfm', 'inputnode.anat2std_xfm'),
            ('outputnode.std2anat_xfm', 'inputnode.std2anat_xfm'),
        ]),
        (anat_seg_wf, anat_derivatives_wf, [
            ('outputnode.anat_dseg', 'inputnode.t1w_dseg'),
            ('outputnode.anat_tpms', 'inputnode.t1w_tpms'),
        ]),
    ])

    if not freesurfer:
        return wf

    # FreeSurfer surfaces
    surface_recon_wf = init_infant_surface_recon_wf(
        age_months=age_months,
        use_aseg=segmentation_atlases is None,
    )

    wf.connect([
        (anat_seg_wf, surface_recon_wf, [
            ('outputnode.anat_aseg', 'inputnode.anat_aseg'),
        ]),
        (inputnode, surface_recon_wf, [
            ('subject_id', 'inputnode.subject_id'),
            ('subjects_dir', 'inputnode.subjects_dir'),
            ('t2w', 'inputnode.t2w'),
        ]),
        (anat_validate, surface_recon_wf, [
            ('out_file', 'inputnode.anat_orig'),
        ]),
        (be_buffer, surface_recon_wf, [
            ('anat_brain', 'inputnode.anat_skullstripped'),
            ('anat_preproc', 'inputnode.anat_preproc'),
        ]),
        (surface_recon_wf, outputnode, [
            ('outputnode.subjects_dir', 'subjects_dir'),
            ('outputnode.subject_id', 'subject_id'),
            ('outputnode.anat2fsnative_xfm', 'anat2fsnative_xfm'),
            ('outputnode.fsnative2anat_xfm', 'fsnative2anat_xfm'),
            ('outputnode.surfaces', 'surfaces'),
            ('outputnode.anat_aparc', 'anat_aparc'),
            ('outputnode.anat_aseg', 'anat_aseg'),
        ]),
        (surface_recon_wf, anat_reports_wf, [
            ('outputnode.subject_id', 'inputnode.subject_id'),
            ('outputnode.subjects_dir', 'inputnode.subjects_dir'),
        ]),
        (surface_recon_wf, anat_derivatives_wf, [
            ('outputnode.anat_aseg', 'inputnode.t1w_fs_aseg'),
            ('outputnode.anat_aparc', 'inputnode.t1w_fs_aparc'),
            ('outputnode.anat2fsnative_xfm', 'inputnode.t1w2fsnative_xfm'),
            ('outputnode.fsnative2anat_xfm', 'inputnode.fsnative2t1w_xfm'),
            ('outputnode.surfaces', 'inputnode.surfaces'),
        ]),
    ])
    # fmt: on
    return wf


def init_t2w_template_wf(longitudinal, omp_nthreads, num_t2w, name="anat_t2w_template_wf"):
    from pkg_resources import resource_filename as pkgr
    from niworkflows.interfaces.images import TemplateDimensions, Conform, ValidateImage

    wf = pe.Workflow(name=name)

    inputnode = pe.Node(niu.IdentityInterface(fields=["t2w"]), name="inputnode")
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=["t2w_ref", "t2w_valid_list", "t2_realign_xfm", "out_report"]),
        name="outputnode",
    )

    # 0. Reorient T1w image(s) to RAS and resample to common voxel space
    t2w_ref_dimensions = pe.Node(TemplateDimensions(), name='t2w_ref_dimensions')
    t2w_conform = pe.MapNode(Conform(), iterfield='in_file', name='t2w_conform')

    wf.connect([
        (inputnode, t2w_ref_dimensions, [('t2w', 't1w_list')]),
        (t2w_ref_dimensions, t2w_conform, [
            ('t1w_valid_list', 'in_file'),
            ('target_zooms', 'target_zooms'),
            ('target_shape', 'target_shape')]),
        (t2w_ref_dimensions, outputnode, [('out_report', 'out_report'),
                                          ('t1w_valid_list', 't2w_valid_list')]),
    ])

    # for now always take the first image
    # TODO: Average when multiple T2w images
    get1st = pe.Node(niu.Select(index=[0]), name='get1st')
    outputnode.inputs.t2w_realign_xfm = [pkgr('smriprep', 'data/itkIdentityTransform.txt')]

    wf.connect([
        (t2w_conform, get1st, [('out_file', 'inlist')]),
        (get1st, outputnode, [('out', 't2w_ref')]),
    ])

    return wf
