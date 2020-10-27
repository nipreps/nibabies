from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl

def init_infant_anat_wf(
    *,
    age_months,
    anatomicals,
    anat_modality,
    bids_root,
    existing_derivatives,
    freesurfer,
    omp_nthreads,
    output_dir,
    segmentation_atlases,
    skull_strip_template,
    sloppy,
    spaces,
    name='infant_anat_wf',
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
    t1w2fsnative_xfm
        LTA-style affine matrix translating from T1w to
        FreeSurfer-conformed subject space
    fsnative2t1w_xfm
        LTA-style affine matrix translating from FreeSurfer-conformed
        subject space to T1w
    surfaces
        GIFTI surfaces (gray/white boundary, midthickness, pial, inflated)
    """

    from smriprep.workflows.anatomical import init_anat_template_wf, _probseg_fast2bids
    from smriprep.workflows.norm import init_anat_norm_wf
    from smriprep.workflows.outputs import init_anat_reports_wf, init_anat_derivatives_wf

    from .brain_extraction import init_infant_brain_extraction_wf
    from .surfaces import init_infant_surface_recon_wf

    # for now, T1w only
    num_anats = len(anatomicals)
    wf = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["t1w", "t2w", "subject_id", "subjects_dir"]), # FLAIR / ROI?
        name='inputnode'
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "anat_preproc", "anat_brain", "anat_mask", "anat_dseg", "anat_tpms",
                "std_preproc", "std_brain", "std_dseg", "std_tpms",
                "subjects_dir", "subject_id", "anat2std_xfm", "std2anat_xfm",
                "anat2fsnative_xfm", "fsnative2anat_xfm", "surfaces",
            ]
        ), name='outputnode')

    # Define output workflows
    anat_reports_wf = init_anat_reports_wf(freesurfer=freesurfer, output_dir=output_dir)
    anat_derivatives_wf = init_anat_derivatives_wf(
        bids_root=bids_root,
        freesurfer=freesurfer,
        num_t1w=num_anats,
        output_dir=output_dir,
        spaces=spaces,
    )

    # Multiple T1w files -> generate average reference
    anat_template_wf = init_anat_template_wf(
        longitudial=longitudial,
        omp_nthreads=omp_nthreads,
        num_t1w=num_anat,
    )

    # INU + Brain Extraction
    brain_extraction_wf = init_infant_brain_extraction_wf(
        age_months=age_months,
        anat_modality=anat_modality,
        ants_affine_init=True,
        skull_strip_template=skull_strip_template,
        template_specs=template_specs,
        omp_nthreads=omp_nthreads,
        output_dir=output_dir,
        sloppy=sloppy,
    )

    # Segmentation - initial implementation should be simple: JLF
    anat_seg_wf = init_anat_seg_wf(
        age_months=age_months,
        anat_modality=anat_modality,
        template_dir=segmentation_atlases,
        sloppy=sloppy,
        omp_nthreads=omp_nthreads,
    )

    # Spatial normalization (requires segmentation)
    anat_norm_wf = init_anat_norm_wf(
        debug=debug,
        omp_nthreads=omp_nthreads,
        templates=spaces.get_spaces(nonstandard=False, dim=(3,)),
    )

    if not freesurfer:  # Flag --fs-no-reconall is set - return
        # This should be handled by anat_seg_wf
        # Brain tissue segmentation - FAST produces: 0 (bg), 1 (wm), 2 (csf), 3 (gm)
        t1w_dseg = pe.Node(fsl.FAST(segments=True, no_bias=True, probability_maps=True),
                           name='t1w_dseg', mem_gb=3)
        lut_t1w_dseg.inputs.lut = (0, 3, 1, 2)  # Maps: 0 -> 0, 3 -> 1, 1 -> 2, 2 -> 3.
        fast2bids = pe.Node(niu.Function(function=_probseg_fast2bids), name="fast2bids",
                            run_without_submitting=True)

        # fmt: off
        workflow.connect([
            (brain_extraction_wf, buffernode, [
                (('outputnode.out_file', _pop), 'anat_brain'),
                ('outputnode.out_mask', 'anat_mask')]),
            (buffernode, t1w_dseg, [('t1w_brain', 'in_files')]),
            (t1w_dseg, lut_t1w_dseg, [('partial_volume_map', 'in_dseg')]),
            (t1w_dseg, fast2bids, [('partial_volume_files', 'inlist')]),
            (fast2bids, anat_norm_wf, [('out', 'inputnode.moving_tpms')]),
            (fast2bids, outputnode, [('out', 'anat_tpms')]),
        ])
        # fmt: on
        return workflow

    # FreeSurfer surfaces
    surface_recon_wf = init_infant_surface_recon_wf(age_months=age_months)
    applyrefined = pe.Node(fsl.ApplyMask(), name='applyrefined')
    # fmt: off
    wf.connect([
        (inputnode, brain_extraction_wf, [
            ('t1w', 'inputnode.t1w'),
            ('t2w', 'inputnode.t2w')]),
        (brain_extraction_wf, outputnode, [
            ('outputnode.out_corrected', 'anat_preproc'),
            ('outputnode.out_brain', 'anat_brain'),
            ('outputnode.out_mask', 'anat_mask')]),
        (brain_extraction_wf, surface_recon_wf, [
            ('outputnode.out_brain', 'inputnode.masked_file')]),
        (surface_recon_wf, outputnode, [
            ('outputnode.subjects_dir', 'subjects_dir')]),
        (surface_recon_wf, anat_reports_wf, [
            ('outputnode.subjects_dir', 'subjects_dir'),
            ('outputnode.subject_id', 'subject_id')]),
    ])
    # fmt: on
    return wf


def init_anat_seg_wf(
    age_months=None,
    anat_modality='T1w',
    template_dir=None,
    sloppy=False,
    omp_nthreads=None,
    name='anat_seg_wf'
):
    """Calculate segmentation from collection of OHSU atlases"""
    from pkg_resources import resource_filename as pkgr_fn
    from nipype.interfaces.ants.segmentation import JointFusion
    from niworkflows.interfaces.fixes import (
        FixHeaderRegistration as Registration,
        FixHeaderApplyTransforms as ApplyTransforms,
    )
    from smriprep.utils.misc import apply_lut

    if template_dir is None:  # TODO: set default
        raise NotImplementedError("No default has been set yet. Please pass segmentations.")
    if anat_modality != 'T1w':
        raise NotImplementedError("Only T1w images are currently accepted for ANTs LabelFusion.")

    wf = pe.Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=["anat_brain"]), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=["anat_dseg"]), name='outputnode')

    tmpl_anats, tmpl_segs = _parse_segmentation_atlases(anat_modality, template_dir)

    # register to templates
    ants_params = "testing" if sloppy else "precise"
    # Register to each subject space
    norm = pe.MapNode(
        Registration(from_file=pkgr_fn(
            "niworkflows.data",
            f"antsBrainExtraction_{ants_params}.json")
        ),
        name="norm",
        iterfield=['fixed_image', 'moving_image'],
        n_procs=omp_nthreads,
        mem_gb=mem_gb,
    )
    norm.inputs.moving_image = tmpl_anats
    norm.inputs.float = True

    apply_atlas = pe.MapNode(
        ApplyTransforms(
            dimension=3,
            interpolation="NearestNeighbor",
            float=True,
        ),
        iterfield=['transforms', 'input_image'],
        name='apply_atlas',
    )
    apply_atlas.inputs.input_image = tmpl_anats

    apply_seg = pe.Node(
        ApplyTransforms(dimension=3, interpolation="MultiLabel"),  # NearestNeighbor?
        name='apply_seg',
        iterfield=['transforms', 'input_image'],
    )
    apply_seg.inputs.input_image = tmpl_segs

    jointfusion = pe.Node(
        JointFusion(
            dimension=3,
            out_label_fusion='fusion_labels.nii.gz',
        ),
        name='jointfusion'
    )

    def _to_list(x):
        return [x]

    # fmt: off
    wf.connect([
        (inputnode, reg, [('anat_brain', 'fixed_image')]),
        (reg, apply_atlas, [('forward_transforms', 'transforms')]),
        (inputnode, apply_atlas, [('anat_brain', 'reference_image')]),
        (reg, apply_seg, [('forward_transforms', 'transforms')]),
        (inputnode, apply_seg, [('anat_brain', 'reference_image')]),
        (inputnode, jointfusion, [(('anat_brain', _to_list), 'target_image')]),
        (apply_atlas, jointfusion, [('output_image', 'atlas_image')]),
        (apply_seg, jointfusion, [('output_image', 'atlas_segmentation_image')]),
    ])
    # fmt: on

    return wf


def _parse_segmentation_atlases(anat_modality, template_dir):
    """
    Parse segmentation templates directory for anatomical and segmentation files.

    This is currently hardcoded to match DCAN lab templates.
    Will need to rethink standardization for more general cases.
    """
    from pathlib import Path

    anats, segs = [], []

    for f in Path(template_dir).glob('**/*.nii*'):
        if 'Segmentation' in f.name:
            segs.append(f.absolute())
        elif anat_modality in f.name:
            anats.append(f.absolute())

    assert anats
    assert segs
    # there should matching files per template
    assert len(anats) == len(segs)

    return sorted(anats), sorted(segs)
