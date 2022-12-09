# Use infant_recon_all to generate subcortical segmentations and cortical parcellations

from nipype.interfaces import fsl
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe

from ...config import DEFAULT_MEMORY_MIN_GB
from ...interfaces.workbench import CreateSignedDistanceVolume


def init_infant_surface_recon_wf(*, age_months, use_aseg=False, name="infant_surface_recon_wf"):
    from nipype.interfaces import freesurfer as fs
    from nipype.interfaces import io as nio
    from niworkflows.engine.workflows import LiterateWorkflow
    from niworkflows.interfaces.freesurfer import PatchedLTAConvert as LTAConvert
    from niworkflows.interfaces.freesurfer import (
        PatchedRobustRegister as RobustRegister,
    )
    from smriprep.workflows.surfaces import init_gifti_surface_wf

    from nibabies.interfaces.freesurfer import InfantReconAll

    # Synchronized inputs to smriprep.workflows.surfaces.init_surface_recon_wf
    wf = LiterateWorkflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "subjects_dir",
                "subject_id",
                "t1w",
                "t2w",
                "flair",
                "skullstripped_t1",
                "corrected_t1",
                "ants_segs",
            ],
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "subjects_dir",
                "subject_id",
                "t1w2fsnative_xfm",
                "fsnative2t1w_xfm",
                "surfaces",
                "out_aseg",
                "out_aparc",
            ]
        ),
        name="outputnode",
    )

    wf.__desc__ = f"""\
Brain surfaces were reconstructed using `infant_recon_all` [FreeSurfer
{fs.Info().looseversion() or "<ver>"}, RRID:SCR_001847, @infantfs],
leveraging the masked, preprocessed T1w and anatomical segmentation.
"""

    gen_recon_outdir = pe.Node(niu.Function(function=_gen_recon_dir), name="gen_recon_outdir")

    # inject the intensity-normalized skull-stripped t1w from the brain extraction workflow
    recon = pe.Node(InfantReconAll(age=age_months), name="reconall")
    fssource = pe.Node(nio.FreeSurferSource(), name='fssource', run_without_submitting=True)

    fsnative2t1w_xfm = pe.Node(
        RobustRegister(auto_sens=True, est_int_scale=True),
        name='fsnative2t1w_xfm',
    )

    t1w2fsnative_xfm = pe.Node(
        LTAConvert(out_lta=True, invert=True),
        name="t1w2fsnative_xfm",
    )

    # convert generated surfaces to GIFTIs
    gifti_surface_wf = init_gifti_surface_wf()

    aparc2nii = pe.Node(fs.MRIConvert(out_type="niigz"), name="aparc2nii")

    if use_aseg:
        wf.connect(inputnode, "ants_segs", recon, "aseg_file")

    # fmt: off
    wf.connect([
        (inputnode, gen_recon_outdir, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id'),
        ]),
        (inputnode, recon, [
            ('skullstripped_t1', 'mask_file'),
            ('subject_id', 'subject_id'),
        ]),
        (gen_recon_outdir, recon, [
            ('out', 'outdir'),
        ]),
        (recon, outputnode, [
            ('subject_id', 'subject_id'),
            (('outdir', _parent), 'subjects_dir'),
        ]),
        (recon, fssource, [
            ('subject_id', 'subject_id'),
            (('outdir', _parent), 'subjects_dir'),
        ]),
        (recon, gifti_surface_wf, [
            ('subject_id', 'inputnode.subject_id'),
            (('outdir', _parent), 'inputnode.subjects_dir'),
        ]),
        (fssource, outputnode, [
            (('aseg', _replace_mgz), 'anat_aseg'),
        ]),
        (inputnode, fsnative2t1w_xfm, [('skullstripped_t1', 'target_file')]),
        (fssource, fsnative2t1w_xfm, [
            (('norm', _replace_mgz), 'source_file'),
        ]),
        (fsnative2t1w_xfm, t1w2fsnative_xfm, [('out_reg_file', 'in_lta')]),
        (fssource, aparc2nii, [
            ('aparc_aseg', 'in_file'),
        ]),
        (aparc2nii, outputnode, [
            ('out_file', 'out_aparc'),
        ]),
        (fssource, outputnode, [
            (('aseg', _replace_mgz), 'out_aseg'),
        ]),
        (fsnative2t1w_xfm, outputnode, [
            ('out_reg_file', 'fsnative2t1w_xfm'),
        ]),
        (t1w2fsnative_xfm, outputnode, [
            ('out_lta', 't1w2fsnative_xfm'),
        ]),
        (fsnative2t1w_xfm, gifti_surface_wf, [
            ('out_reg_file', 'inputnode.fsnative2t1w_xfm')]),
        (gifti_surface_wf, outputnode, [
            ('outputnode.surfaces', 'surfaces'),
        ]),
    ])
    # fmt: on
    return wf


def init_anat_ribbon_wf(name="anat_ribbon_wf"):
    # 0, 1 = wm; 2, 3 = pial; 6, 7 = mid
    # note that order of lh / rh within each surf type is not guaranteed due to use
    # of unsorted glob by FreeSurferSource prior, but we can do a sort
    # to ensure consistent ordering
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "surfaces",  # anat_giftis,
                "t1w_mask",
            ]
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "anat_ribbon",
            ]
        ),
        name="outputnode",
    )

    select_wm = pe.Node(
        niu.Select(index=[0, 1]),
        name="select_wm",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    select_pial = pe.Node(
        niu.Select(index=[2, 3]),
        name="select_pial",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    select_midthick = pe.Node(
        niu.Select(index=[6, 7]),
        name="select_midthick",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    create_wm_distvol = pe.MapNode(
        CreateSignedDistanceVolume(),
        iterfield=["surf_file"],
        name="create_wm_distvol",
    )

    create_pial_distvol = pe.MapNode(
        CreateSignedDistanceVolume(),
        iterfield=["surf_file"],
        name="create_pial_distvol",
    )

    thresh_wm_distvol = pe.MapNode(
        fsl.maths.MathsCommand(args="-thr 0 -bin -mul 255"),
        iterfield=["in_file"],
        name="thresh_wm_distvol",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    uthresh_pial_distvol = pe.MapNode(
        fsl.maths.MathsCommand(args="-uthr 0 -abs -bin -mul 255"),
        iterfield=["in_file"],
        name="uthresh_pial_distvol",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    bin_wm_distvol = pe.MapNode(
        fsl.maths.UnaryMaths(operation="bin"),
        iterfield=["in_file"],
        name="bin_wm_distvol",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    bin_pial_distvol = pe.MapNode(
        fsl.maths.UnaryMaths(operation="bin"),
        iterfield=["in_file"],
        name="bin_pial_distvol",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    split_wm_distvol = pe.Node(
        niu.Split(splits=[1, 1]),
        name="split_wm_distvol",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    merge_wm_distvol_no_flatten = pe.Node(
        niu.Merge(2),
        no_flatten=True,
        name="merge_wm_distvol_no_flatten",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    make_ribbon_vol = pe.MapNode(
        fsl.maths.MultiImageMaths(op_string="-mas %s -mul 255"),
        iterfield=["in_file", "operand_files"],
        name="make_ribbon_vol",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    bin_ribbon_vol = pe.MapNode(
        fsl.maths.UnaryMaths(operation="bin"),
        iterfield=["in_file"],
        name="bin_ribbon_vol",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    split_squeeze_ribbon_vol = pe.Node(
        niu.Split(splits=[1, 1], squeeze=True),
        name="split_squeeze_ribbon",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    combine_ribbon_vol_hemis = pe.Node(
        fsl.maths.BinaryMaths(operation="add"),
        name="combine_ribbon_vol_hemis",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    # make HCP-style ribbon volume in T1w space
    workflow.connect(
        [
            (inputnode, select_wm, [("surfaces", "inlist")]),
            (inputnode, select_pial, [("surfaces", "inlist")]),
            (inputnode, select_midthick, [("surfaces", "inlist")]),
            (select_wm, create_wm_distvol, [(("out", _sorted_by_basename), "surf_file")]),
            (inputnode, create_wm_distvol, [("t1w_mask", "ref_file")]),
            (select_pial, create_pial_distvol, [(("out", _sorted_by_basename), "surf_file")]),
            (inputnode, create_pial_distvol, [("t1w_mask", "ref_file")]),
            (create_wm_distvol, thresh_wm_distvol, [("out_file", "in_file")]),
            (create_pial_distvol, uthresh_pial_distvol, [("out_file", "in_file")]),
            (thresh_wm_distvol, bin_wm_distvol, [("out_file", "in_file")]),
            (uthresh_pial_distvol, bin_pial_distvol, [("out_file", "in_file")]),
            (bin_wm_distvol, split_wm_distvol, [("out_file", "inlist")]),
            (split_wm_distvol, merge_wm_distvol_no_flatten, [("out1", "in1")]),
            (split_wm_distvol, merge_wm_distvol_no_flatten, [("out2", "in2")]),
            (bin_pial_distvol, make_ribbon_vol, [("out_file", "in_file")]),
            (merge_wm_distvol_no_flatten, make_ribbon_vol, [("out", "operand_files")]),
            (make_ribbon_vol, bin_ribbon_vol, [("out_file", "in_file")]),
            (bin_ribbon_vol, split_squeeze_ribbon_vol, [("out_file", "inlist")]),
            (split_squeeze_ribbon_vol, combine_ribbon_vol_hemis, [("out1", "in_file")]),
            (split_squeeze_ribbon_vol, combine_ribbon_vol_hemis, [("out2", "operand_file")]),
            (combine_ribbon_vol_hemis, outputnode, [("out_file", "anat_ribbon")]),
        ]
    )
    return workflow


def _parent(p):
    from pathlib import Path

    return str(Path(p).parent)


def _gen_recon_dir(subjects_dir, subject_id):
    from pathlib import Path

    p = Path(subjects_dir) / subject_id
    p.mkdir(parents=True, exist_ok=True)
    return str(p)


def _replace_mgz(in_file):
    return in_file.replace('.mgz', '.nii.gz')


def _sorted_by_basename(inlist):
    from os.path import basename

    return sorted(inlist, key=lambda x: str(basename(x)))
