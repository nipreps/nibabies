# Use infant_recon_all to generate subcortical segmentations and cortical parcellations

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, freesurfer as fs
from niworkflows.interfaces.nibabel import ApplyMask
from niworkflows.interfaces.freesurfer import (
    PatchedRobustRegister as RobustRegister,
    PatchedLTAConvert as LTAConvert,
)
from smriprep.workflows.surfaces import init_gifti_surface_wf

from ...interfaces.freesurfer import InfantReconAll


def init_infant_surface_recon_wf(*, age_months, name="infant_surface_recon_wf"):
    wf = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "subjects_dir",
                "subject_id",
                "anat_orig",
                "anat_skullstripped",
                "anat_preproc",
                "anat_seg",
                "t2w",
            ],
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "subjects_dir",
                "subject_id",
                "anat2fsnative_xfm",
                "fsnative2anat_xfm",
                "surfaces",
                "out_aparc",
            ]
        ),
        name="outputnode",
    )

    # inject the intensity-normalized skull-stripped t1w from the brain extraction workflow
    recon = pe.Node(InfantReconAll(age=age_months), name="reconall")

    # these files are created by babyFS, but transforms are for masked anatomicals
    # https://github.com/freesurfer/freesurfer/blob/
    # 8b40551f096294cc6603ce928317b8df70bce23e/infant/infant_recon_all#L744
    # TODO: calculate full anat -> fsnative transform?
    anat2fsnative_xfm = pe.Node(
        niu.Function(function=_get_talairch_xfm), name="anat2fsnative_xfm"
    )
    fsnative2anat_xfm = pe.Node(
        LTAConvert(out_lta=True, invert=True), name="fsnative2anat_xfm"
    )

    # convert generated surfaces to GIFTIs
    gifti_surface_wf = init_gifti_surface_wf()

    get_aparc = pe.Node(niu.Function(function=_get_aparc), name="get_aparc")
    aparc2nii = pe.Node(fs.MRIConvert(out_type="niigz"), name="aparc2nii")

    # fmt: off
    wf.connect([
        (inputnode, recon, [
            ('anat_skullstripped', 'mask_file'),
            ('subjects_dir', 'outdir'),
            ('subject_id', 'subject_id'),
            # ('anat_seg', 'aseg_file'),  # TODO: Add precomputed segmentation when upgrading version
        ]),
        (recon, outputnode, [
            ('outdir', 'subjects_dir'),
        ]),
        (recon, gifti_surface_wf, [
            ('subject_id', 'inputnode.subject_id'),
            ('outdir', 'inputnode.subjects_dir'),
        ]),
        (recon, get_aparc, [
            ('subject_id', 'subject_id'),
            ('outdir', 'subjects_dir'),
        ]),
        (get_aparc, outputnode, [
            ('out', 'out_aparc'),
        ]),
        (recon, anat2fsnative_xfm, [
            ('subject_id', 'subject_id'),
            ('outdir', 'subjects_dir'),
        ]),
        (anat2fsnative_xfm, outputnode, [
            ('out', 'anat2fsnative_xfm'),
        ]),
        (anat2fsnative_xfm, fsnative2anat_xfm, [
            ('out', 'in_lta'),
        ]),
        (fsnative2anat_xfm, outputnode, [
            ('out_lta', 'fsnative2anat_xfm'),
        ]),
        (fsnative2anat_xfm, gifti_surface_wf, [
            ('out_lta', 'inputnode.fsnative2t1w_xfm')]),
        (gifti_surface_wf, outputnode, [
            ('outputnode.surfaces', 'surfaces'),
        ]),
    ])
    # fmt: on
    return wf


def _get_talairch_xfm(subject_id, subjects_dir):
    """Fetch pre-computed transform from infant_recon_all"""
    from pathlib import Path

    xfm = Path(subjects_dir) / subject_id / "mri" / "transforms" / "talairach.xfm"
    if not xfm.exists():
        raise FileNotFoundError("Could not find talairach transform.")
    return str(xfm.absolute())


def _get_aparc(subject_id, subjects_dir):
    """Fetch infant_recon_all's aparc+aseg"""
    from pathlib import Path

    aparc = Path(subjects_dir) / subject_id / "mri" / "aparc+aseg.mgz"
    if not aparc.exists():
        raise FileNotFoundError("Could not find aparc.")
    return str(aparc)
