# Use infant_recon_all to generate subcortical segmentations and cortical parcellations

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, freesurfer as fs
from niworkflows.interfaces.freesurfer import (
    PatchedLTAConvert as LTAConvert,
)
from smriprep.workflows.surfaces import init_gifti_surface_wf

from ...interfaces.freesurfer import InfantReconAll


def init_infant_surface_recon_wf(*, age_months, use_aseg=False, name="infant_surface_recon_wf"):
    wf = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "subjects_dir",
                "subject_id",
                "anat_orig",
                "anat_skullstripped",
                "anat_preproc",
                "anat_aseg",
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
                "anat_aseg",
                "anat_aparc",
            ]
        ),
        name="outputnode",
    )

    gen_recon_outdir = pe.Node(niu.Function(function=_gen_recon_dir), name="gen_recon_outdir")

    # inject the intensity-normalized skull-stripped t1w from the brain extraction workflow
    recon = pe.Node(InfantReconAll(age=age_months), name="reconall")

    # these files are created by babyFS, but transforms are for masked anatomicals
    # https://github.com/freesurfer/freesurfer/blob/
    # 8b40551f096294cc6603ce928317b8df70bce23e/infant/infant_recon_all#L744
    # TODO: calculate full anat -> fsnative transform?
    get_tal_lta = pe.Node(
        niu.Function(function=_get_talairch_lta),
        name="get_tal_xfm",
    )
    fsnative2anat_xfm = pe.Node(
        LTAConvert(out_lta=True, invert=True),
        name="fsnative2anat_xfm",
    )

    # convert generated surfaces to GIFTIs
    gifti_surface_wf = init_gifti_surface_wf()

    get_aseg = pe.Node(niu.Function(function=_get_aseg), name="get_aseg")
    get_aparc = pe.Node(niu.Function(function=_get_aparc), name="get_aparc")
    aparc2nii = pe.Node(fs.MRIConvert(out_type="niigz"), name="aparc2nii")

    if use_aseg:
        # TODO: Add precomputed segmentation upon new babyFS rel
        wf.connect(inputnode, "anat_aseg", recon, "aseg_file")

    # fmt: off
    wf.connect([
        (inputnode, gen_recon_outdir, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id'),
        ]),
        (inputnode, recon, [
            ('anat_skullstripped', 'mask_file'),
            ('subject_id', 'subject_id'),
        ]),
        (gen_recon_outdir, recon, [
            ('out', 'outdir'),
        ]),
        (recon, outputnode, [
            ('subject_id', 'subject_id'),
            (('outdir', _parent), 'subjects_dir'),
        ]),
        (recon, gifti_surface_wf, [
            ('subject_id', 'inputnode.subject_id'),
            (('outdir', _parent), 'inputnode.subjects_dir'),
        ]),
        (recon, get_aparc, [
            ('outdir', 'fs_subject_dir'),
        ]),
        (recon, get_aseg, [
            ('outdir', 'fs_subject_dir'),
        ]),
        (get_aseg, outputnode, [
            ('out', 'anat_aseg'),
        ]),
        (get_aparc, aparc2nii, [
            ('out', 'in_file'),
        ]),
        (aparc2nii, outputnode, [
            ('out_file', 'anat_aparc'),
        ]),
        (recon, get_tal_lta, [
            ('outdir', 'fs_subject_dir'),
        ]),
        (get_tal_lta, outputnode, [
            ('out', 'anat2fsnative_xfm'),
        ]),
        (get_tal_lta, fsnative2anat_xfm, [
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


def _parent(p):
    from pathlib import Path

    return str(Path(p).parent)


def _gen_recon_dir(subjects_dir, subject_id):
    from pathlib import Path

    p = Path(subjects_dir) / subject_id
    p.mkdir(parents=True, exist_ok=True)
    return str(p)


def _get_talairch_lta(fs_subject_dir):
    """Fetch pre-computed transform from infant_recon_all"""
    from pathlib import Path

    xfm = Path(fs_subject_dir) / "mri" / "transforms" / "niftyreg_affine.lta"
    if not xfm.exists():
        raise FileNotFoundError("Could not find talairach transform.")
    return str(xfm.absolute())


def _get_aseg(fs_subject_dir):
    """Fetch infant_recon_all's aparc+aseg"""
    from pathlib import Path

    aseg = Path(fs_subject_dir) / "mri" / "aseg.nii.gz"
    if not aseg.exists():
        raise FileNotFoundError("Could not find aseg.")
    return str(aseg)


def _get_aparc(fs_subject_dir):
    """Fetch infant_recon_all's aparc+aseg"""
    from pathlib import Path

    aparc = Path(fs_subject_dir) / "mri" / "aparc+aseg.mgz"
    if not aparc.exists():
        raise FileNotFoundError("Could not find aparc.")
    return str(aparc)
