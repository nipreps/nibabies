# Use infant_recon_all to generate subcortical segmentations and cortical parcellations


def init_infant_surface_recon_wf(*, age_months, use_aseg=False, name="infant_surface_recon_wf"):
    from nipype.interfaces import freesurfer as fs
    from nipype.interfaces import utility as niu
    from nipype.pipeline import engine as pe
    from niworkflows.engine.workflows import LiterateWorkflow
    from niworkflows.interfaces.freesurfer import PatchedLTAConvert as LTAConvert
    from niworkflows.interfaces.freesurfer import (
        PatchedRobustRegister as RobustRegister,
    )
    from smriprep.workflows.surfaces import init_gifti_surface_wf

    from nibabies.interfaces.freesurfer import InfantReconAll

    wf = LiterateWorkflow(name=name)
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

    wf.__desc__ = f"""\
Brain surfaces were reconstructed using `infant_recon_all` [FreeSurfer
{fs.Info().looseversion() or "<ver>"}, RRID:SCR_001847, @infantfs],
leveraging the masked, preprocessed T1w and anatomical segmentation.
"""

    gen_recon_outdir = pe.Node(niu.Function(function=_gen_recon_dir), name="gen_recon_outdir")

    # inject the intensity-normalized skull-stripped t1w from the brain extraction workflow
    recon = pe.Node(InfantReconAll(age=age_months), name="reconall")

    fsnative2anat_xfm = pe.Node(
        niu.Function(function=_create_identity_lta),
        name="fsnative2anat_xfm",
    )

    anat2fsnative_xfm = pe.Node(
        LTAConvert(out_lta=True, invert=True),
        name="anat2fsnative_xfm",
    )

    # convert generated surfaces to GIFTIs
    gifti_surface_wf = init_gifti_surface_wf()

    get_aseg = pe.Node(niu.Function(function=_get_aseg), name="get_aseg")
    get_aparc = pe.Node(niu.Function(function=_get_aparc), name="get_aparc")
    aparc2nii = pe.Node(fs.MRIConvert(out_type="niigz"), name="aparc2nii")

    if use_aseg:
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
        (inputnode, fsnative2anat_xfm, [('anat_skullstripped', 'in_file')]),
        (fsnative2anat_xfm, anat2fsnative_xfm, [('out', 'in_lta')]),
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
        (fsnative2anat_xfm, outputnode, [
            ('out', 'fsnative2anat_xfm'),
        ]),
        (anat2fsnative_xfm, outputnode, [
            ('out_lta', 'anat2fsnative_xfm'),
        ]),
        (fsnative2anat_xfm, gifti_surface_wf, [
            ('out', 'inputnode.fsnative2t1w_xfm')]),
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


def _create_identity_lta(in_file):
    """`infant_recon_all` will use the masked T1w as fsnative"""
    from pathlib import Path

    import nitransforms as nt
    import numpy as np

    outfile = Path('identity.lta')
    xfm = nt.Affine(np.eye(4), in_file)
    xfm.to_filename(str(outfile), fmt='lta')
    return outfile


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
