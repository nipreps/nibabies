# Use infant_recon_all to generate subcortical segmentations and cortical parcellations


def init_infant_surface_recon_wf(*, age_months, use_aseg=False, name="infant_surface_recon_wf"):
    from nipype.interfaces import freesurfer as fs
    from nipype.interfaces import io as nio
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
    fssource = pe.Node(nio.FreeSurferSource(), name='fssource', run_without_submitting=True)

    fsnative2anat_xfm = pe.Node(
        RobustRegister(auto_sens=True, est_int_scale=True),
        name='fsnative2anat_xfm',
    )

    anat2fsnative_xfm = pe.Node(
        LTAConvert(out_lta=True, invert=True),
        name="anat2fsnative_xfm",
    )

    # convert generated surfaces to GIFTIs
    gifti_surface_wf = init_gifti_surface_wf()

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
        (inputnode, fsnative2anat_xfm, [('anat_skullstripped', 'target_file')]),
        (fssource, fsnative2anat_xfm, [
            (('norm', _replace_mgz), 'source_file'),
        ]),
        (fsnative2anat_xfm, anat2fsnative_xfm, [('out_reg_file', 'in_lta')]),
        (fssource, aparc2nii, [
            ('aparc_aseg', 'in_file'),
        ]),
        (aparc2nii, outputnode, [
            ('out_file', 'anat_aparc'),
        ]),
        (fsnative2anat_xfm, outputnode, [
            ('out_reg_file', 'fsnative2anat_xfm'),
        ]),
        (anat2fsnative_xfm, outputnode, [
            ('out_lta', 'anat2fsnative_xfm'),
        ]),
        (fsnative2anat_xfm, gifti_surface_wf, [
            ('out_reg_file', 'inputnode.fsnative2t1w_xfm')]),
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


def _replace_mgz(in_file):
    return in_file.replace('.mgz', '.nii.gz')
