# Use infant_recon_all to generate subcortical segmentations and cortical parcellations

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from niworkflows.interfaces.nibabel import ApplyMask
from niworkflows.interfaces.freesurfer import (
    PatchedRobustRegister as RobustRegister,
    PatchedLTAConvert as LTAConvert,
)
from smriprep.workflows.surfaces import init_gifti_surface_wf

from ...interfaces.freesurfer import InfantReconAll


def init_infant_surface_recon_wf(
    *,
    age_months,
    name='infant_surface_recon_wf'
):
    wf = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "subjects_dir",
                "subject_id",
                "anat_orig",
                "anat_skullstripped",
                "anat_preproc",
            ],
        ),
        name="inputnode",
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["subjects_dir", "subject_id", "surfaces"]),
        name="outputnode",
    )

    # inject the intensity-normalized skull-stripped t1w from the brain extraction workflow
    recon = pe.Node(
        InfantReconAll(age=age_months),
        name='reconall'
    )

    # these files are created by babyFS, but transforms are for masked anatomicals
    # https://github.com/freesurfer/freesurfer/blob/8b40551f096294cc6603ce928317b8df70bce23e/infant/infant_recon_all#L743-L788
    # calculate anat -> fsnative transform
    fsnative2anat_xfm = pe.Node(RobustRegister(auto_sens=True, est_int_scale=True),
                            name='fsnative2t1w_xfm')
    anat2fsnative_xfm = pe.Node(LTAConvert(out_lta=True, invert=True),
                               name='t1w2fsnative_xfm')

    # convert generated surfaces to GIFTIs
    gifti_surface_wf = init_gifti_surface_wf()

    # fmt: off
    wf.connect([
        (inputnode, recon, [
            ('anat_skullstripped', 'mask_file'),
            ('subjects_dir', 'outdir'),
            ('subject_id', 'subject_id'),
            # ('in_seg', 'aseg_file'),  # TODO: Add precomputed segmentation when upgrading version
        ]),
        (recon, outputnode, [
            ('outdir', 'subjects_dir'),
        ]),
        (inputnode, fsnative2t1w_xfm, [
            ('anat_orig', 'target_file'),
        ]),
        (inputnode, gifti_surface_wf, [
            ('subject_id', 'inputnode.subject_id'),
            ('subjects_dir', 'inputnode.subjects_dir'),
        ]),
        (gifti_surface_wf, outputnode, [
            ('outputnode.surfaces', 'surfaces'),
        ]),
    ])
    # fmt: on
    return wf
