# Use infant_recon_all to generate subcortical segmentations and cortical parcellations

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from niworkflows.interfaces.nibabel import ApplyMask
from niworkflows.interfaces.freesurfer import PatchedRobustRegister as RobustRegister
from smriprep.workflows.surfaces import init_gifti_surface_wf

from ..interfaces.freesurfer import InfantReconAll


def init_infant_surface_recon_wf(
    *,
    age_months,
    output_dir,
    subject_id,
    name='infant_surface_recon_wf'
):
    wf = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=["masked_file", "t1_file"]), name="inputnode"
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["subjects_dir"]), name="outputnode"
    )

    # we will use the intensity-normalized t1w from the brain extraction workflow with
    # the brainmask applied, and then feed that into baby freesurfer
    # TO TEST: injecting another segmentation to replace baby FS's aseg.
    recon = pe.Node(
        InfantReconAll(age=age_months, outdir=output_dir, subject_id=subject_id),
        name='reconall'
    )
    wf.connect([
        (inputnode, recon, [
            ('masked_file', 'mask_file'),
            # ('in_seg', 'aseg_file'),f
        ]),
        (recon, outputnode, [('outdir', 'subjects_dir')])
    ])

    # convert generated surfaces to GIFTIs
    # gifti_surface_wf = init_gifti_surface_wf()
    # fsnative2t1w_xfm = pe.Node(RobustRegister(auto_sens=True, est_int_scale=True),
    #                            name='fsnative2t1w_xfm')
    return wf
