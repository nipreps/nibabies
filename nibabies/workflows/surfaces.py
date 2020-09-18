# Use infant_recon_all to generate subcortical segmentations and cortical parcellations

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from niworkflows.interfaces.nibabel import ApplyMask
from niworkflows.interfaces.freesurfer import PatchedRobustRegister as RobustRegister
from smriprep.workflows.surfaces import init_gifti_surface_wf

from ..interfaces.freesurfer import InfantReconAll


def init_infant_surface_recon_wf(age_months):
    inputnode = pe.Node(
        niu.IdentityInterface(fields=["in_masked", "t1_file", "in_seg"]), name="inputnode"
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=["t1w_aseg", "t1w_aparc"]), name="outputnode"
    )

    # we will use the intensity-normalized t1w from the brain extraction workflow with
    # the brainmask applied, and then feed that into baby freesurfer
        # down the line, we might want to inject another segmentation to replace baby FS's aseg.
    recon = pe.Node(InfantReconAll(age=age_months), name='reconall')
    wf.connect([
        (inputnode, recon, [('t1_file', 't1_file'),
                            ('masked_file', 'mask_file')]),
    ])

    fsnative2t1w_xfm = pe.Node(RobustRegister(auto_sens=True, est_int_scale=True),
                               name='fsnative2t1w_xfm')



