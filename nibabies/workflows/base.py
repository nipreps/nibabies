from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from .brain_extraction import init_infant_brain_extraction_wf
from .surfaces import init_infant_surface_recon_wf

def init_infant_anat_wf():
    inputnode = pe.Node(niu.IdentityInterface(fields=["in_files"]), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=["anat_corrected", "anat_brain", "anat_mask", "surfaces", "anat_aseg", "anat_aparc"]
    ), name='outputnode')