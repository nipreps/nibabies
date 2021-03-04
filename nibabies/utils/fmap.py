# Patch to avoid any SDC correction
# TODO: Incorporate SDCFlows 2.0 into base workflow and remove this

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu
from nipype import logging

from niworkflows.engine.workflows import LiterateWorkflow as Workflow


def init_sdc_estimate_wf(fmaps, epi_meta, omp_nthreads=1, debug=False):
    """A legacy workflow to omit susceptibility distortion correction"""
    workflow = Workflow(name='sdc_bypass_wf')
    inputnode = pe.Node(niu.IdentityInterface(
        fields=['epi_file', 'epi_brain', 'epi_mask', 't1w_brain', 'std2anat_xfm']),
        name='inputnode')

    outputnode = pe.Node(niu.IdentityInterface(
        fields=['epi_corrected', 'epi_mask', 'epi_brain',
                'out_warp', 'syn_ref', 'method']),
        name='outputnode')

    # No fieldmaps - forward inputs to outputs
    workflow.__postdesc__ = """\
Susceptibility distortion correction (SDC) was omitted.
"""
    outputnode.inputs.method = 'None'
    outputnode.inputs.out_warp = 'identity'
    workflow.connect([
        (inputnode, outputnode, [('epi_file', 'epi_corrected'),
                                    ('epi_mask', 'epi_mask'),
                                    ('epi_brain', 'epi_brain')]),
    ])
    return workflow
