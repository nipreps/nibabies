"""
Subcortical alignment into MNI space
"""

from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl


def gen_subcortical_alignment_wf(name='subcortical_alignment_wf'):
    """
    Align individual subcortical structures into MNI space.

    This is a nipype workflow port of the DCAN infant pipeline.
    https://github.com/DCAN-Labs/dcan-infant-pipeline/blob/247e19e5441cc814cea2f23720caeeb6c6aeadf8/fMRISurface/scripts/SubcorticalAlign_ROIs.sh

    Parameters
    ----------
    rois : :obj:`list`
        ROIs to align
    name : :obj:`str`
        Name of the workflow

    Inputs
    ------
    std_xfm : :obj:`str`
        File containing transform to the standard (MNI) space

    Outputs
    -------

    """
    pass