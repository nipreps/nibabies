# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
#
# Copyright The NiPreps Developers <nipreps@gmail.com>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
"""
Session-level BOLD workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: init_bold_coreg_runs_wf

"""

from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.func.util import init_skullstrip_bold_wf


def init_coreg_bolds_wf(
    *,
    num_bold_runs: int,
    unbiased: bool = False,
    omp_nthreads: int = 1,
    name: str = 'coreg_bolds_wf',
) -> Workflow:
    """
    Register all BOLD runs of a session to a single session-reference.

    Parameters
    ----------
    num_bolds : :obj:`int`
        Number of BOLD runs
    omp_nthreads : :obj:`int`
        Number of threads.
    unbiased : :obj:`bool`
        Whether to use an unbiased registration strategy (default: False).

    Inputs
    ------
    boldref_files
        List of BOLD reference files to be coregistered.

    Outputs
    -------
    boldref
        The computed session-level BOLD reference.
    bold_mask
        Brain mask for the session-level BOLD reference.
    boldref_files
        List of BOLD reference files (same as input).
    orig2boldref_xfms
        Transforms from each run's original space to the session boldref

    """
    from niworkflows.interfaces.freesurfer import StructuralReference
    from niworkflows.interfaces.nitransforms import ConvertAffine

    workflow = Workflow(name=name)
    workflow.__desc__ = 'All BOLD runs were coregistered to create a session-level BOLD reference.'

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['boldref_files']),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=['boldref', 'bold_mask', 'orig2boldref_xfms', 'boldref_files']
        ),
        name='outputnode',
    )

    # TODO?: Do we want to denoise as well?
    # https://neurostars.org/t/ants-denoiseimage-for-fmri-epis/3091/2

    # if len(bold_runs) == 1:
    #     from niworkflows.data import load as niw_load

    #     identity_xfm = niw_load('identity_xfm')
    #     outputnode.inputs.orig2boldref_xfms = [identity_xfm]
    #     return workflow

    boldref_template = pe.Node(
        StructuralReference(
            auto_detect_sensitivity=True,
            initial_timepoint=1,
            intensity_scaling=False,
            subsample_threshold=200,
            fixed_timepoint=not unbiased,
            no_iteration=not unbiased,
            transform_outputs=True,
            out_file='boldref_template.nii.gz',
        ),
        mem_gb=2 * num_bold_runs - 1,
        name='boldref_template',
        n_procs=omp_nthreads,
    )

    to_itk = pe.MapNode(
        ConvertAffine(in_fmt='fs', out_fmt='itk'),
        iterfield=['in_xfm'],
        name='to_itk',
    )

    skullstrip_boldref_wf = init_skullstrip_bold_wf(name='skullstrip_boldref_wf')

    workflow.connect([
        (inputnode, boldref_template, [
            ('boldref_files', 'in_files'),
        ]),
        (boldref_template, outputnode, [
            ('out_file', 'boldref'),
        ]),
        (boldref_template, to_itk, [
            ('transform_outputs', 'in_xfm'),
        ]),
        (to_itk, outputnode, [
            ('out_xfm', 'orig2boldref_xfms'),
        ]),
        (inputnode, outputnode, [
            ('boldref_files', 'boldref_files'),
        ]),
        (boldref_template, skullstrip_boldref_wf, [('out_file', 'inputnode.in_file')]),
        (skullstrip_boldref_wf, outputnode, [('outputnode.mask_file', 'bold_mask')]),
    ])  # fmt:skip

    return workflow
