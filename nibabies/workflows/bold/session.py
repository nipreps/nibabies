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
    bold_runs: list[str],
    unbiased: bool = False,
    omp_nthreads: int = 1,
    name: str = 'coreg_bolds_wf',
) -> pe.Workflow:
    """
    Register all BOLD runs of a session to a single session-reference.

    Parameters
    ----------
    bold_runs : :obj:`list` of :obj:`str`
        List of paths to BOLD files (used for names/metadata).
    omp_nthreads : :obj:`int`
        Number of threads.

    Inputs
    ------
    boldrefs
        List of coregistered boldref images (one per run).
    bold_masks
        List of brain masks corresponding to each boldref. Unused if only one run.

    Outputs
    -------
    boldref
        The calculated session-level BOLD reference.
    bold_mask
        The brain mask for the session-level BOLD reference.
    run2session_xfms
        Transforms from each run's boldref to the session boldref.

    """
    from niworkflows.interfaces.freesurfer import StructuralReference

    workflow = Workflow(name=name)
    workflow.__desc__ = 'All BOLD runs were coregistered to create a session-level BOLD reference.'

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['boldrefs', 'bold_masks'],
        ),
        name='inputnode',
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['boldref', 'bold_mask', 'run2session_xfms']),
        name='outputnode',
    )

    # TODO?: Do we want to denoise as well?
    # https://neurostars.org/t/ants-denoiseimage-for-fmri-epis/3091/2

    if len(bold_runs) == 1:
        from niworkflows.data import load as niw_load

        identity_xfm = niw_load('identity_xfm')
        outputnode.inputs.run2session_xfms = [identity_xfm]
        workflow.connect([
            (inputnode, outputnode, [
                (('boldrefs', _first), 'boldref'),
                (('bold_masks', _first), 'bold_mask'),
            ]),
        ])  # fmt:skip
        return workflow

    boldref_template = pe.Node(
        StructuralReference(
            auto_detect_sensitivity=True,
            initial_timepoint=1,
            intensity_scaling=False,
            subsample_threshold=200,
            fixed_timepoint=not unbiased,
            no_iteration=not unbiased,
            transform_outputs=True,
        ),
        mem_gb=2 * len(bold_runs) - 1,
        name='boldref_template',
        n_procs=omp_nthreads,
    )

    skullstrip_sboldref_wf = init_skullstrip_bold_wf(name='skullstrip_sboldref_wf')

    workflow.connect([
        (inputnode, boldref_template, [('boldrefs', 'in_files')]),
        (boldref_template, outputnode, [
            ('out_file', 'session_boldref'),
            ('transform_outputs', 'boldref2session_xfms'),
        ]),
        (boldref_template, skullstrip_sboldref_wf, [('out_file', 'inputnode.bold_file')]),
        (skullstrip_sboldref_wf, outputnode, [('outputnode.bold_mask', 'bold_mask')]),
    ])  # fmt:skip

    return workflow


def _first(inlist):
    return inlist[0]
