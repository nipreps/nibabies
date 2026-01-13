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
from niworkflow.data import load as load_nwf_data
from niworkflows.engine.workflows import LiterateWorkflow as Workflow


def init_bold_coreg_runs_wf(
    *,
    bold_runs: list[str],
    unbiased: bool = False,
    omp_nthreads: int = 1,
    name: str = 'bold_coreg_runs_wf',
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

    Outputs
    -------
    session_boldref
        The calculated session-level BOLD reference.
    boldref2session_xfms
        Transforms from each run's boldref to the session boldref.

    """
    from niworkflows.interfaces.freesurfer import StructuralReference

    workflow = Workflow(name=name)
    workflow.__desc__ = 'All BOLD runs were coregistered to create a session-level BOLD reference.'

    inputnode = pe.Node(niu.IdentityInterface(fields=['boldrefs']), name='inputnode')

    outputnode = pe.Node(
        niu.IdentityInterface(fields=['session_boldref', 'boldref2session_xfms']),
        name='outputnode',
    )

    # TODO?: Do we want to denoise as well?
    # https://neurostars.org/t/ants-denoiseimage-for-fmri-epis/3091/2

    boldref_iterable = pe.Node(niu.IdentityInterface(fields=['boldref']), name='boldref_iterable')
    boldref_iterable.iterables = [('boldref', bold_runs)]

    ds_session_boldref_wf = init_ds_registration_wf(
        name='ds_session_boldref_wf',
        source='rboldref',  # run-boldref
        dest='sboldref',  # session-boldref
    )

    boldref_join = pe.JoinNode(
        niu.IdentityInterface(fields=['boldref2session_xfm']),
        joinsource='boldref_iterable',
        joinfield='boldref2session_xfm',
        name='boldref_join',
    )

    # Use identity matrix for single runs
    if len(bold_runs) == 1:
        identity_xfm = load_nwf_data('itkIdentityTransform.txt')
        outputnode.inputs.boldref2session_xfms = identity_xfm

        workflow.connect([
            (inputnode, outputnode, [
                (('boldrefs', _first), 'session_boldref'),
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

    workflow.connect([
        (inputnode, boldref_template, [('boldrefs', 'in_files')]),
        (boldref_template, outputnode, [
            ('out_file', 'session_boldref'),
            ('transform_outputs', 'boldref2session_xfms'),
        ]),
    ])  # fmt:skip

    return workflow


def _first(inlist):
    return inlist[0]
