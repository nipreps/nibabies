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
# We support and encourage derived works from this project, please read
# about our expectations at
#
#     https://www.nipreps.org/community/licensing/
#
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow as Workflow
from niworkflows.interfaces.header import ValidateImage
from niworkflows.utils.misc import pass_dummy_scans
from nibabies.interfaces.reference import DetectReferenceFrame

DEFAULT_MEMORY_MIN_GB = 0.01


def init_raw_boldref_wf(
    bold_file: str | None = None,
    multiecho: bool = False,
    ref_frame_start: int = 16,
    name: str = 'raw_boldref_wf',
):
    """
    Build a workflow that generates reference BOLD images for a series.

    The raw reference image is the target of :abbr:`HMC (head motion correction)`, and a
    contrast-enhanced reference is the subject of distortion correction, as well as
    boundary-based registration to anat and template spaces.

    This workflow assumes only one BOLD file has been passed.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.bold.reference import init_raw_boldref_wf
            wf = init_raw_boldref_wf()

    Parameters
    ----------
    bold_file : :obj:`str`
        BOLD series NIfTI file
    multiecho : :obj:`bool`
        If multiecho data was supplied, data from the first echo will be selected
    ref_frame_start: :obj:`int`
        BOLD frame to start creating the reference map from.
    estimate_good_refframe: :obj:`bool`
        If True, use a heuristic to find a single low-motion BOLD reference frame out of each timeseries
        instead of running RobustAverage over all frames after ref_frame_start.
    name : :obj:`str`
        Name of workflow (default: ``raw_boldref_wf``)

    Inputs
    ------
    bold_file : str
        BOLD series NIfTI file
    dummy_scans : int or None
        Number of non-steady-state volumes specified by user at beginning of ``bold_file``

    Outputs
    -------
    bold_file : str
        Validated BOLD series NIfTI file
    boldref : str
        Reference image to which BOLD series is motion corrected
    skip_vols : int
        Number of non-steady-state volumes selected at beginning of ``bold_file``
    algo_dummy_scans : int
        Number of non-steady-state volumes agorithmically detected at
        beginning of ``bold_file``

    """
    from niworkflows.interfaces.images import RobustAverage

    workflow = Workflow(name=name)
    workflow.__desc__ = f"""\
First, a reference volume was generated{' from the shortest echo of the BOLD run' * multiecho},
using a custom methodology of *NiBabies*, for use in head motion correction.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['bold_file', 'dummy_scans']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'bold_file',
                'boldref',
                'skip_vols',
                'algo_dummy_scans',
                'validation_report',
            ]
        ),
        name='outputnode',
    )

    # Simplify manually setting input image
    if bold_file is not None:
        inputnode.inputs.bold_file = bold_file

    validation_and_dummies_wf = init_validation_and_dummies_wf()

    workflow.connect([
        (inputnode, validation_and_dummies_wf, [
            ('bold_file', 'inputnode.bold_file'),
            ('dummy_scans', 'inputnode.dummy_scans'),
        ]),
        (validation_and_dummies_wf, outputnode, [
            ('outputnode.bold_file', 'bold_file'),
            ('outputnode.skip_vols', 'skip_vols'),
            ('outputnode.algo_dummy_scans', 'algo_dummy_scans'),
            ('outputnode.validation_report', 'validation_report'),
        ])
    ])  # fmt:skip
    # Drop frames to avoid startle when MRI begins acquiring
    if not estimate_good_refframe:
        select_frames = pe.Node(
            niu.Function(function=_select_frames, output_names=['start_frame', 't_mask']),
            name='select_frames',
        )
        select_frames.inputs.ref_frame_start = ref_frame_start

        gen_avg = pe.Node(RobustAverage(), name='gen_avg', mem_gb=1)
        workflow.connect([
            (validation_and_dummies_wf, gen_avg, [
                ('outputnode.bold_file', 'in_file'),
            ]),
            (validation_and_dummies_wf, select_frames, [
                ('outputnode.bold_file', 'in_file'),
            ]),
            (inputnode, select_frames, [('dummy_scans', 'dummy_scans')]),
            (select_frames, gen_avg, [('t_mask', 't_mask')]),
            (gen_avg, outputnode, [('out_file', 'boldref')]),
        ])  # fmt:skip
    else:  # Select a single low-motion frame
        detect_referenece_frame = pe.Node(DetectReferenceFrame(), name='detect_referenece_frame')
        detect_referenece_frame.inputs.ref_frame_start = ref_frame_start
        workflow.connect([
            (validation_and_dummies_wf, detect_referenece_frame, [
                ('outputnode.bold_file', 'in_file'),
            ]),
            (inputnode, detect_referenece_frame, [('dummy_scans', 'dummy_scans')]),
            (detect_referenece_frame, outputnode, [('out_file', 'boldref')]),
        ])  # fmt:skip
    return workflow


def _select_frames(
    in_file: str, ref_frame_start: int, dummy_scans: int | None
) -> tuple[int, list]:
    import warnings

    import nibabel as nb
    import numpy as np

    img = nb.load(in_file)
    img_len = img.shape[3]

    # Ensure start index is the largest of the two
    # Will usually be `ref_frame_start`
    start_frame = max(ref_frame_start, dummy_scans) if dummy_scans else ref_frame_start

    if start_frame >= img_len:
        warnings.warn(
            f'Caculating the BOLD reference starting on frame {start_frame} but only {img_len} '
            'volumes in BOLD file, so using last volume.',
            stacklevel=1,
        )
        start_frame = img_len - 1

    t_mask = np.array([False] * img_len, dtype=bool)
    t_mask[start_frame:] = True
    return start_frame, list(t_mask)





def init_validation_and_dummies_wf(
    bold_file: str | None = None,
    name: str = 'validation_and_dummies_wf',
):
    """
    Build a workflow that validates a BOLD image and detects non-steady-state volumes.

    Workflow Graph
        .. workflow::
            :graph2use: orig
            :simple_form: yes

            from fmriprep.workflows.bold.reference import init_validation_and_dummies_wf
            wf = init_validation_and_dummies_wf()

    Parameters
    ----------
    bold_file : :obj:`str`
        BOLD series NIfTI file
    name : :obj:`str`
        Name of workflow (default: ``validation_and_dummies_wf``)

    Inputs
    ------
    bold_file : str
        BOLD series NIfTI file
    dummy_scans : int or None
        Number of non-steady-state volumes specified by user at beginning of ``bold_file``

    Outputs
    -------
    bold_file : str
        Validated BOLD series NIfTI file
    skip_vols : int
        Number of non-steady-state volumes selected at beginning of ``bold_file``
    algo_dummy_scans : int
        Number of non-steady-state volumes agorithmically detected at
        beginning of ``bold_file``

    """
    from niworkflows.interfaces.bold import NonsteadyStatesDetector

    workflow = Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['bold_file', 'dummy_scans']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'bold_file',
                'skip_vols',
                'algo_dummy_scans',
                't_mask',
                'validation_report',
            ]
        ),
        name='outputnode',
    )

    # Simplify manually setting input image
    if bold_file is not None:
        inputnode.inputs.bold_file = bold_file

    val_bold = pe.Node(
        ValidateImage(),
        name='val_bold',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    get_dummy = pe.Node(NonsteadyStatesDetector(), name='get_dummy')

    calc_dummy_scans = pe.Node(
        niu.Function(function=pass_dummy_scans, output_names=['skip_vols_num']),
        name='calc_dummy_scans',
        run_without_submitting=True,
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    workflow.connect([
        (inputnode, val_bold, [('bold_file', 'in_file')]),
        (val_bold, outputnode, [
            ('out_file', 'bold_file'),
            ('out_report', 'validation_report'),
        ]),
        (inputnode, get_dummy, [('bold_file', 'in_file')]),
        (inputnode, calc_dummy_scans, [('dummy_scans', 'dummy_scans')]),
        (get_dummy, calc_dummy_scans, [('n_dummy', 'algo_dummy_scans')]),
        (get_dummy, outputnode, [
            ('n_dummy', 'algo_dummy_scans'),
            ('t_mask', 't_mask'),
        ]),
        (calc_dummy_scans, outputnode, [('skip_vols_num', 'skip_vols')]),
    ])  # fmt:skip

    return workflow
