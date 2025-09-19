# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Within-baby registration of a T1w into a T2w image."""

from __future__ import annotations

from collections import defaultdict

from nipype.interfaces import (
    utility as niu,
)
from nipype.interfaces.ants.base import Info as ANTsInfo
from nipype.pipeline import engine as pe
from niworkflows.engine import Workflow, tag
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from smriprep.workflows.fit.registration import (
    TemplateDesc,
    _fmt_cohort,
    get_metadata,
    tf_ver,
)


@tag('anat.coreg')
def init_coregistration_wf(
    *,
    bspline_fitting_distance: int = 200,
    mem_gb: float = 3.0,
    omp_nthreads: int | None = None,
    sloppy: bool = False,
    debug: bool = False,
    t1w_mask: bool = False,
    probmap: bool = True,
    name: str = 'coregistration_wf',
):
    """
    Set-up a T2w-to-T1w within-baby co-registration framework.

    See the ANTs' registration config file (under ``nibabies/data``) for further
    details.
    The main surprise in it is that, for some participants, accurate registration
    requires extra degrees of freedom (one affine level and one SyN level) to ensure
    that the T1w and T2w images align well.
    I attribute this requirement to the following potential reasons:

      * The T1w image and the T2w image were acquired in different sessions, apart in
        time enough for growth to happen.
        Although this is, in theory possible, it doesn't seem the images we have tested
        on are acquired on different sessions.
      * The skull is still so malleable that a change of position of the baby inside the
        coil made an actual change on the overall shape of their head.
      * Nonlinear distortions of the T1w and T2w images are, for some reason, more notorious
        for babies than they are for adults.
        We would need to look into each sequence's details to confirm this.

    Parameters
    ----------
    bspline_fitting_distance : :obj:`float`
        Distance in mm between B-Spline control points for N4 INU estimation.
    mem_gb : :obj:`float`
        Base memory fingerprint unit.
    name : :obj:`str`
        This particular workflow's unique name (Nipype requirement).
    omp_nthreads : :obj:`int`
        The number of threads for individual processes in this workflow.
    sloppy : :obj:`bool`
        Run in *sloppy* mode.
    debug : :obj:`bool`
        Produce intermediate registration files
    t1w_mask : :obj:`bool`
        A precomputed mask for the T1w is available. In this case, generate a
        quick mask to assist in coregistration, but use the precomputed mask
        as the final output.
    probmap: :obj:`bool`
        A probabilistic brainmask is present in T2w space.

    Inputs
    ------
    in_t1w : :obj:`str`
        The preprocessed input T1w image (Denoising/INU/Clipping)
    in_t2w : :obj:`str`
        The preprocessed input T2w image (Denoising/INU/Clipping)
    in_mask : :obj:`str`
        The brainmask.
        If `t1w_mask` is False, will be in T2w space.
        If `t1w_mask` is True, will be in T1w space.
    in_probmap : :obj:`str`
        The probabilistic brainmask, as obtained in T2w space.

    Outputs
    -------
    t1w_coreg : :obj:`str`
        The preprocessed T1w image (INU and clipping), in its native space.
    t2w_coreg : :obj:`str`
        The preprocessed T2w image (INU and clipping), aligned into the T1w's space.
    t1w_brain : :obj:`str`
        The preprocessed, brain-extracted T1w image.
    t1w_mask : :obj:`str`
        The binary brainmask in T1w space.
    t1w2t2w_xfm : :obj:`str`
        The T1w-to-T2w mapping.

    """
    from nipype.interfaces.ants import N4BiasFieldCorrection
    from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
    from niworkflows.interfaces.fixes import FixHeaderRegistration as Registration
    from niworkflows.interfaces.nibabel import ApplyMask, Binarize, BinaryDilation

    from nibabies.utils.misc import get_file

    workflow = pe.Workflow(name)

    inputnode = pe.Node(
        niu.IdentityInterface(fields=['in_t1w', 'in_t2w', 'in_mask', 'in_probmap']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                't1w_preproc',
                't1w_brain',
                't1w_mask',
                't1w2t2w_xfm',
                't2w2t1w_xfm',
                't2w_preproc',
            ]
        ),
        name='outputnode',
    )

    # Dilate t2w mask for easier t1->t2 registration
    fixed_masks_arg = pe.Node(niu.Merge(3), name='fixed_masks_arg', run_without_submitting=True)
    reg_mask = pe.Node(BinaryDilation(radius=8, iterations=3), name='reg_mask')
    refine_mask = pe.Node(BinaryDilation(radius=8, iterations=1), name='refine_mask')

    # Set up T1w -> T2w within-subject registration
    coreg = pe.Node(
        Registration(from_file=get_file('nibabies', 'data/t1-t2-coreg.json')),
        name='coreg',
        n_procs=omp_nthreads,
        mem_gb=mem_gb,
    )
    coreg.inputs.float = sloppy
    if debug:
        coreg.inputs.args = '--write-interval-volumes 5'
        coreg.inputs.output_inverse_warped_image = sloppy
        coreg.inputs.output_warped_image = sloppy

    final_n4 = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            bspline_fitting_distance=bspline_fitting_distance,
            save_bias=True,
            copy_header=True,
            n_iterations=[50] * 5,
            convergence_threshold=1e-7,
            rescale_intensities=True,
            shrink_factor=4,
        ),
        n_procs=omp_nthreads,
        name='final_n4',
    )
    # Move the T2w into T1w space, and apply the mask to the T1w
    map_t2w = pe.Node(ApplyTransforms(interpolation='BSpline'), name='map_t2w', mem_gb=1)
    apply_mask = pe.Node(ApplyMask(), name='apply_mask')

    # fmt: off
    workflow.connect([
        (inputnode, final_n4, [('in_t1w', 'input_image')]),
        (inputnode, coreg, [('in_t1w', 'moving_image'),
                            ('in_t2w', 'fixed_image')]),
        (reg_mask, fixed_masks_arg, [
            ('out_file', 'in1'),
            ('out_file', 'in2')]),
        (refine_mask, fixed_masks_arg, [('out_file', 'in3')]),
        (inputnode, map_t2w, [
            ('in_t1w', 'reference_image'),
            ('in_t2w', 'input_image')]),
        (fixed_masks_arg, coreg, [('out', 'fixed_image_masks')]),
        (coreg, map_t2w, [
            ('inverse_composite_transform', 'transforms'),
        ]),
        (final_n4, apply_mask, [('output_image', 'in_file')]),
        (final_n4, outputnode, [('output_image', 't1w_preproc')]),
        (map_t2w, outputnode, [('output_image', 't2w_preproc')]),
        (apply_mask, outputnode, [('out_file', 't1w_brain')]),
        (coreg, outputnode, [('composite_transform', 't1w2t2w_xfm')]),
        (coreg, outputnode, [('inverse_composite_transform', 't2w2t1w_xfm')]),
    ])
    # fmt: on

    if t1w_mask:
        # The input mask is already in T1w space.
        # Generate a quick, rough mask of the T2w to be used to facilitate co-registration.
        from sdcflows.interfaces.brainmask import BrainExtraction

        masker = pe.Node(BrainExtraction(), name='t2w_masker')
        # fmt:off
        workflow.connect([
            (inputnode, masker, [('in_t2w', 'in_file')]),
            (masker, reg_mask, [('out_mask', 'in_file')]),
            (masker, refine_mask, [('out_mask', 'in_file')]),
            (inputnode, apply_mask, [('in_mask', 'in_mask')]),
            (inputnode, outputnode, [('in_mask', 't1w_mask')]),
        ])
        # fmt:on
        return workflow

    if probmap:
        # The T2w mask from the brain extraction workflow will be mapped to T1w space
        map_mask = pe.Node(ApplyTransforms(interpolation='Gaussian'), name='map_mask', mem_gb=1)
        thr_mask = pe.Node(Binarize(thresh_low=0.80), name='thr_mask')
        # fmt:off
        workflow.connect([
            (inputnode, reg_mask, [('in_mask', 'in_file')]),
            (inputnode, refine_mask, [('in_mask', 'in_file')]),
            (inputnode, map_mask, [
                ('in_t1w', 'reference_image'),
                ('in_probmap', 'input_image')]),
            (coreg, map_mask, [
                ('inverse_composite_transform', 'transforms'),
            ]),
            (map_mask, thr_mask, [('output_image', 'in_file')]),
            (map_mask, final_n4, [('output_image', 'weight_image')]),
            (thr_mask, outputnode, [('out_mask', 't1w_mask')]),
            (thr_mask, apply_mask, [('out_mask', 'in_mask')]),
        ])
        # fmt:on
        return workflow

    # A precomputed T2w mask was provided
    map_precomp_mask = pe.Node(
        ApplyTransforms(interpolation='MultiLabel'), name='map_precomp_mask'
    )
    # fmt:off
    workflow.connect([
        (inputnode, reg_mask, [('in_mask', 'in_file')]),
        (inputnode, refine_mask, [('in_mask', 'in_file')]),
        (inputnode, map_precomp_mask, [
            ('in_t1w', 'reference_image'),
            ('in_mask', 'input_image')]),
        (coreg, map_precomp_mask, [
            ('inverse_composite_transform', 'transforms'),
        ]),
        (map_precomp_mask, final_n4, [('output_image', 'weight_image')]),
        (map_precomp_mask, outputnode, [('output_image', 't1w_mask')]),
        (map_precomp_mask, apply_mask, [('output_image', 'in_mask')]),
    ])
    # fmt:on
    return workflow


@tag('anat.coreg-derivatives')
def init_coregister_derivatives_wf(
    *, t1w_mask: bool, t1w_aseg: bool, t2w_aseg: bool, name: str = 'coregister_derivatives_wf'
):
    """Move derivatives from T1w / T2w space."""
    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=['t1w_ref', 't2w_ref', 't1w2t2w_xfm', 't1w_mask', 't1w_aseg', 't2w_aseg']
        ),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['t2w_mask', 't1w_aseg', 't2w_aseg']), name='outputnode'
    )

    if t1w_mask:
        t1wmask2t2w = pe.Node(ApplyTransforms(interpolation='MultiLabel'), name='t1wmask2t2w')
        # fmt:off
        workflow.connect([
            (inputnode, t1wmask2t2w, [
                ('t1w_mask', 'input_image'),
                ('t1w2t2w_xfm', 'transforms'),
                ('t2w_ref', 'reference_image')]),
            (t1wmask2t2w, outputnode, [('output_image', 't2w_mask')])
        ])
        # fmt:on
    if t1w_aseg:
        # fmt:off
        t1waseg2t2w = pe.Node(ApplyTransforms(interpolation='MultiLabel'), name='t1waseg2t2w')
        workflow.connect([
            (inputnode, t1waseg2t2w, [
                ('t1w_aseg', 'input_image'),
                ('t1w2t2w_xfm', 'transforms'),
                ('t2w_ref', 'reference_image')]),
            (t1waseg2t2w, outputnode, [('output_image', 't2w_aseg')])
        ])
        # fmt:on
    if t2w_aseg:
        # fmt:off
        t2waseg2t1w = pe.Node(ApplyTransforms(interpolation='MultiLabel'), name='t2waseg2t1w')
        t2waseg2t1w.inputs.invert_transform_flags = [True, False]
        workflow.connect([
            (inputnode, t2waseg2t1w, [
                ('t2w_aseg', 'input_image'),
                ('t1w2t2w_xfm', 'transforms'),
                ('t1w_ref', 'reference_image')]),
            (t2waseg2t1w, outputnode, [('output_image', 't1w_aseg')])
        ])
        # fmt:on
    return workflow


@tag('anat.concat-reg')
def init_concat_registrations_wf(
    *,
    templates,
    name='concat_registrations_wf',
):
    """
    Concatenate two transforms to produce a single composite transform from native to template.

    Parameters
    ----------
    templates : :obj:`list` of :obj:`str`
        List of standard space fullnames (e.g., ``MNI152NLin6Asym``
        or ``MNIPediatricAsym:cohort-4``) which are targets for spatial
        normalization.

    Inputs
    ------
    anat_preproc
        The anatomical reference, to be used as a reference image for std2anat_xfm
    intermediate
        Standard space fullname (usually ``MNIInfant:cohort-X``) which serves as
        the intermediate space between native and *template*
    anat2std_xfm
        The incoming anat2std transform (from native to MNIInfant)
    std2anat_xfm
        The incoming std2anat transform (from MNIInfant to native)

    Outputs
    -------
    anat2std_xfm
        Anatomical -> MNIInfant -> template transform.
    std2anat_xfm
        The template -> MNIInfant -> anatomical transform.
    template
        Template name extracted from the input parameter ``template``, for further
        use in downstream nodes.
    template_spec
        Template specifications extracted from the input parameter ``template``, for
        further use in downstream nodes.

    """
    from nibabies.interfaces.patches import CompositeTransformUtil

    ntpls = len(templates)
    workflow = Workflow(name=name)

    if templates:
        workflow.__desc__ = """\
Volume-based spatial normalization to {targets} ({targets_id}) was performed by
concatenating two registrations with `antsRegistration` (ANTs {ants_ver}). First, the
anatomical reference was registered to the Infant MNI templates (@mniinfant). Separately,
the Infant MNI templates were registered to {targets}, with the saved transform to template
stored for reuse and accessed with *TemplateFlow* [{tf_ver}, @templateflow]:
""".format(
            ants_ver=ANTsInfo.version() or '(version unknown)',
            targets='{} standard space{}'.format(
                defaultdict('several'.format, {1: 'one', 2: 'two', 3: 'three', 4: 'four'})[ntpls],
                's' * (ntpls != 1),
            ),
            targets_id=', '.join(templates),
            tf_ver=tf_ver,
        )

        # Append template citations to description
        for template in templates:
            template_meta = get_metadata(template.split(':')[0])
            template_refs = ['@{}'.format(template.split(':')[0].lower())]

            if template_meta.get('RRID', None):
                template_refs += [f'RRID:{template_meta["RRID"]}']

            workflow.__desc__ += """\
*{template_name}* [{template_refs}; TemplateFlow ID: {template}]""".format(
                template=template,
                template_name=template_meta['Name'],
                template_refs=', '.join(template_refs),
            )
            workflow.__desc__ += '.\n' if template == templates[-1] else ', '

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'template',  # template identifier (name[+cohort])
                'intermediate',  # intermediate space (name[+cohort])
                'anat2int_xfm',  # anatomical -> intermediate
                'int2anat_xfm',  # intermediate -> anatomical
            ]
        ),
        name='inputnode',
    )
    inputnode.inputs.template = templates

    out_fields = [
        'anat2std_xfm',
        'std2anat_xfm',
        'template',
        'template_spec',
    ]
    outputnode = pe.Node(niu.IdentityInterface(fields=out_fields), name='outputnode')

    fetch_tf_xfms = pe.MapNode(
        niu.Function(function=_get_intermediate_xfms, output_names=['int2tgt_xfm', 'tgt2int_xfm']),
        name='fetch_tf_xfms',
        iterfield=['target'],
        run_without_submitting=True,
    )

    split_desc = pe.MapNode(
        TemplateDesc(), run_without_submitting=True, iterfield='template', name='split_desc'
    )

    fmt_cohort = pe.MapNode(
        niu.Function(function=_fmt_cohort, output_names=['template', 'spec']),
        name='fmt_cohort',
        run_without_submitting=True,
        iterfield=['template', 'spec'],
    )

    # Disassemble each composite transform individually for readability
    dis_anat2int = pe.Node(
        CompositeTransformUtil(process='disassemble', output_prefix='anat2int'),
        name='dis_anat2int',
    )

    dis_int2std = pe.Node(
        CompositeTransformUtil(process='disassemble', output_prefix='int2std'),
        name='dis_int2std',
    )

    dis_std2int = pe.Node(
        CompositeTransformUtil(process='disassemble', output_prefix='std2int'),
        name='dis_std2int',
    )

    dis_int2anat = pe.Node(
        CompositeTransformUtil(process='disassemble', output_prefix='int2anat'),
        name='dis_int2anat',
    )

    order_anat2std = pe.Node(niu.Merge(4), name='order_anat2std', run_without_submitting=True)
    order_std2anat = pe.Node(niu.Merge(4), name='order_std2anat', run_without_submitting=True)

    assemble_anat2std = pe.Node(
        CompositeTransformUtil(process='assemble', out_file='anat2std.h5'),
        name='assemble_anat2std',
    )
    assemble_std2anat = pe.Node(
        CompositeTransformUtil(process='assemble', out_file='std2anat.h5'),
        name='assemble_std2anat',
    )

    workflow.connect([
        # Transform concatenation
        (inputnode, dis_anat2int, [('anat2int_xfm', 'in_file')]),
        (inputnode, dis_int2anat, [('int2anat_xfm', 'in_file')]),
        (inputnode, fetch_tf_xfms, [('intermediate', 'intermediate'),
                                    ('template', 'target')]),
        (fetch_tf_xfms, dis_int2std, [('int2tgt_xfm', 'in_file')]),
        (fetch_tf_xfms, dis_std2int, [('tgt2int_xfm', 'in_file')]),
        (dis_anat2int, order_anat2std, [
            ('affine_transform', 'in1'),
            ('displacement_field', 'in2'),
        ]),
        (dis_int2std, order_anat2std, [
            ('affine_transform', 'in3'),
            ('displacement_field', 'in4'),
        ]),
        # Because std2anat are inverse transforms, warp is first
        (dis_std2int, order_std2anat, [
            ('affine_transform', 'in2'),
            ('displacement_field', 'in1'),
        ]),
        (dis_int2anat, order_std2anat, [
            ('affine_transform', 'in4'),
            ('displacement_field', 'in3'),
        ]),
        (order_anat2std, assemble_anat2std, [('out', 'in_file')]),
        (order_std2anat, assemble_std2anat, [('out', 'in_file')]),
        (assemble_anat2std, outputnode, [('out_file', 'anat2std_xfm')]),
        (assemble_std2anat, outputnode, [('out_file', 'std2anat_xfm')]),

        # Template name wrangling
        (inputnode, split_desc, [('template', 'template')]),
        (split_desc, fmt_cohort, [
            ('name', 'template'),
            ('spec', 'spec'),
        ]),
        (fmt_cohort, outputnode, [
            ('template', 'template'),
            ('spec', 'template_spec'),
        ]),
    ])  # fmt:skip

    return workflow


def _get_intermediate_xfms(intermediate, target):
    import templateflow.api as tf

    # Native -> MNIInfant:cohort-X (int) -> Target (std)
    ispace, _, icohort = intermediate.partition(':cohort-')
    ispaceid = f'{ispace}+{icohort}' if icohort else ispace

    int2std = tf.get(
        target,
        suffix='xfm',
        **{'from': ispaceid},
        raise_empty=True,
    )

    std2int = tf.get(
        ispace,
        cohort=icohort or None,
        suffix='xfm',
        **{'from': target},
        raise_empty=True,
    )

    return int2std, std2int
