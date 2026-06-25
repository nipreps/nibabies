"""Anatomical surface projections"""

from nipype.interfaces import freesurfer as fs
from nipype.interfaces import io as nio
from nipype.interfaces import utility as niu
from nipype.interfaces.ants import N4BiasFieldCorrection
from nipype.pipeline import engine as pe
from niworkflows.engine import Workflow, tag
from niworkflows.interfaces.freesurfer import (
    PatchedLTAConvert as LTAConvert,
)
from niworkflows.interfaces.freesurfer import (
    PatchedRobustRegister as RobustRegister,
)
from niworkflows.interfaces.morphology import BinaryDilation
from niworkflows.interfaces.patches import FreeSurferSource
from smriprep.interfaces.freesurfer import MakeMidthickness
from smriprep.workflows.surfaces import _extract_fs_fields

SURFACE_INPUTS = [
    't1w',
    't2w',
    'flair',
    'skullstripped_t1',
    'subjects_dir',
    'subject_id',
    # Customize aseg
    'in_aseg',
    'in_mask',
]
SURFACE_OUTPUTS = [
    'subjects_dir',
    'subject_id',
    'anat2fsnative_xfm',
    'fsnative2anat_xfm',
]


@tag('anat.recon')
def init_mcribs_surface_recon_wf(
    *,
    omp_nthreads: int,
    use_aseg: bool,
    precomputed: dict,
    mcribs_dir: str | None = None,
    name: str = 'mcribs_surface_recon_wf',
):
    """
    Reconstruct cortical surfaces using the M-CRIB-S pipeline.

    This workflow injects a precomputed segmentation into the M-CRIB-S pipeline, bypassing the
    DrawEM segmentation step that is normally performed.
    """
    from niworkflows.interfaces.nibabel import MapLabels, ReorientImage

    from nibabies.interfaces.mcribs import MCRIBReconAll

    if not use_aseg:
        raise NotImplementedError(
            'A previously computed segmentation is required for the M-CRIB-S workflow.'
        )

    inputnode = pe.Node(niu.IdentityInterface(fields=SURFACE_INPUTS), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=SURFACE_OUTPUTS), name='outputnode')

    workflow = Workflow(name=name)
    workflow.__desc__ = (
        'Brain surfaces were reconstructed with a modified `MCRIBReconAll` [M-CRIB-S, @mcribs]'
        'workflow, using the reference T2w and a pre-computed anatomical segmentation'
    )

    # mapping of labels from FS to M-CRIB-S (DrawEM)
    fs2mcribs = {
        1: 21,
        2: 51,
        3: 21,
        4: 49,
        5: 0,
        6: 17,
        7: 17,
        8: 17,
        10: 43,
        11: 41,
        12: 47,
        13: 47,
        14: 114,
        15: 0,
        16: 19,
        17: 1,
        18: 3,
        24: 83,
        26: 41,
        28: 45,
        30: 51,
        31: 49,
        41: 52,
        42: 20,
        43: 50,
        44: 0,
        46: 18,
        47: 18,
        49: 42,
        50: 40,
        51: 46,
        52: 46,
        53: 2,
        54: 4,
        58: 40,
        60: 44,
        62: 25,
        63: 50,
        77: 51,
        85: 21,
        172: 172,
        253: 48,
    }
    fs_to_mcribs = pe.Node(MapLabels(mappings=fs2mcribs), name='fs_to_mcribs')

    t2w_las = pe.Node(ReorientImage(target_orientation='LAS'), name='t2w_las')
    seg_las = pe.Node(ReorientImage(target_orientation='LAS'), name='seg_las')

    # dilated mask and use in recon-neonatal-cortex
    mask_dil = pe.Node(BinaryDilation(radius=3), name='mask_dil')
    mask_las = pe.Node(ReorientImage(target_orientation='LAS'), name='mask_las')

    # N4BiasCorrection occurs in MCRIBTissueSegMCRIBS (which is skipped)
    # Run it (with mask to rescale intensities) prior injection
    n4_mcribs = pe.Node(
        N4BiasFieldCorrection(
            dimension=3,
            bspline_fitting_distance=200,
            save_bias=True,
            copy_header=True,
            n_iterations=[50] * 5,
            convergence_threshold=1e-7,
            rescale_intensities=True,
            shrink_factor=4,
        ),
        name='n4_mcribs',
    )

    mcribs_recon = pe.Node(
        MCRIBReconAll(
            surfrecon=True,
            surfrecon_method='Deformable',
            join_thresh=1.0,
            fast_collision=True,
            nthreads=omp_nthreads,
            outdir=mcribs_dir,
        ),
        name='mcribs_recon',
        mem_gb=5,
    )
    mcribs_recon.config = {'execution': {'remove_unnecessary_outputs': False}}

    mcribs_postrecon = pe.Node(
        MCRIBReconAll(autorecon_after_surf=True, nthreads=omp_nthreads),
        name='mcribs_postrecon',
        mem_gb=5,
    )
    mcribs_postrecon.config = {'execution': {'remove_unnecessary_outputs': False}}

    fssource = pe.Node(FreeSurferSource(), name='fssource', run_without_submitting=True)
    midthickness_wf = init_make_midthickness_wf(omp_nthreads=omp_nthreads)

    workflow.connect([
        (inputnode, t2w_las, [('t2w', 'in_file')]),
        (inputnode, fs_to_mcribs, [('in_aseg', 'in_file')]),
        (inputnode, mask_dil, [('in_mask', 'in_mask')]),
        (mask_dil, mask_las, [('out_mask', 'in_file')]),
        (mask_las, mcribs_recon, [('out_file', 'mask_file')]),
        (fs_to_mcribs, seg_las, [('out_file', 'in_file')]),
        (inputnode, mcribs_recon, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id')]),
        (t2w_las, n4_mcribs, [('out_file', 'input_image')]),
        (mask_las, n4_mcribs, [('out_file', 'mask_image')]),
        (n4_mcribs, mcribs_recon, [('output_image', 't2w_file')]),
        (seg_las, mcribs_recon, [('out_file', 'segmentation_file')]),
        (inputnode, mcribs_postrecon, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id')]),
        (mcribs_recon, mcribs_postrecon, [('mcribs_dir', 'outdir')]),
        (mcribs_postrecon, fssource, [('subjects_dir', 'subjects_dir')]),
        (mcribs_postrecon, midthickness_wf, [('subjects_dir', 'inputnode.subjects_dir')]),
        (inputnode, fssource, [('subject_id', 'subject_id')]),
        (inputnode, midthickness_wf, [('subject_id', 'inputnode.subject_id')]),
        (fssource, midthickness_wf, [
            ('white', 'inputnode.white'),
            ('graymid', 'inputnode.graymid'),
        ]),
        (midthickness_wf, outputnode, [
            ('outputnode.subjects_dir', 'subjects_dir'),
            ('outputnode.subject_id', 'subject_id'),
        ]),
    ])  # fmt:skip

    if 'fsnative' not in precomputed.get('transforms', {}):
        fsnative2anat_xfm = pe.Node(
            RobustRegister(auto_sens=True, est_int_scale=True),
            name='fsnative2anat_xfm',
        )
        anat2fsnative_xfm = pe.Node(
            LTAConvert(out_lta=True, invert=True),
            name='anat2fsnative_xfm',
        )
        workflow.connect([
            (inputnode, fsnative2anat_xfm, [('t2w', 'target_file')]),
            (fssource, fsnative2anat_xfm, [('T2', 'source_file')]),
            (fsnative2anat_xfm, outputnode, [('out_reg_file', 'fsnative2anat_xfm')]),
            (fsnative2anat_xfm, anat2fsnative_xfm, [('out_reg_file', 'in_lta')]),
            (anat2fsnative_xfm, outputnode, [('out_lta', 'anat2fsnative_xfm')]),
        ])  # fmt:skip

    return workflow


@tag('anat.recon')
def init_infantfs_surface_recon_wf(
    *,
    age_months: int,
    precomputed: dict,
    omp_nthreads: int,
    use_aseg: bool = False,
    name: str = 'infantfs_surface_recon_wf',
):
    from nibabies.interfaces.freesurfer import InfantReconAll

    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=SURFACE_INPUTS), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=SURFACE_OUTPUTS), name='outputnode')

    desc = (
        'Brain surfaces were reconstructed using `infant_recon_all` [FreeSurfer '
        f'{fs.Info().looseversion() or "<ver>"}, RRID:SCR_001847, @infantfs], '
        'using the reference T1w'
    )
    desc += '.' if not use_aseg else ' and a pre-computed anatomical segmentation.'
    workflow.__desc__ = desc

    gen_recon_outdir = pe.Node(niu.Function(function=_gen_recon_dir), name='gen_recon_outdir')

    # inject the intensity-normalized skull-stripped t1w from the brain extraction workflow
    recon = pe.Node(InfantReconAll(age=age_months), name='reconall')
    if use_aseg:
        workflow.connect(inputnode, 'in_aseg', recon, 'aseg_file')

    fssource = pe.Node(FreeSurferSource(), name='fssource', run_without_submitting=True)
    midthickness_wf = init_make_midthickness_wf(omp_nthreads=omp_nthreads)

    workflow.connect([
        (inputnode, gen_recon_outdir, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id'),
        ]),
        (inputnode, recon, [
            ('skullstripped_t1', 'mask_file'),
            ('subject_id', 'subject_id'),
        ]),
        (gen_recon_outdir, recon, [
            ('out', 'outdir'),
        ]),
        (recon, fssource, [
            ('subject_id', 'subject_id'),
            (('outdir', _parent), 'subjects_dir'),
        ]),
        (recon, midthickness_wf, [
            ('subject_id', 'inputnode.subject_id'),
            (('outdir', _parent), 'inputnode.subjects_dir'),
        ]),
        (fssource, midthickness_wf, [
            ('white', 'inputnode.white'),
            ('graymid', 'inputnode.graymid'),
        ]),
        (midthickness_wf, outputnode, [
            ('outputnode.subjects_dir', 'subjects_dir'),
            ('outputnode.subject_id', 'subject_id'),
        ])
    ])  # fmt:skip

    if 'fsnative' not in precomputed.get('transforms', {}):
        fsnative2anat_xfm = pe.Node(
            RobustRegister(auto_sens=True, est_int_scale=True),
            name='fsnative2anat_xfm',
        )
        anat2fsnative_xfm = pe.Node(
            LTAConvert(out_lta=True, invert=True),
            name='anat2fsnative_xfm',
        )
        workflow.connect([
            (inputnode, fsnative2anat_xfm, [('skullstripped_t1', 'target_file')]),
            (fssource, fsnative2anat_xfm, [
                (('norm', _replace_mgz), 'source_file'),
            ]),
            (fsnative2anat_xfm, anat2fsnative_xfm, [('out_reg_file', 'in_lta')]),
            (fsnative2anat_xfm, outputnode, [
                ('out_reg_file', 'fsnative2anat_xfm'),
            ]),
            (anat2fsnative_xfm, outputnode, [
                ('out_lta', 'anat2fsnative_xfm'),
            ]),
        ])  # fmt:skip

    return workflow


@tag('anat.midthickness')
def init_make_midthickness_wf(
    *, omp_nthreads: int, name: str = 'make_midthickness_wf'
) -> pe.Workflow:
    """
    Standalone workflow to create and save cortical midthickness, derived from
    the generated white / graymid surfaces.
    """

    workflow = pe.Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['subject_id', 'subjects_dir', 'white', 'graymid']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['subject_id', 'subjects_dir']),
        name='outputnode',
    )

    midthickness = pe.MapNode(
        MakeMidthickness(thickness=True, distance=0.5, out_name='midthickness'),
        iterfield='in_file',
        name='midthickness',
        n_procs=min(omp_nthreads, 12),
    )
    save_midthickness = pe.Node(nio.DataSink(parameterization=False), name='save_midthickness')

    sync = pe.Node(
        niu.Function(
            function=_extract_fs_fields,
            output_names=['subjects_dir', 'subject_id'],
        ),
        name='sync',
    )

    workflow.connect([
        (inputnode, midthickness, [
            ('white', 'in_file'),
            ('graymid', 'graymid'),
        ]),
        (midthickness, save_midthickness, [('out_file', 'surf.@graymid')]),
        (inputnode, save_midthickness, [
            ('subjects_dir', 'base_directory'),
            ('subject_id', 'container'),
        ]),
        (save_midthickness, sync, [('out_file', 'filenames')]),
        (sync, outputnode, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id'),
        ]),
    ])  # fmt:skip
    return workflow


@tag('anat.fslr-reg')
def init_mcribs_fsLR_reg_wf(*, name: str = 'mcribs_fsLR_reg_wf') -> Workflow:
    """
    Register M-CRIB-S native surfaces to fsLR via M-CRIB-S's own deformation.

    M-CRIB-S surfaces are calibrated to M-CRIB-S's bundled ``fsaverage``/``fs_LR_32k``
    templates, and should be used to improve cortical folding registration to fsLR.

    Inputs
    ------
    sphere_reg
        Native sphere registered to the M-CRIB-S fsaverage frame (``sphere.reg2``)

    Outputs
    -------
    sphere_reg_fsLR
        Native sphere mesh with fsLR coordinates
    """
    from smriprep.interfaces.workbench import SurfaceSphereProjectUnproject

    from nibabies.data import load as load_data

    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(['sphere_reg']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(['sphere_reg_fsLR']), name='outputnode')

    project_unproject = pe.MapNode(
        SurfaceSphereProjectUnproject(),
        iterfield=['sphere_in', 'sphere_project_to', 'sphere_unproject_from'],
        name='project_unproject',
    )
    project_unproject.inputs.sphere_project_to = [
        str(load_data('mcribs/L.sphere.surf.gii')),
        str(load_data('mcribs/R.sphere.surf.gii')),
    ]
    project_unproject.inputs.sphere_unproject_from = [
        str(load_data('mcribs/L.sphere.reg.fsLR.surf.gii')),
        str(load_data('mcribs/R.sphere.reg.fsLR.surf.gii')),
    ]

    workflow.connect([
        (inputnode, project_unproject, [('sphere_reg', 'sphere_in')]),
        (project_unproject, outputnode, [('sphere_out', 'sphere_reg_fsLR')]),
    ])  # fmt:skip

    return workflow


def _parent(p):
    from pathlib import Path

    return str(Path(p).parent)


def _gen_recon_dir(subjects_dir, subject_id):
    from pathlib import Path

    p = Path(subjects_dir) / subject_id
    p.mkdir(parents=True, exist_ok=True)
    return str(p)


def _replace_mgz(in_file):
    return in_file.replace('.mgz', '.nii.gz')
