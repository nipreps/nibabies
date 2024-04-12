"""Anatomical surface projections"""


from nipype.interfaces import freesurfer as fs
from nipype.interfaces import io as nio
from nipype.interfaces import utility as niu
from nipype.pipeline import engine as pe
from niworkflows.engine.workflows import LiterateWorkflow
from niworkflows.interfaces.freesurfer import PatchedLTAConvert as LTAConvert
from niworkflows.interfaces.freesurfer import PatchedRobustRegister as RobustRegister
from smriprep.interfaces.freesurfer import MakeMidthickness
from smriprep.workflows.surfaces import _extract_fs_fields

from nibabies.config import DEFAULT_MEMORY_MIN_GB
from nibabies.data import load as load_data

SURFACE_INPUTS = [
    't1w',
    't2w',
    'flair',
    'skullstripped_t1',
    'subjects_dir',
    'subject_id',
    # Customize aseg
    'in_aseg',
]
SURFACE_OUTPUTS = [
    'subjects_dir',
    'subject_id',
    'anat2fsnative_xfm',
    'fsnative2anat_xfm',
]


def init_mcribs_surface_recon_wf(
    *,
    omp_nthreads: int,
    use_aseg: bool,
    use_mask: bool,
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

    inputnode = pe.Node(
        niu.IdentityInterface(fields=SURFACE_INPUTS + ['anat_mask']), name='inputnode'
    )
    outputnode = pe.Node(niu.IdentityInterface(fields=SURFACE_OUTPUTS), name='outputnode')

    workflow = LiterateWorkflow(name=name)
    workflow.__desc__ = (
        'Brain surfaces were reconstructed with a modified `MCRIBReconAll` [M-CRIB-S, @mcribs]'
        'workflow, using the reference T2w and a pre-computed anatomical segmentation'
    )

    # mapping of labels from FS to M-CRIB-S
    fs2mcribs = {
        2: 51,
        3: 21,
        4: 49,
        5: 0,
        7: 17,
        8: 17,
        10: 43,
        11: 41,
        12: 47,
        13: 47,
        14: 0,
        15: 0,
        16: 19,
        17: 1,
        18: 3,
        26: 41,
        28: 45,
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
        63: 50,
        253: 48,
    }
    map_labels = pe.Node(MapLabels(mappings=fs2mcribs), name='map_labels')

    t2w_las = pe.Node(ReorientImage(target_orientation='LAS'), name='t2w_las')
    seg_las = t2w_las.clone(name='seg_las')

    mcribs_recon = pe.Node(
        MCRIBReconAll(
            surfrecon=True,
            surfrecon_method='Deformable',
            join_thresh=1.0,
            fast_collision=True,
            nthreads=omp_nthreads,
        ),
        name='mcribs_recon',
        mem_gb=5,
    )
    if mcribs_dir:
        mcribs_recon.inputs.outdir = mcribs_dir
        mcribs_recon.config = {'execution': {'remove_unnecessary_outputs': False}}

    if use_mask:
        # If available, dilated mask and use in recon-neonatal-cortex
        from niworkflows.interfaces.morphology import BinaryDilation

        mask_dil = pe.Node(BinaryDilation(radius=3), name='mask_dil')
        mask_las = t2w_las.clone(name='mask_las')
        workflow.connect([
            (inputnode, mask_dil, [('anat_mask', 'in_mask')]),
            (mask_dil, mask_las, [('out_mask', 'in_file')]),
            (mask_las, mcribs_recon, [('out_file', 'mask_file')]),
        ])  # fmt:skip

    mcribs_postrecon = pe.Node(
        MCRIBReconAll(autorecon_after_surf=True, nthreads=omp_nthreads),
        name='mcribs_postrecon',
        mem_gb=5,
    )

    fssource = pe.Node(nio.FreeSurferSource(), name='fssource', run_without_submitting=True)
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
        (inputnode, t2w_las, [('t2w', 'in_file')]),
        (inputnode, map_labels, [('in_aseg', 'in_file')]),
        (map_labels, seg_las, [('out_file', 'in_file')]),
        (inputnode, mcribs_recon, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id')]),
        (t2w_las, mcribs_recon, [('out_file', 't2w_file')]),
        (seg_las, mcribs_recon, [('out_file', 'segmentation_file')]),
        (inputnode, mcribs_postrecon, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id')]),
        (mcribs_recon, mcribs_postrecon, [('mcribs_dir', 'outdir')]),
        (mcribs_postrecon, fssource, [('subjects_dir', 'subjects_dir')]),
        (inputnode, fssource, [('subject_id', 'subject_id')]),
        (fssource, midthickness, [
            ('white', 'in_file'),
            ('graymid', 'graymid'),
        ]),
        (midthickness, save_midthickness, [('out_file', 'surf.@graymid')]),
        (save_midthickness, sync, [('out_file', 'filenames')]),
        (sync, outputnode, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id'),
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


def init_mcribs_sphere_reg_wf(*, name='mcribs_sphere_reg_wf'):
    """
    Generate GIFTI registration sphere files from MCRIBS template to dHCP42 (32k).

    TODO: Clarify any distinction with fsLR
    """
    from smriprep.interfaces.surf import FixGiftiMetadata
    from smriprep.interfaces.workbench import SurfaceSphereProjectUnproject

    workflow = LiterateWorkflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(['subjects_dir', 'subject_id']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(['sphere_reg', 'sphere_reg_fsLR']),
        name='outputnode',
    )

    get_spheres = pe.Node(
        niu.Function(function=_get_dhcp_spheres),
        name='get_spheres',
        run_without_submitting=True,
    )

    # Via FreeSurfer2CaretConvertAndRegisterNonlinear.sh#L270-L273
    #
    # See https://github.com/DCAN-Labs/DCAN-HCP/tree/9291324
    sphere_gii = pe.MapNode(
        fs.MRIsConvert(out_datatype='gii'), iterfield='in_file', name='sphere_gii'
    )

    fix_meta = pe.MapNode(FixGiftiMetadata(), iterfield='in_file', name='fix_meta')

    # load template files
    atlases = load_data.cached('atlases')

    # SurfaceSphereProjectUnProject
    # project to 41k dHCP atlas sphere
    #   - sphere-in: Individual native sphere in surf directory registered to 41k atlas sphere
    #   - sphere-to: the 41k atlas sphere, in the fsaverage directory
    #   - sphere-unproject-from: 41k atlas sphere registered to dHCP 42wk sphere,
    #                            in the fsaverage directory
    #   - sphere-out: lh.sphere.reg2.dHCP42.native.surf.gii
    project_unproject = pe.MapNode(
        SurfaceSphereProjectUnproject(),
        iterfield=['sphere_in', 'sphere_project_to', 'sphere_unproject_from'],
        name='project_unproject',
    )
    project_unproject.inputs.sphere_project_to = [
        atlases / 'tpl-fsaverage_hemi-L_den-41k_desc-reg_sphere.surf.gii',
        atlases / 'tpl-fsaverage_hemi-R_den-41k_desc-reg_sphere.surf.gii',
    ]
    project_unproject.inputs.sphere_unproject_from = [
        atlases / 'tpl-dHCP_space-fsaverage_hemi-L_den-41k_desc-reg_sphere.surf.gii',
        atlases / 'tpl-dHCP_space-fsaverage_hemi-R_den-41k_desc-reg_sphere.surf.gii',
    ]

    # fmt:off
    workflow.connect([
        (inputnode, get_spheres, [
            ('subjects_dir', 'subjects_dir'),
            ('subject_id', 'subject_id'),
        ]),
        (get_spheres, sphere_gii, [(('out', _sorted_by_basename), 'in_file')]),
        (sphere_gii, fix_meta, [('converted', 'in_file')]),
        (fix_meta, project_unproject, [('out_file', 'sphere_in')]),
        (sphere_gii, outputnode, [('converted', 'sphere_reg')]),
        (project_unproject, outputnode, [('sphere_out', 'sphere_reg_fsLR')]),
    ])
    # fmt:on

    return workflow


def init_infantfs_surface_recon_wf(
    *,
    age_months: int,
    precomputed: dict,
    use_aseg: bool = False,
    name: str = 'infantfs_surface_recon_wf',
):
    from nibabies.interfaces.freesurfer import InfantReconAll

    workflow = LiterateWorkflow(name=name)
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
    fssource = pe.Node(nio.FreeSurferSource(), name='fssource', run_without_submitting=True)
    if use_aseg:
        workflow.connect(inputnode, 'in_aseg', recon, 'aseg_file')

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
        (recon, outputnode, [
            ('subject_id', 'subject_id'),
            (('outdir', _parent), 'subjects_dir'),
        ]),
        (recon, fssource, [
            ('subject_id', 'subject_id'),
            (('outdir', _parent), 'subjects_dir'),
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


def init_anat_ribbon_wf(name='anat_ribbon_wf'):
    from nipype.interfaces import fsl

    from nibabies.interfaces.workbench import CreateSignedDistanceVolume

    # 0, 1 = wm; 2, 3 = pial; 6, 7 = mid
    # note that order of lh / rh within each surf type is not guaranteed due to use
    # of unsorted glob by FreeSurferSource prior, but we can do a sort
    # to ensure consistent ordering
    workflow = pe.Workflow(name=name)

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'surfaces',  # anat_giftis,
                't1w_mask',
            ]
        ),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'anat_ribbon',
            ]
        ),
        name='outputnode',
    )

    select_wm = pe.Node(
        niu.Select(index=[0, 1]),
        name='select_wm',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    select_pial = pe.Node(
        niu.Select(index=[2, 3]),
        name='select_pial',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    select_midthick = pe.Node(
        niu.Select(index=[6, 7]),
        name='select_midthick',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    create_wm_distvol = pe.MapNode(
        CreateSignedDistanceVolume(),
        iterfield=['surf_file'],
        name='create_wm_distvol',
    )

    create_pial_distvol = pe.MapNode(
        CreateSignedDistanceVolume(),
        iterfield=['surf_file'],
        name='create_pial_distvol',
    )

    thresh_wm_distvol = pe.MapNode(
        fsl.maths.MathsCommand(args='-thr 0 -bin -mul 255'),
        iterfield=['in_file'],
        name='thresh_wm_distvol',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    uthresh_pial_distvol = pe.MapNode(
        fsl.maths.MathsCommand(args='-uthr 0 -abs -bin -mul 255'),
        iterfield=['in_file'],
        name='uthresh_pial_distvol',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    bin_wm_distvol = pe.MapNode(
        fsl.maths.UnaryMaths(operation='bin'),
        iterfield=['in_file'],
        name='bin_wm_distvol',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    bin_pial_distvol = pe.MapNode(
        fsl.maths.UnaryMaths(operation='bin'),
        iterfield=['in_file'],
        name='bin_pial_distvol',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    split_wm_distvol = pe.Node(
        niu.Split(splits=[1, 1]),
        name='split_wm_distvol',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    merge_wm_distvol_no_flatten = pe.Node(
        niu.Merge(2),
        no_flatten=True,
        name='merge_wm_distvol_no_flatten',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    make_ribbon_vol = pe.MapNode(
        fsl.maths.MultiImageMaths(op_string='-mas %s -mul 255'),
        iterfield=['in_file', 'operand_files'],
        name='make_ribbon_vol',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    bin_ribbon_vol = pe.MapNode(
        fsl.maths.UnaryMaths(operation='bin'),
        iterfield=['in_file'],
        name='bin_ribbon_vol',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    split_squeeze_ribbon_vol = pe.Node(
        niu.Split(splits=[1, 1], squeeze=True),
        name='split_squeeze_ribbon',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
        run_without_submitting=True,
    )

    combine_ribbon_vol_hemis = pe.Node(
        fsl.maths.BinaryMaths(operation='add'),
        name='combine_ribbon_vol_hemis',
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    # make HCP-style ribbon volume in T1w space
    workflow.connect(
        [
            (inputnode, select_wm, [('surfaces', 'inlist')]),
            (inputnode, select_pial, [('surfaces', 'inlist')]),
            (inputnode, select_midthick, [('surfaces', 'inlist')]),
            (select_wm, create_wm_distvol, [(('out', _sorted_by_basename), 'surf_file')]),
            (inputnode, create_wm_distvol, [('t1w_mask', 'ref_file')]),
            (select_pial, create_pial_distvol, [(('out', _sorted_by_basename), 'surf_file')]),
            (inputnode, create_pial_distvol, [('t1w_mask', 'ref_file')]),
            (create_wm_distvol, thresh_wm_distvol, [('out_file', 'in_file')]),
            (create_pial_distvol, uthresh_pial_distvol, [('out_file', 'in_file')]),
            (thresh_wm_distvol, bin_wm_distvol, [('out_file', 'in_file')]),
            (uthresh_pial_distvol, bin_pial_distvol, [('out_file', 'in_file')]),
            (bin_wm_distvol, split_wm_distvol, [('out_file', 'inlist')]),
            (split_wm_distvol, merge_wm_distvol_no_flatten, [('out1', 'in1')]),
            (split_wm_distvol, merge_wm_distvol_no_flatten, [('out2', 'in2')]),
            (bin_pial_distvol, make_ribbon_vol, [('out_file', 'in_file')]),
            (merge_wm_distvol_no_flatten, make_ribbon_vol, [('out', 'operand_files')]),
            (make_ribbon_vol, bin_ribbon_vol, [('out_file', 'in_file')]),
            (bin_ribbon_vol, split_squeeze_ribbon_vol, [('out_file', 'inlist')]),
            (split_squeeze_ribbon_vol, combine_ribbon_vol_hemis, [('out1', 'in_file')]),
            (split_squeeze_ribbon_vol, combine_ribbon_vol_hemis, [('out2', 'operand_file')]),
            (combine_ribbon_vol_hemis, outputnode, [('out_file', 'anat_ribbon')]),
        ]
    )
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


def _sorted_by_basename(inlist):
    from os.path import basename

    return sorted(inlist, key=lambda x: str(basename(x)))


def _get_dhcp_spheres(subject_id: str, subjects_dir: str) -> list:
    from pathlib import Path

    out = []
    for hemi in 'lr':
        sphere = Path(subjects_dir) / subject_id / 'surf' / f'{hemi}h.sphere.reg2'
        if not sphere.exists():
            raise OSError('MCRIBS spherical registration not found.')
        out.append(str(sphere))
    return out
