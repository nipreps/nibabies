import typing as ty

import nipype.interfaces.utils as niu
import nipype.pipeline.engine as pe
import templateflow.api as tf
from niworkflows.engine.workflows import LiterateWorkflow
from smriprep.interfaces.workbench import SurfaceResample
from smriprep.workflows.surfaces import (
    _collate,
    _sorted_by_basename,
    init_morph_grayords_wf,
)

from nibabies.config import DEFAULT_MEMORY_MIN_GB
from nibabies.data import load_resource
from nibabies.interfaces.utils import CiftiSelect


def init_anat_fsLR_resampling_wf(
    grayord_density: ty.Literal["91k"], mcribs: bool, name="anat_fsLR_resampling_wf"
) -> LiterateWorkflow:
    """Resample the surfaces into fsLR space"""
    workflow = LiterateWorkflow(name=name)
    fslr_density = "32k" if grayord_density == "91k" else "59k"

    workflow.__desc__ = """\
The BOLD time-series were resampled onto the left/right-symmetric template
"fsLR" [@hcppipelines].
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'subject_id',
                'subjects_dir',
                'surfaces',
                'morphometrics',
                'sphere_reg_fsLR',
            ]
        ),
        name='inputnode',
    )

    itersource = pe.Node(
        niu.IdentityInterface(fields=['hemi']),
        name='itersource',
        iterables=[('hemi', ['L', 'R'])],
    )

    outputnode = pe.Node(niu.IdentityInterface(fields=['fsLR_midthickness']), name='outputnode')

    # select white, midthickness and pial surfaces based on hemi
    select_surfaces = pe.Node(CiftiSelect(), name='select_surfaces')

    if mcribs:
        atlases = load_resource('atlases')
        # use dHCP 32k fsLR instead
        select_surfaces.inputs.template_spheres = [
            str(atlases / 'dHCP' / 'dHCP.week42.L.sphere.surf.gii'),
            str(atlases / 'dHCP' / 'dHCP.week42.R.sphere.surf.gii'),
        ]
    else:
        select_surfaces.inputs.template_spheres = [
            str(sphere)
            for sphere in tf.get(
                template='fsLR',
                density=fslr_density,
                suffix='sphere',
                space=None,
                extension='.surf.gii',
            )
        ]

    # Line 393 of FreeSurfer2CaretConvertAndRegisterNonlinear.sh
    downsampled_midthickness = pe.Node(
        SurfaceResample(method="BARYCENTRIC"),
        name="downsampled_midthickness",
        mem_gb=DEFAULT_MEMORY_MIN_GB,
    )

    joinnode = pe.JoinNode(
        niu.IdentityInterface(fields=['fsLR_midthickness']),
        name='joinnode',
        joinsource='itersource',
    )

    # resample surfaces / morphometrics to 32k
    if mcribs:
        morph_grayords_wf = init_mcribs_morph_grayords_wf(grayord_density)
        workflow.connect(
            joinnode,
            "fsLR_midthickness",
            morph_grayords_wf,
            "inputnode.fsLR_midthickness",
        )
    else:
        morph_grayords_wf = init_morph_grayords_wf(grayord_density)

    workflow.connect(
        [
            (
                inputnode,
                select_surfaces,
                [("surfaces", "surfaces"), ("sphere_reg_fsLR", "spherical_registrations")],
            ),
            (itersource, select_surfaces, [("hemi", "hemi")]),
            # Downsample midthickness to fsLR density
            (
                select_surfaces,
                downsampled_midthickness,
                [
                    ("midthickness", "surface_in"),
                    ("sphere_reg", "current_sphere"),
                    ("template_sphere", "new_sphere"),
                ],
            ),
            (downsampled_midthickness, joinnode, [("surface_out", "fsLR_midthickness")]),
            (joinnode, outputnode, [("surface_out", "fsLR_midthickness")]),
            # resample surfaces
            (
                inputnode,
                morph_grayords_wf,
                [
                    ("subject_id", "inputnode.subject_id"),
                    ("subjects_dir", "inputnode.subjects_dir"),
                ],
            ),
            (
                morph_grayords_wf,
                outputnode,
                [
                    ("outputnode.cifti_morph", "cifti_morph"),
                    ("outputnode.cifti_metadata", "cifti_metadata"),
                ],
            ),
        ]
    )
    return workflow


def init_mcribs_morph_grayords_wf(
    grayord_density: ty.Literal['91k', '170k'],
    name: str = "morph_grayords_wf",
):
    """
    Sample Grayordinates files onto the fsLR atlas.

    If `mcribs` is disabled (default), the fsaverage sphere will be resampled to fsLR.
    If `mcribs` is enabled, the M-CRIB-S sphere will be resampled to dHCP 42 week.

    Outputs are in CIFTI2 format.

    Workflow Graph
        .. workflow::
            :graph2use: colored
            :simple_form: yes

            from smriprep.workflows.surfaces import init_morph_grayords_wf
            wf = init_morph_grayords_wf(grayord_density="91k")

    Parameters
    ----------
    grayord_density : :obj:`str`
        Either `91k` or `170k`, representing the total of vertices or *grayordinates*.
    name : :obj:`str`
        Unique name for the subworkflow (default: ``"morph_grayords_wf"``)

    Inputs
    ------
    subject_id : :obj:`str`
        FreeSurfer subject ID
    subjects_dir : :obj:`str`
        FreeSurfer SUBJECTS_DIR

    Outputs
    -------
    cifti_morph : :obj:`list` of :obj:`str`
        Paths of CIFTI dscalar files
    cifti_metadata : :obj:`list` of :obj:`str`
        Paths to JSON files containing metadata corresponding to ``cifti_morph``

    """
    import templateflow.api as tf
    from nipype.interfaces.io import FreeSurferSource
    from nipype.interfaces.workbench import MetricResample
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from smriprep.interfaces.cifti import GenerateDScalar

    workflow = Workflow(name=name)
    workflow.__desc__ = f"""\
*Grayordinate* "dscalar" files [@hcppipelines] containing {grayord_density} samples were
also generated using the highest-resolution ``fsaverage`` as an intermediate standardized
surface space.
"""

    fslr_density = "32k" if grayord_density == "91k" else "59k"

    inputnode = pe.Node(
        niu.IdentityInterface(fields=["subject_id", "subjects_dir", "fsLR_midthickness"]),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=["cifti_morph", "cifti_metadata"]),
        name="outputnode",
    )

    get_surfaces = pe.Node(FreeSurferSource(), name="get_surfaces")

    surfmorph_list = pe.Node(
        niu.Merge(3, ravel_inputs=True),
        name="surfmorph_list",
        run_without_submitting=True,
    )

    # Setup Workbench command. LR ordering for hemi can be assumed, as it is imposed
    # by the iterfield of the MapNode in the surface sampling workflow above.
    resample = pe.MapNode(
        MetricResample(method="ADAP_BARY_AREA", area_metrics=True),
        name="resample",
        iterfield=[
            "in_file",
            "out_file",
            "new_sphere",
            "new_area",
            "current_sphere",
            "current_area",
        ],
    )

    atlases = load_resource('atlases')
    resample.inputs.current_sphere = [
        str(atlases / 'mcribs' / 'lh.sphere.reg.dHCP42.surf.gii'),
        str(atlases / 'mcribs' / 'rh.sphere.reg.dHCP42.surf.gii'),
    ] * 3
    # current area: FreeSurfer directory midthickness
    resample.inputs.new_sphere = [
        str(atlases / 'dHCP' / 'dHCP.week42.L.sphere.surf.gii'),
        str(atlases / 'dHCP' / 'dHCP.week42.R.sphere.surf.gii'),
    ] * 3
    # new area: dHCP midthickness

    resample.inputs.out_file = [
        f"space-fsLR_hemi-{h}_den-{grayord_density}_{morph}.shape.gii"
        # Order: curv-L, curv-R, sulc-L, sulc-R, thickness-L, thickness-R
        for morph in ('curv', 'sulc', 'thickness')
        for h in "LR"
    ]

    gen_cifti = pe.MapNode(
        GenerateDScalar(
            grayordinates=grayord_density,
        ),
        iterfield=['scalar_name', 'scalar_surfs'],
        name="gen_cifti",
    )
    gen_cifti.inputs.scalar_name = ['curv', 'sulc', 'thickness']

    # fmt: off
    workflow.connect([
        (inputnode, get_surfaces, [
            ('subject_id', 'subject_id'),
            ('subjects_dir', 'subjects_dir'),
        ]),
        (get_surfaces, surfmorph_list, [
            (('curv', _sorted_by_basename), 'in1'),
            (('sulc', _sorted_by_basename), 'in2'),
            (('thickness', _sorted_by_basename), 'in3'),
        ]),
        # (surfmorph_list, surf2surf, [('out', 'source_file')]),
        # (surf2surf, resample, [('out_file', 'in_file')]),
        (resample, gen_cifti, [
            (("out_file", _collate), "scalar_surfs")]),
        (gen_cifti, outputnode, [("out_file", "cifti_morph"),
                                 ("out_metadata", "cifti_metadata")]),
    ])
    # fmt: on

    return workflow
