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

    outputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                'fsLR_midthickness',
                'cifti_morph',
                'cifti_metadata',
            ]
        ),
        name='outputnode',
    )

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

    # fmt:off
    workflow.connect([
        (inputnode, select_surfaces, [
            ("surfaces", "surfaces"),
            ("sphere_reg_fsLR", "spherical_registrations")]),
        (itersource, select_surfaces, [("hemi", "hemi")]),
        # Downsample midthickness to fsLR density
        (select_surfaces, downsampled_midthickness, [
            ("midthickness", "surface_in"),
            ("sphere_reg", "current_sphere"),
            ("template_sphere", "new_sphere")]),
        (downsampled_midthickness, joinnode, [("surface_out", "fsLR_midthickness")]),
        (joinnode, outputnode, [("surface_out", "fsLR_midthickness")]),
        # resample morphometrics to fsLR 32k
        (inputnode, morph_grayords_wf, [
            ("subject_id", "inputnode.subject_id"),
            ("subjects_dir", "inputnode.subjects_dir")]),
        (morph_grayords_wf, outputnode, [
            ("outputnode.cifti_morph", "cifti_morph"),
            ("outputnode.cifti_metadata", "cifti_metadata")]),
    ])
    # fmt:on
    return workflow


def init_mcribs_morph_grayords_wf(
    grayord_density: ty.Literal['91k'],  # Only 91k supported ATM
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
    from nipype.interfaces.workbench import MetricResample
    from niworkflows.engine.workflows import LiterateWorkflow as Workflow
    from smriprep.interfaces.cifti import GenerateDScalar

    workflow = Workflow(name=name)
    workflow.__desc__ = f"""\
*Grayordinate* "dscalar" files [@hcppipelines] containing {grayord_density} samples were
also generated using `M-CRIB-S` as an intermediate standardized
surface space.
"""

    inputnode = pe.Node(
        niu.IdentityInterface(
            fields=[
                "subject_id",
                "subjects_dir",
                "surfaces",
                "morphometrics",
                "fsLR_midthickness",
            ]
        ),
        name="inputnode",
    )

    outputnode = pe.Node(
        niu.IdentityInterface(fields=["cifti_morph", "cifti_metadata"]),
        name="outputnode",
    )

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
    resample.inputs.new_sphere = [
        str(atlases / 'dHCP' / 'dHCP.week42.L.sphere.surf.gii'),
        str(atlases / 'dHCP' / 'dHCP.week42.R.sphere.surf.gii'),
    ] * 3
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
        (inputnode, resample, [
            ("fsLR_midthickness", "new_area"),
            (('surfaces', _get_surf, "midthickness", 3), "current_area")]),
        (inputnode, surfmorph_list, [
            (('morphometrics', _get_surf, "curv"), "in1"),
            (('morphometrics', _get_surf, "sulc"), "in2"),
            (('morphometrics', _get_surf, "thickness"), "in3"),
        ]),
        (resample, gen_cifti, [
            (("out_file", _collate), "scalar_surfs")]),
        (gen_cifti, outputnode, [("out_file", "cifti_morph"),
                                 ("out_metadata", "cifti_metadata")]),
    ])
    # fmt: on
    return workflow


def _get_surf(surfaces, name, mult=1):
    "Select a specific surface by name, and optionally multiple it."
    return [surf for surf in _sorted_by_basename(surfaces) if name in surf] * mult
