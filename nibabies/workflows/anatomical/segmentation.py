import sys
import typing as ty
from pathlib import Path

from nipype.interfaces import ants, fsl
from nipype.interfaces import utility as niu
from nipype.interfaces.ants.segmentation import JointFusion
from nipype.pipeline import engine as pe
from niworkflows.data import load as load_nwf
from niworkflows.engine import Workflow, tag
from niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from niworkflows.interfaces.fixes import FixHeaderRegistration as Registration
from niworkflows.utils.connections import listify
from smriprep.utils.misc import apply_lut
from smriprep.workflows.anatomical import (
    _aseg_to_three,
    _probseg_fast2bids,
    _split_segments,
)

from nibabies import config

LOGGER = config.loggers.workflow


@tag('anat.segmentation')
def init_segmentation_wf(
    *,
    sloppy: bool,
    method: ty.Literal['fast', 'jlf'] = 'fast',
    image_type: ty.Literal['T1w', 'T2w'] = 'T2w',
    jlf_template_dir: Path | None = None,
    omp_nthreads: int = 1,
    has_aseg: bool = False,
    name: str = 'segmentation_wf',
):
    workflow = Workflow(name=name)
    inputnode = pe.Node(
        niu.IdentityInterface(fields=['anat_brain', 'anat_aseg']),
        name='inputnode',
    )
    outputnode = pe.Node(
        niu.IdentityInterface(fields=['anat_dseg', 'anat_tpms', 'anat_aseg']),
        name='outputnode',
    )
    aseg_buffer = pe.Node(niu.IdentityInterface(fields=['anat_aseg']), name='aseg_buffer')

    to_dseg = pe.Node(
        niu.Function(function=apply_lut, output_names=['out_dseg']),
        name='to_dseg',
    )

    if has_aseg:
        LOGGER.info('ANAT Segmentation: Using existing segmentation')
        workflow.connect(inputnode, 'anat_aseg', aseg_buffer, 'anat_aseg')

    elif method == 'fast':
        workflow.__desc__ = (
            'Brain tissue segmentation of cerebrospinal fluid (CSF), white-matter (WM), and '
            f'gray-matter (GM) was performed on the brain-extracted {image_type} using FSL '
            f'FAST, distributed with {fsl.Info.version() or "version unknown"}'
        )
        # From FAST docs:
        # Number of classes to be segmented.
        # Normally you will want 3 (Grey Matter, White Matter and CSF).
        # However, if there is very poor grey/white contrast you may want to reduce this to 2;
        # alternatively, if there are strong lesions showing up as a fourth class, you may want
        # to increase this. Also, if you are segmenting T2-weighted images, you may need to
        # select 4 classes so that dark non-brain matter is processed correctly
        # (this is not a problem with T1-weighted as CSF and dark non-brain matter look similar).

        # fast_img_type = 1 if image_type == 'T1w' else 2
        # fast_n_classes = 3
        fast = pe.Node(
            fsl.FAST(
                segments=True,
                no_bias=True,
                probability_maps=True,
                # img_type=1 if image_type == 'T1w' else 2,
                # number_classes=fast_n_classes,
            ),
            name='fast',
            mem_gb=3,
        )

        to_dseg.inputs.lut = (0, 3, 1, 2)  # Maps: 0 -> 0, 3 -> 1, 1 -> 2, 2 -> 3.
        fast2bids = pe.Node(
            niu.Function(function=_probseg_fast2bids),
            name='fast2bids',
            run_without_submitting=True,
        )
        workflow.connect([
            (inputnode, fast, [('anat_brain', 'in_files')]),
            (fast, to_dseg, [('partial_volume_map', 'in_dseg')]),
            (to_dseg, outputnode, [('out_dseg', 'anat_dseg')]),
            (fast, fast2bids, [('partial_volume_files', 'inlist')]),
            (fast2bids, outputnode, [('out', 'anat_tpms')]),
        ])  # fmt:skip
        return workflow  # NOTE: no aseg will be output

    elif method == 'jlf':
        if not jlf_template_dir or not Path(jlf_template_dir).exists():
            raise RuntimeError('JLF requires a template directory.')

        jlf_wf = init_jlf_wf(
            jlf_template_dir=jlf_template_dir,
            sloppy=sloppy,
            image_type=image_type,
            omp_nthreads=omp_nthreads,
        )

        workflow.connect([
            (inputnode, jlf_wf, [('inputnode.anat_brain', 'inputnode.anat_brain')]),
            (jlf_wf, aseg_buffer, [('outputnode.anat_aseg', 'anat_aseg')]),
        ])  # fmt:skip

    to_dseg.inputs.lut = _aseg_to_three()
    split_seg = pe.Node(
        niu.Function(function=_split_segments, output_names=['out_tpms']),
        name='split_seg',
    )

    workflow.connect([
        (aseg_buffer, outputnode, [('anat_aseg', 'anat_aseg')]),
        (aseg_buffer, to_dseg, [('anat_aseg', 'in_dseg')]),
        (to_dseg, outputnode, [('out_dseg', 'anat_dseg')]),
        (to_dseg, split_seg, [('out_dseg', 'in_file')]),
        (split_seg, outputnode, [('out_tpms', 'anat_tpms')]),
    ])  # fmt:skip
    return workflow


@tag('anat.segmentation.jlf')
def init_jlf_wf(
    jlf_template_dir: Path,
    sloppy: bool,
    image_type: ty.Literal['T1w', 'T2w'] = 'T2w',
    omp_nthreads: int = 1,
    max_templates: int | None = None,
    name: str = 'jlf_wf',
):
    workflow = Workflow(name=name)
    inputnode = pe.Node(niu.IdentityInterface(fields=['anat_brain']), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(fields=['anat_aseg']))

    jlf_templates = _parse_jlf_templates(
        jlf_template_dir,
        image_type=image_type,
        max_templates=max_templates,
    )
    segmentations = jlf_templates.keys()
    references = jlf_templates.values()

    workflow.__desc__ = (
        f'The {image_type} image was registered to {len(segmentations)} templates for '
        f'JointFusion, distributed with ANTs {ants.base.Info.version() or "version unknown"}, '
        'for image segmentation. Brain tissue segmentation of cerebrospinal fluid (CSF), '
        'white-matter (WM), and gray-matter (GM) were derived from the label fused image.'
    )

    precision = 'testing' if sloppy else 'precise'
    norm_templates = pe.MapNode(
        Registration(from_file=load_nwf(f'antsBrainExtraction_{precision}.json')),
        name='norm_templates',
        iterfield=['moving_image'],
        n_procs=omp_nthreads,
        mem_gb=config.DEFAULT_MEMORY_MIN_GB,
    )
    norm_templates.inputs.moving_image = references
    norm_templates.inputs.float = True

    apply_template = pe.MapNode(
        ApplyTransforms(
            dimension=3,
            interpolation='NearestNeighbor',
            float=True,
        ),
        iterfield=['transform', 'input_image'],
        name='apply_template',
    )
    apply_template.inputs.input_image = norm_templates

    apply_seg = pe.MapNode(
        ApplyTransforms(dimension=3, interpolation='MultiLabel'),
        name='apply_seg',
        iterfield=['transforms', 'input_image'],
    )
    apply_seg.inputs.input_image = segmentations

    jointfusion = pe.Node(
        JointFusion(
            dimension=3,
            out_label_fusion='fusion_labels.nii.gz',
            num_threads=omp_nthreads,
        ),
        name='jointfusion',
    )
    clean_label_file = pe.Node(
        niu.Function(function=_to_dtype, output_names=['out_file']), name='clean_label_file'
    )
    workflow.connect([
        (inputnode, norm_templates, [('anat_brain', 'fixed_image')]),
        (norm_templates, apply_template, [('forward_transforms', 'transforms')]),
        (inputnode, apply_template, [('anat_brain', 'reference_image')]),
        (norm_templates, apply_seg, [('forward_transforms', 'transforms')]),
        (inputnode, apply_seg, [('anat_brain', 'reference_image')]),
        (inputnode, jointfusion, [(('anat_brain', listify), 'target_image')]),
        (apply_template, jointfusion, [('output_image', 'atlas_image')]),
        (apply_seg, jointfusion, [('output_image', 'atlas_segmentation_image')]),
        (jointfusion, clean_label_file, [('out_label_fusion', 'in_file')]),
        (clean_label_file, outputnode, [('out_file', 'anat_aseg')]),
    ])  # fmt:skip
    return workflow


def _parse_jlf_templates(
    templates_dir: Path | str,
    image_type: ty.Literal['T1w', 'T2w'] = 'T2w',
    max_templates: int | None = None,
):
    """
    Parse segmentation templates directory for anatomical and segmentation files.
    The segmentations are expected to follow the FreeSurfer LUT, and the anatomicals
    should be masked.

    This is compatible with the DCAN layout::

    jlf-templates/
    ├── Template01
    │   ├── Segmentation.nii.gz
    │   ├── T1w_brain.nii.gz
    │   └── T2w_brain.nii.gz
    ├── Template02
    ...

    And the BIDS layout::

    Templates/
    ├── dataset_description.json
    ├── sub-01
    │   ├── sub-01_desc-aseg_dseg.nii.gz
    │   ├── sub-01_T1w.json
    │   ├── sub-01_T1w.nii.gz
    │   ├── sub-01_T2w.json
    │   └── sub-01_T2w.nii.gz
    ├── sub-02
    ...

    """
    segmentations = {}
    templates_dir = Path(templates_dir)
    templates = [template.name for template in templates_dir.iterdir() if template.is_dir()]
    if not max_templates:
        max_templates = len(templates)

    if not templates:
        raise FileNotFoundError('JLF requested but no templates found.')

    for template in templates[:max_templates]:
        files = sorted((templates_dir / template).iterdir())
        seg = None
        anat = None
        for fl in files:
            if 'Segmentation' in fl.name or '_dseg' in fl.name:
                seg = str(fl)
            elif image_type in fl.name:
                anat = str(fl)
        if seg is None or anat is None:
            print(
                f'No anatomical or segmentation found for JLF template: {template}',
                file=sys.stderr,
            )
            continue
        segmentations[seg] = anat

    if len(segmentations) == 0:
        raise FileNotFoundError('JLF requested but anatomicals / segmentations were not found.')
    return segmentations


def _to_dtype(in_file, dtype='uint8'):
    """
    Freesurfer's ``mri_convert`` complains about unsigned 32-bit integers.
    Since we may use the JLF segmentation with FreeSurfer tools
    better to make this change now.
    """
    from pathlib import Path

    import nibabel as nb
    import numpy as np

    img = nb.load(in_file)
    out_file = Path(f'labels{"".join(Path(in_file).suffixes)}').absolute()

    new_data = np.asanyarray(img.get_fdata(), dtype=dtype)
    img.set_data_dtype(dtype)
    img.__class__(new_data, img.affine, img.header).to_filename(out_file)
    return str(out_file)
