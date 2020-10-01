from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from .brain_extraction import init_infant_brain_extraction_wf
from .surfaces import init_infant_surface_recon_wf


def init_infant_anat_wf(
    *,
    template_name,
    template_specs,
    age_months,
    mri_scheme,
    omp_nthreads,
    output_dir,
    subject_id,
    name='infant_anat_wf',
):

    wf = pe.Workflow(name=name)
    # inputnode = pe.Node(niu.IdentityInterface(fields=["in_files", "in_seg"]), name='inputnode')
    inputnode = pe.Node(niu.IdentityInterface(fields=["in_files"]), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=["anat_corrected", "anat_brain", "anat_mask", "subjects_dir"]
    ), name='outputnode')

    brain_extraction_wf = init_infant_brain_extraction_wf(
        ants_affine_init=True,
        in_template=template_name,
        template_specs=template_specs,
        mri_scheme=mri_scheme,
        omp_nthreads=omp_nthreads,
        output_dir=output_dir,
    )

    wf.connect([
        (inputnode, brain_extraction_wf, [
            ('in_files', 'inputnode.in_files'),
            # ('in_seg', 'inputnode.in_seg')
        ]),
        (brain_extraction_wf, outputnode, [
            ('outputnode.out_corrected', 'anat_corrected'),
            ('outputnode.out_brain', 'anat_brain'),
            ('outputnode.out_mask', 'anat_mask')
        ]),
    ])

    surface_recon_wf = init_infant_surface_recon_wf(
        age_months=age_months,
        output_dir=output_dir,
        subject_id=subject_id,
    )

    wf.connect([
        # (inputnode, surface_recon_wf, [('in_seg', 'inputnode.in_seg')]),
        (brain_extraction_wf, surface_recon_wf, [('outputnode.out_brain', 'inputnode.masked_file')]),
        (surface_recon_wf, outputnode, [('outputnode.subjects_dir', 'subjects_dir')])
    ])

    return wf