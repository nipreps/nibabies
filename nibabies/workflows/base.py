from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu

from .brain_extraction import init_infant_brain_extraction_wf
from .surfaces import init_infant_surface_recon_wf


def init_infant_anat_wf(
    *,
    in_template,
    template_specs,
    age_months,
    mri_scheme,
    omp_nthreads,
    output_dir
):
    inputnode = pe.Node(niu.IdentityInterface(fields=["in_files"], ["in_seg"]), name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=["anat_corrected", "anat_brain", "anat_mask", "surfaces", "anat_aseg", "anat_aparc"]
    ), name='outputnode')

    brain_extraction_wf = init_infant_brain_extraction_wf(
        ants_affine_init=True,
        in_template=opts.template,
        template_specs=template_specs,
        mri_scheme=opts.mri_scheme,
        omp_nthreads=opts.omp_nthreads,
        output_dir=opts.output_dir,
    )

    wf.connect([
        (inputnode, brain_extraction_wf, [('in_files', 'inputnode.in_files'),
                                          ('in_seg', 'inputnode.in_seg')]),
    ])

    surface_recon_wf = init_infant_surface_recon_wf(
        age_months=age_months,
    )

    wf.connect([
        (inputnode, surface_recon_wf, [('in_files', 'inputnode.t1_file'),
                                       ('in_seg', 'inputnode.in_seg')]),
        (brain_extraction_wf, surface_recon_wf, [('out_brain', 'inputnode.masked_file')]),
    ])
