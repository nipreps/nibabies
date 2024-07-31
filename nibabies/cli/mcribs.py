import os
from argparse import ArgumentParser

import nipype.pipeline.engine as pe
from niworkflows.interfaces.nibabel import MapLabels

from nibabies.interfaces.mcribs import MCRIBReconAll


def _parser():
    parser = ArgumentParser(description='Test script for MCRIBS surfaces')
    parser.add_argument('subject', help='Subject ID')
    parser.add_argument('t2w', type=os.path.abspath, help='Input T2w (radioisotropic)')
    parser.add_argument(
        'segmentation', type=os.path.abspath, help='Input anatomical segmentation in T2w space'
    )
    parser.add_argument(
        '--outdir', type=os.path.abspath, help='Output directory to persist MCRIBS output'
    )
    parser.add_argument('--nthreads', type=int, help='Number of threads to parallelize tasks')
    return parser


def main(argv: list = None):
    pargs = _parser().parse_args(argv)

    t2w_file = _check_file(pargs.t2w)
    seg_file = _check_file(pargs.segmentation)

    aseg2mcrib = {
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
    map_labels = pe.Node(MapLabels(in_file=seg_file, mappings=aseg2mcrib), name='map_labels')

    recon = pe.Node(
        MCRIBReconAll(subject_id=pargs.subject, t2w_file=t2w_file), name='mcribs_recon'
    )
    if pargs.outdir:
        recon.inputs.outdir = pargs.outdir
    if pargs.nthreads:
        recon.inputs.nthreads = pargs.nthreads

    wf = pe.Workflow(f'MRA_{pargs.subject}')
    wf.connect(map_labels, 'out_file', recon, 'segmentation_file')
    wf.run()


def _check_file(fl: str) -> str:
    import nibabel as nb
    import numpy as np

    img = nb.load(fl)
    if len(img.shape) != 3:
        raise ValueError('Image {fl} is not 3 dimensional.')

    voxdims = img.header['pixdim'][1:4]
    if not np.allclose(voxdims, voxdims[1]):
        raise ValueError(f'Image {fl} is not isotropic: {voxdims}.')

    ornt = nb.io_orientation(img.affine)
    axcodes = nb.orientations.ornt2axcodes(ornt)
    if ''.join(axcodes) != 'LAS':
        las = nb.orientations.axcodes2ornt('LAS')
        transform = nb.orientations.ornt_transform(ornt, las)
        reornt = img.as_reoriented(transform)
        outfl = os.path.abspath(f'LASornt_{os.path.basename(fl)}')
        print(f'Creating reorientated image {outfl}')
        reornt.to_filename(outfl)
        return outfl
    return fl


if __name__ == '__main__':
    main()
