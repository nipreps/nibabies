import argparse
from collections import defaultdict
import json
from pathlib import Path
import warnings
warnings.simplefilter("ignore", UserWarning)

import nibabel as nb
import nitransforms as nt
import numpy as np


DEFAULT_OUTPUT_SPACES = {"MNIInfant",}

CHECKMARK = u' \u2713'


def get_parser():

    # class AddSpaces(argparse.Action):
    #     _spaces = {"MNIInfant", }

    #     def __call__(self, parser, args, values, option_string=None):
    #         for value in values:
    #             if "cohort" in value:  # skip cohorts for now
    #                 continue
    #             self._spaces.add(value)
    #         setattr(args, self.dest, self._spaces)

    parser = argparse.ArgumentParser(description="Script to verify NiBabies outputs")
    parser.add_argument("output_dir", type=Path, help="Directory containing nibabies outputs")
    return parser


def main(argv=None):
    parser = get_parser()
    pargs = parser.parse_args(argv)
    test_output_layout(pargs.output_dir)
    for anats in pargs.output_dir.glob("sub-*/**/anat"):
        if not anats.is_dir():
            continue
        test_anatomicals(anats)

    for funcs in pargs.output_dir.glob("sub-*/**/func"):
        if not funcs.is_dir():
            continue
        test_functionals(funcs)
    print(f"No errors found for output directory {str(pargs.output_dir)}")


def test_output_layout(output_dir):
    assert output_dir.exists()
    dd = output_dir / 'dataset_description.json'
    assert dd.exists()
    assert 'NiBabies' in json.loads(dd.read_text())['Name']
    for ext in ('bib', 'html', 'md', 'tex'):
        assert (output_dir / "logs" / f"CITATION.{ext}").exists()


def test_anatomicals(anat_dir):
    print(f"** Checking anatomical derivatives ({str(anat_dir)}):")
    outputs = defaultdict(list)
    for fl in anat_dir.glob("*"):
        outtype = None
        if 'desc-brain_mask' in fl.name:
            outtype = 'masks'
        elif 'desc-preproc_T1w' in fl.name:
            outtype = 'preprocs'
        elif '_dseg.' in fl.name or '_probseg.' in fl.name:
            outtype = 'segs'
        elif '_xfm.' in fl.name:
            outtype = 'xfms'
        elif 'surf.gii' in fl.name:
            outtype = 'surfs'

        if outtype:
            outputs[outtype].append(fl)

    _check_masks(outputs['masks'])
    _check_segs(outputs['segs'])
    _check_xfms(outputs['xfms'])
    _check_surfs(outputs['surfs'])
    _check_t1w_preprocs(outputs['preprocs'])


def test_functionals(func_dir):
    print(f"** Checking functional derivatives ({str(func_dir)}):")
    outputs = defaultdict(list)
    for fl in func_dir.glob("*"):
        outtype = None
        if 'desc-brain_mask' in fl.name:
            outtype = 'masks'
        elif 'desc-preproc_bold' in fl.name:
            outtype = 'preprocs'
        elif '_boldref.' in fl.name:
            outtype = 'boldrefs'
        elif 'desc-aseg_dseg' in fl.name or '_probseg.' in fl.name:
            outtype = 'segs'
        elif '_xfm.' in fl.name:
            outtype = 'xfms'
        elif '.surf.gii' in fl.name:
            outtype = 'surfs'
        elif 'bold.dtseries' in fl.name:
            outtype = 'ciftis'

        if outtype:
            outputs[outtype].append(fl)

    _check_bold_preprocs(outputs['preprocs'])
    _check_boldrefs(outputs['boldrefs'])
    _check_masks(outputs['masks'])
    _check_segs(outputs['segs'])
    _check_xfms(outputs['xfms'])
    _check_surfs(outputs['surfs'])
    _check_ciftis(outputs['ciftis'])


def _check_masks(masks):
    for mask in masks:
        print(str(mask), end='')
        if mask.name.endswith('.json'):
            metadata = json.loads(mask.read_text())
            assert metadata
            # assert metadata["Type"] == "Brain"
        else:
            img = nb.load(mask)
            assert img.dataobj.dtype == np.uint8
            assert np.all(np.unique(img.dataobj) == [0, 1])
        print(CHECKMARK)


def _check_segs(segs):
    for seg in segs:
        print(str(seg), end='')
        img = nb.load(seg)
        if 'desc-aseg' in seg.name:
            assert img.dataobj.dtype == np.int16
            labels = set(np.unique(img.dataobj))
            # check for common subcortical labels
            subcor = {26, 58, 18, 54, 16, 11, 50, 8, 47, 28, 60, 17, 53, 13, 52, 12, 51, 10, 49}
            missing = subcor.difference(labels)
            if missing:
                print(f"Missing labels {missing}")
        elif 'desc-aparcaseg' in seg.name:
            assert img.dataobj.dtype == np.int16
            labels = set(np.unique(img.dataobj))
            # check for common cortical labels
            cort = {1000, 1007, 1022, 2000, 2007, 2022}
            missing = cort.difference(labels)
            if missing:
                print(f"Missing labels {missing}")
        elif '_probseg.' in seg.name:
            assert img.dataobj.dtype == np.float32
            assert np.max(img.dataobj) == 1
            assert np.min(img.dataobj) == 0
        print(CHECKMARK)


def _check_xfms(xfms):
    for xfm in xfms:
        print(str(xfm), end='')
        if xfm.name.endswith(".txt"):
            assert nt.linear.load(xfm, fmt='itk')
        elif xfm.name.endswith('.h5'):
            assert nt.manip.load(xfm, fmt='itk')
        elif xfm.name.endswith('.json'):
            meta = json.loads(xfm.read_text())
            assert 'Sources' in meta
        else:
            raise NotImplementedError
        print(CHECKMARK)


def _check_surfs(surfs):
    for surf in surfs:
        print(str(surf), end='')
        if surf.name.endswith('.surf.gii'):
            img = nb.load(surf)
            assert img.numDA == 2
            da0, da1 = img.darrays
            assert da0.intent == 1008  # NIFTI_INTENT_POINTSET
            assert da1.intent == 1009  # NIFTI_INTENT_TRIANGLE
        print(CHECKMARK)


def _check_t1w_preprocs(preprocs):
    for preproc in preprocs:
        print(str(preproc), end='')
        if '.json' in preproc.name:
            metadata = json.loads(preproc.read_text())
            assert 'SkullStripped' in metadata
        else:
            img = nb.load(preproc)
            assert len(img.shape) == 3
        print(CHECKMARK)


def _check_bold_preprocs(preprocs):
    for preproc in preprocs:
        print(str(preproc), end='')
        if '.json' in preproc.name:
            metadata = json.loads(preproc.read_text())
            assert metadata['SkullStripped'] is False
        else:
            img = nb.load(preproc)
            assert len(img.shape) == 4
        print(CHECKMARK)


def _check_boldrefs(boldrefs):
    for boldref in boldrefs:
        print(str(boldref), end='')
        if '.json' in boldref.name:
            metadata = json.loads(boldref.read_text())
            assert 'Sources' in metadata
        else:
            bimg = nb.load(boldref)
            if len(bimg.shape) > 3:
                assert len(bimg.shape) == 4 and bimg.shape[3] == 1
        print(CHECKMARK)


def _check_ciftis(ciftis):
    for cifti in ciftis:
        print(str(cifti), end='')
        if '.json' in cifti.name:
            metadata = json.loads(cifti.read_text())
            assert metadata
        else:
            img = nb.load(cifti)
            assert len(img.shape) == 2
            matrix = img.header.matrix
            assert matrix.mapped_indices == [0, 1]
            series_map = matrix.get_index_map(0)
            bm_map = matrix.get_index_map(1)
            assert series_map.indices_map_to_data_type == 'CIFTI_INDEX_TYPE_SERIES'
            assert bm_map.indices_map_to_data_type == 'CIFTI_INDEX_TYPE_BRAIN_MODELS'
        print(CHECKMARK)


if __name__ == "__main__":
    main()
