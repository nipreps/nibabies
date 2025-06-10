import json
import os
from pathlib import Path

import pooch
from nipype.interfaces.base import (
    DynamicTraitedSpec,
    File,
    SimpleInterface,
    TraitedSpec,
    traits,
)

import nibabies


class _RetrievePoochFilesInputSpec(DynamicTraitedSpec):
    intermediate = traits.Str(required=True, desc='the intermediate space')
    target = traits.Str(required=True, desc='the target space')


class _RetrievePoochFilesOutputSpec(TraitedSpec):
    int2tgt_xfm = File(desc='Intermediate to target transform')
    tgt2int_xfm = File(desc='Target to intermediate transform')


class RetrievePoochFiles(SimpleInterface):
    input_spec = _RetrievePoochFilesInputSpec
    output_spec = _RetrievePoochFilesOutputSpec

    def _run_interface(self, runtime):
        int2tgt, tgt2int = _retrieve_xfms(self.inputs.intermediate, self.inputs.target)
        self._results['int2tgt_xfm'] = int2tgt
        self._results['tgt2int_xfm'] = tgt2int
        return runtime


def _retrieve_xfms(
    intermediate: str,
    target: str,
):
    """Fetch transforms from the OSF repository (https://osf.io/y763j/)."""

    manifest = json.loads(nibabies.data.load('xfm_manifest.json').read_text())

    def sanitize(space):
        # MNIInfant:cohort-1 -> MNIInfant+1
        return space.replace(':cohort-', '+')

    intmd = sanitize(intermediate)
    tgt = sanitize(target)

    cache_dir = Path(os.getenv('NIBABIES_POOCH_DIR', Path.cwd()))

    int2std_name = f'from-{intmd}_to-{tgt}_xfm.h5'
    int2std_meta = manifest[int2std_name]
    int2std = pooch.retrieve(
        url=int2std_meta['url'],
        path=cache_dir,
        known_hash=int2std_meta['hash'],
        fname=int2std_name,
    )

    std2int_name = f'from-{tgt}_to-{intmd}_xfm.h5'
    std2int_meta = manifest[std2int_name]
    std2int = pooch.retrieve(
        url=std2int_meta['url'],
        path=cache_dir,
        known_hash=std2int_meta['hash'],
        fname=std2int_name,
    )

    return int2std, std2int
