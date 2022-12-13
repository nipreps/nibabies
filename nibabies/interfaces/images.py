#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""
Image tools interfaces
~~~~~~~~~~~~~~~~~~~~~~


"""
import numpy as np
import nibabel as nb
from textwrap import indent

from nipype import logging
from nipype.utils.filemanip import fname_presuffix
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec, SimpleInterface,
    File)
LOGGER = logging.getLogger('nipype.interface')

class RemapLabelsInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="input volume")
    dict = traits.Dict(exists=True, mandatory=True, key_trait=traits.BaseInt, value_trait=traits.BaseInt, desc="key-value pairs for remapping")
class RemapLabelsOutputSpec(TraitedSpec):
    out_file = File(Exists=True, desc="output volume")

class RemapLabels(SimpleInterface):
    """
    Remap values of the input NIfTI volume, using the key-value integer pairs 
    from the input dictionary.  
    """
    
    input_spec = RemapLabelsInputSpec
    output_spec = RemapLabelsOutputSpec
    
    def _run_interface(self, runtime):
        in_file = self.inputs.in_file
        dict = self.inputs.dict   
        self._results['out_file'] = fname_presuffix(self.inputs.in_file,
                                                    suffix='_remapped',
                                                    newpath=runtime.cwd)

        in_img = nb.load(in_file)
        in_data = in_img.get_fdata()
        out_data = np.vectorize(dict.get)(in_data)
        out_data[out_data == None] = 0
        nb.Nifti1Image(out_data.astype(int), in_img.affine, in_img.header).to_filename(self._results['out_file'])
        return runtime
        