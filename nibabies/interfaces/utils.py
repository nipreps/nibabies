import os
import re

from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    File,
    InputMultiObject,
    OutputMultiObject,
    SimpleInterface,
    TraitedSpec,
    traits,
)


class CiftiSelectInputSpec(BaseInterfaceInputSpec):
    hemi = traits.Enum('L', 'R', desc='Hemisphere')
    surfaces = InputMultiObject(File(exists=True), desc='Surfaces')
    morphometrics = InputMultiObject(File(exists=True), desc='Surface morphometrics')
    spherical_registrations = InputMultiObject(
        File(exists=True), desc='Spherical registration to fsLR'
    )
    template_spheres = InputMultiObject(File(exists=True), desc='fsLR sphere')
    template_surfaces = InputMultiObject(File(exists=True), desc='fsLR midthickness')
    template_rois = InputMultiObject(File(exists=True), desc='fsLR ROIs')


class CiftiSelectOutputSpec(TraitedSpec):
    white = OutputMultiObject(File, desc='white surface')
    pial = OutputMultiObject(File, desc='pial surface')
    midthickness = OutputMultiObject(File, desc='midthickness surface')
    thickness = OutputMultiObject(File, desc='thickness surface')
    sphere_reg = OutputMultiObject(File, desc='fsLR spherical registration')
    template_sphere = OutputMultiObject(File, desc='fsLR sphere')
    template_surface = OutputMultiObject(File, desc='fsLR surface (midthickness)')
    template_roi = OutputMultiObject(File, desc='fsLR ROIs')


class CiftiSelect(SimpleInterface):
    input_spec = CiftiSelectInputSpec
    output_spec = CiftiSelectOutputSpec

    def _run_interface(self, runtime):
        idx = 0 if self.inputs.hemi == 'L' else 1
        all_surfaces = (self.inputs.surfaces or []) + (self.inputs.morphometrics or [])
        container = {
            'white': [],
            'pial': [],
            'midthickness': [],
            'thickness': [],
            'sphere_reg': self.inputs.spherical_registrations or [],
            'template_sphere': self.inputs.template_spheres or [],
            'template_surface': self.inputs.template_surfaces or [],
            'template_roi': self.inputs.template_rois or [],
        }
        find_name = re.compile(r'(?:^|[^d])(?P<name>white|pial|midthickness|thickness)')
        for surface in all_surfaces:
            match = find_name.search(os.path.basename(surface))
            if match:
                container[match.group('name')].append(surface)

        for name, vals in container.items():
            if vals:
                self._results[name] = sorted(vals, key=os.path.basename)[idx]
        return runtime
