from nipype.interfaces.base import (
    TraitedSpec,
    BaseInterfaceInputSpec,
    File,
    SimpleInterface,
)

from sdcflows.utils.epimanip import epi_mask


class _EPIMaskInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc="File to be masked")


class _EPIMaskOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="Masked EPI")


class EPIMask(SimpleInterface):
    """Mask an EPI image"""

    input_spec = _EPIMaskInputSpec
    output_spec = _EPIMaskOutputSpec

    def _run_interface(self, runtime):
        self._results["out_file"] = epi_mask(self.inputs.in_file)
