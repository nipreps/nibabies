from __future__ import annotations

from fmriprep import config
from fmriprep.extensions.descriptor import ExtensionDescriptor

from nibabies._version import __version__
from nibabies.workflows.anatomical.fit import (
    init_infant_anat_fit_wf,
    init_infant_single_anat_fit_wf,
)

_MCRIBS_RECOMMEND_AGE_CAP = 3


class NibabiesExtension(ExtensionDescriptor):
    name = 'nibabies'
    version = __version__
    fmriprep_compat = '>=26,<27'
    contracts = {'anat_fit'}

    def cli_extend(self, parser) -> None:
        from nibabies.cli import extend_parser

        extend_parser(parser)

    def cli_populate(self, opts) -> None:
        from nibabies.cli import populate_config

        populate_config(opts, self)

    def init_config(self) -> None:
        super().init_config()
        if self.get('age_months') is None:
            from fmriprep.extensions.exceptions import ExtensionConfigError

            raise ExtensionConfigError('nibabies', '--age-months is required')

    def init_anat_fit_wf(self, **kwargs):
        age_months: int = self.get('age_months')
        recon_method = self.get('recon_method', 'auto')
        segmentation_atlases = self.get('segmentation_atlases')
        reference_anat = self.get('reference_anat')
        t1w: list = kwargs.get('t1w', [])
        t2w: list = kwargs.get('t2w', [])
        precomputed: dict = kwargs.get('precomputed', {})

        if recon_method == 'auto':
            if age_months <= _MCRIBS_RECOMMEND_AGE_CAP and precomputed.get('t2w_aseg'):
                recon_method = 'mcribs'
            elif age_months <= 24:
                recon_method = 'infantfs'
            else:
                recon_method = 'freesurfer'

        if not (t1w and t2w):
            single_anat = True
            derived_ref = 'T1w' if t1w else 'T2w'
            if reference_anat is not None and reference_anat != derived_ref:
                raise ValueError(
                    f'Requested --reference-anat {reference_anat} '
                    f'but only {derived_ref} images are available.'
                )
            reference_anat = derived_ref
        else:
            single_anat = False
            if reference_anat is None:
                reference_anat = (
                    'T2w'
                    if any(
                        (
                            recon_method is None and age_months <= _MCRIBS_RECOMMEND_AGE_CAP,
                            recon_method == 'mcribs',
                        )
                    )
                    else 'T1w'
                )

        nibabies_kwargs = dict(
            **kwargs,
            age_months=age_months,
            recon_method=recon_method,
            segmentation_atlases=segmentation_atlases,
            reference_anat=reference_anat,
            cifti_output=config.workflow.cifti_output,
        )

        if single_anat:
            return init_infant_single_anat_fit_wf(**nibabies_kwargs)
        return init_infant_anat_fit_wf(**nibabies_kwargs)
