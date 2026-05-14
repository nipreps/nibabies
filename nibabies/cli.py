"""CLI argument extension for the NiBabies fmriprep extension."""

from __future__ import annotations

from argparse import ArgumentTypeError, BooleanOptionalAction
from pathlib import Path


def _str_none(val: str) -> str | None:
    """Convert the literal strings 'none' / 'None' to Python ``None``."""
    if not isinstance(val, str):
        return val
    return None if val.lower() == 'none' else val


def _dir_not_empty(value: str) -> Path:
    path = Path(value)
    if not path.is_dir():
        raise ArgumentTypeError(f'Path is not a directory: <{value}>')
    if not any(path.iterdir()):
        raise ArgumentTypeError(f'Directory is empty: <{value}>')
    return path


def _pos_int_or_auto(value: str) -> str | int:
    if value.lower() == 'auto':
        return 'auto'
    try:
        v = int(value)
        if v < 0:
            raise ValueError
        return v
    except ValueError as err:
        raise ArgumentTypeError(
            f"--hmc-bold-frame must be either 'auto' or a positive integer. Received {value}."
        ) from err


def extend_parser(parser) -> None:
    """Add a *NiBabies options* argument group to *parser* in-place.

    Called by :meth:`~nibabies.descriptor.NibabiesExtension.cli_extend`
    with fmriprep's fully-built :class:`~argparse.ArgumentParser`.
    """
    g = parser.add_argument_group('NiBabies options')
    g.add_argument(
        '--age-months',
        dest='age_months',
        type=int,
        help='Age of the participant in months.',
    )
    g.add_argument(
        '--segmentation-atlases-dir',
        dest='segmentation_atlases_dir',
        type=_dir_not_empty,
        default=None,
        help='Directory containing precalculated segmentations for JointLabelFusion.',
    )
    g.add_argument(
        '--fd-radius',
        dest='fd_radius',
        type=float,
        default=45,
        help='Head radius in mm for framewise displacement calculation.',
    )
    g.add_argument(
        '--surface-recon-method',
        dest='surface_recon_method',
        choices=('auto', 'infantfs', 'freesurfer', 'mcribs', None),
        type=_str_none,
        default='auto',
        help="Method to use for surface reconstruction. Use 'auto' to select based on age (default).",
    )
    g.add_argument(
        '--reference-anat',
        '--reference-anatomical',
        dest='reference_anat',
        choices=('T1w', 'T2w'),
        default=None,
        help='Override which anatomical to use as the structural reference. '
        'Determined automatically from age and recon-method when not set.',
    )
    g.add_argument(
        '--hmc-bold-frame',
        dest='hmc_bold_frame',
        type=_pos_int_or_auto,
        default=16,
        metavar='{auto,FRAME_NUMBER}',
        help='Frame to start head motion estimation on BOLD. '
        '``auto`` chooses this frame based on a sum-of-least-squares heuristic.',
    )
    g.add_argument(
        '--norm-csf',
        dest='norm_csf',
        action='store_true',
        help='Replace low intensity voxels in CSF mask with average.',
    )
    g.add_argument(
        '--multi-step-reg',
        dest='multi_step_reg',
        action=BooleanOptionalAction,
        default=True,
        help='For certain adult templates (MNI152NLin6Asym), perform two-step '
        'registrations (native -> MNIInfant -> template) and concatenate into a single xfm.',
    )


def populate_config(opts, ext) -> None:
    """Write parsed CLI values into the extension's config namespace.

    Called by :meth:`~nibabies.descriptor.NibabiesExtension.cli_populate`
    after argparse runs. ``opts`` is the parsed
    :class:`~argparse.Namespace`; ``ext`` is the
    :class:`~nibabies.descriptor.NibabiesExtension` instance.

    The ``--surface-recon-method`` flag is stored under the key
    ``recon_method`` to match the parameter name used by the nibabies
    workflow builders.
    """
    ext.set('age_months', getattr(opts, 'age_months', None))
    ext.set('segmentation_atlases', getattr(opts, 'segmentation_atlases_dir', None))
    ext.set('fd_radius', getattr(opts, 'fd_radius', 45))
    ext.set('recon_method', getattr(opts, 'surface_recon_method', 'auto'))
    ext.set('reference_anat', getattr(opts, 'reference_anat', None))
    ext.set('hmc_bold_frame', getattr(opts, 'hmc_bold_frame', 16))
    ext.set('norm_csf', getattr(opts, 'norm_csf', False))
    ext.set('multi_step_reg', getattr(opts, 'multi_step_reg', True))
