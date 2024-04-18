import json
from pathlib import Path

from niworkflows.data import load as nwf_load

from nibabies.data import load


def collect_anatomical_derivatives(
    derivatives_dir: Path | str,
    subject_id: str,
    std_spaces: list,
    spec: dict | None = None,
    patterns: list | None = None,
):
    """
    Collect outputs from across processing stages.

    Potential files:
    - T1w preproc
    - T2w preproc
    - T1w mask
    - T2w mask


    """
    from bids.layout import BIDSLayout

    if spec is None or patterns is None:
        _spec, _patterns = tuple(json.loads(load('io_spec_anat.json').read_text()).values())

        if spec is None:
            spec = _spec
        if patterns is None:
            patterns = _patterns

    deriv_config = nwf_load('nipreps.json')
    layout = BIDSLayout(derivatives_dir, config=deriv_config, validate=False)
    derivs_cache = {}

    for key, qry in spec['baseline'].items():
        qry['subject'] = subject_id
        item = layout.get(return_type='filename', **qry)
        if not item:
            continue

        derivs_cache[key] = item[0] if len(item) == 1 else item

    for key, qry in spec['coreg'].items():  # T1w->T2w, T2w->T1w
        qry['subject'] = subject_id
        item = layout.get(return_type='filename', **qry)
        if not item:
            continue
        derivs_cache[key] = item[0] if len(item) == 1 else item

    transforms = derivs_cache.setdefault('transforms', {})
    for _space in std_spaces:
        space = _space.replace(':cohort-', '+')
        for key, qry in spec['transforms'].items():
            qry = qry.copy()
            qry['subject'] = subject_id
            qry['from'] = qry['from'] or space
            qry['to'] = qry['to'] or space
            item = layout.get(return_type='filename', **qry)
            if not item:
                continue
            transforms.setdefault(_space, {})[key] = item[0] if len(item) == 1 else item

    for key, qry in spec['surfaces'].items():
        qry['subject'] = subject_id
        item = layout.get(return_type='filename', **qry)
        if not item or len(item) != 2:
            continue

        derivs_cache[key] = sorted(item)

    return derivs_cache
