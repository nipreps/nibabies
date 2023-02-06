import configparser
import os

from ..mirtk import ReconNeonatalCortex


def test_recon_neonatal_cortex_config(tmp_path, data_dir):
    os.chdir(tmp_path)

    ctx = ReconNeonatalCortex()
    config_file = ctx._create_config()
    config = _load_config(config_file)
    for section in [
        'recon-neonatal-cortex',
        'recon-neonatal-cortex white_model',
        'recon-neonatal-cortex pial_model',
    ]:
        assert section in config.keys()
    assert not config['recon-neonatal-cortex'].get('input_t1w_image')

    ctx.inputs.t1w_file = data_dir / 'sub-01_run-01_echo-1_bold.nii.gz'
    ctx.inputs.t2w_file = data_dir / 'sub-01_run-01_echo-2_bold.nii.gz'
    config_file = ctx._create_config()
    config = _load_config(config_file)
    assert config['recon-neonatal-cortex']['input_t1w_image']
    assert config['recon-neonatal-cortex']['input_t2w_image']


def _load_config(config_file: str) -> dict:
    cp = configparser.ConfigParser()
    with open(config_file) as fp:
        cp.read_file(fp)
    config = {s: dict(cp.items(s)) for s in cp.sections()}
    return config
