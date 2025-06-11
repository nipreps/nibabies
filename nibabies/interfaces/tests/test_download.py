from pathlib import Path

from nibabies.interfaces.download import RetrievePoochFiles


def test_RetrievePoochFiles(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    getter = RetrievePoochFiles(intermediate='MNIInfant:cohort-1', target='MNI152NLin6Asym')
    outputs = getter.run().outputs
    assert Path(outputs.int2tgt_xfm).exists()
    assert Path(outputs.tgt2int_xfm).exists()

    cache = tmp_path / 'mycache'
    monkeypatch.setenv('NIBABIES_POOCH_DIR', cache)
    getter = RetrievePoochFiles(intermediate='MNIInfant:cohort-1', target='MNI152NLin6Asym')
    outputs = getter.run().outputs

    assert Path(outputs.int2tgt_xfm) == cache / 'from-MNIInfant+1_to-MNI152NLin6Asym_xfm.h5'
    assert Path(outputs.tgt2int_xfm) == cache / 'from-MNI152NLin6Asym_to-MNIInfant+1_xfm.h5'
