# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Base module variables."""

try:
    from ._version import __version__
except ImportError:
    __version__ = '0+unknown'

__org__ = 'nipreps'
__packagename__ = 'nibabies'
__copyright__ = 'Copyright 2023, Center for Reproducible Neuroscience, Stanford University'
__credits__ = (
    'Contributors: please check the ``.zenodo.json`` file at the top-level folder'
    'of the repository'
)
__url__ = f'https://github.com/{__org__}/{__packagename__}'

DOWNLOAD_URL = f'https://github.com/{__org__}/{__packagename__}/archive/{__version__}.tar.gz'
