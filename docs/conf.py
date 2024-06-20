# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime, timezone

from packaging.version import Version, parse
from sphinx import __version__ as sphinxversion

import nibabies

# -- Path setup --------------------------------------------------------------
here = os.path.dirname(__file__)
# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.append(os.path.join(here, 'sphinxext'))
sys.path.insert(0, os.path.join(here, '..', 'wrapper'))

# this is only available after sphinxext to PATH
from github_link import make_linkcode_resolve  # noqa: E402

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '1.5.3'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.intersphinx',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.linkcode',
    'sphinx.ext.napoleon',
    'sphinxcontrib.bibtex',
    'sphinxarg.ext',  # argparse extension
    'nipype.sphinxext.plot_workflow',
    'myst_parser',  # allow markdown
    # 'sphinx-togglebutton',  # collapse admonitions
]

bibtex_bibfiles = ['../nibabies/data/boilerplate.bib']

autodoc_mock_imports = ['numpy', 'nibabel', 'nilearn']
if parse(sphinxversion) >= parse('1.7.0'):
    autodoc_mock_imports = [
        'pandas',
        'nilearn',
        'seaborn',
    ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


source_suffix = ['.rst', '.md']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'shibuya'

# Options specific to theme
html_theme_options = {
    'color_mode': 'light',
    'dark_code': True,
    'github_url': 'https://github.com/nipreps/nibabies',
    'nav_links': [
        {
            'title': 'NiPreps Homepage',
            'url': 'https://nipreps.org',
            'external': True,
        },
        {
            'title': 'Docker Hub',
            'url': 'https://hub.docker.com/r/nipreps/nibabies',
            'external': True,
        },
    ],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# -- Napoleon parameters -----------------------------------------------------
# Accept custom section names to be parsed for numpy-style docstrings
# of parameters.
# Requires pinning sphinxcontrib-napoleon to a specific commit while
# https://github.com/sphinx-contrib/napoleon/pull/10 is merged.
napoleon_use_param = False
napoleon_custom_sections = [
    ('Inputs', 'Parameters'),
    ('Outputs', 'Parameters'),
]

# -- MyST parameters ---------------------------------------------------------

myst_heading_anchors = 3
myst_enable_extensions = [
    'colon_fence',
    'substitution',
]

linkcode_resolve = make_linkcode_resolve(
    'nibabies',
    'https://github.com/nipreps/' 'nibabies/blob/{revision}/' '{package}/{path}#L{lineno}',
)

project = 'NiBabies'
author = 'The NiPreps developers'

copyright = f'2021-{datetime.now(tz=timezone.utc)}, {author}'

nibabies_ver = Version(nibabies.__version__)
release = 'version' if nibabies_ver.is_prerelease else nibabies_ver.public

# to avoid Python highlighting in literal text
highlight_language = 'none'
