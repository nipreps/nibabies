# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime
from sphinx import __version__ as sphinxversion
from packaging.version import Version

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.append(os.path.abspath("sphinxext"))
sys.path.insert(0, os.path.abspath("../wrapper"))

from github_link import make_linkcode_resolve  # this is only available after sphinxext to PATH

# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = "1.5.3"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "sphinx.ext.intersphinx",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.linkcode",
    "sphinx.ext.napoleon",
    "sphinxarg.ext",  # argparse extension
    "nipype.sphinxext.plot_workflow",
    "myst_nb",  # stop segregating rst/md
]

autodoc_mock_imports = [
    "numpy",
    "nibabel",
    "nilearn"
]
if Version.parse(sphinxversion) >= Version.parse("1.7.0"):
    autodoc_mock_imports += [
        "pandas",
        "nilearn",
        "seaborn",
    ]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


source_suffix = [".rst", ".md"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

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
    ("Inputs", "Parameters"),
    ("Outputs", "Parameters"),
]

# -- MyST parameters ---------------------------------------------------------

myst_heading_anchors = 3
myst_enable_extensions = [
    "colon_fence",
    "substitution",
]

project = "NiBabies"
author = "The NiPreps developers"
copyright = "2021-%s, %s" % (datetime.now().year, author)

import nibabies

nibabies_ver = Version(nibabies.__version__)
release = "version" if nibabies_ver.is_prerelease else nibabies_ver.public

myst_substitutions = {
    "release": release,
    "version": str(nibabies_ver),
    "dockerbuild": "docker pull nipreps/nibabies:{{ release }}",
    "singbuild": (
        "singularity build nibabies-{{ release }}.sif docker://nipreps/nibabies:{{ release }}"
    ),
}

# to avoid Python highlighting in literal text
highlight_language = "none"
