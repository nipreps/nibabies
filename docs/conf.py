# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
from datetime import datetime


# -- Project information -----------------------------------------------------

project = 'NiBabies'
start_year = 2020

author = f'The {project} Developers'
copyright = f'{start_year}-{datetime.now().year}, {author}'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    # "sphinx.ext.autodoc",
    # "nipype.sphinxext.plot_workflow",
    # "sphinxcontrib.napoleon",
    "sphinxarg.ext",  # argparse extension
    "myst_nb",
]

autodoc_mock_imports = [
    "numpy",
    "nibabel",
    "nilearn"
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

myst_heading_anchors = 2
