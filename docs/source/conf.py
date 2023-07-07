# Configuration file for the Sphinx documentation builder.
#
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('../../bin'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AmesCAP'
copyright = '2023, Alex Kling, Courtney Batterson, & Victoria Hartwick (Mars Climate Modeling Center | NASA Ames Research Center)'
author = 'Alex Kling, Courtney Batterson, & Victoria Hartwick (Mars Climate Modeling Center | NASA Ames Research Center)'
release = '1.0'
master_doc = 'index'
root_doc = 'index'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'autoapi.extension',
]

autoapi_type = 'python'
autoapi_dirs = ['../../bin']
autoapi_add_toctree_entry = False
autoapi_keep_files = True
# autoapi_options = [ 'members', 'inherited-members', 'show-inheritance', 'show-inheritance-diagram', 'show-module-summary', 'imported-members' ]
autoapi_options = [ 'members', 'inherited-members', 'show-inheritance', 'show-module-summary', 'imported-members' ]

pygments_style = 'sas'

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

# html_theme = 'furo'
html_theme = "sphinx_rtd_theme"

# Options: alabaster, classic, sphinxdoc, bizstyle, furo

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

