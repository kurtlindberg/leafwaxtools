# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html


import leafwaxtools as leafwax
import sys
import os

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
package_path = os.path.abspath('../..')
os.environ['PYTHONPATH']=':'.join(((package_path), os.environ.get('PYTHONPATH','')))
sys.path.insert(0,os.path.abspath('../leafwaxtools'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information


project = 'leafwaxtools'
copyright = '2025, Kurt Lindberg'
author = 'Kurt Lindberg'
release = cat.__version__
version = cat.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions  = [
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.viewcode',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'numpydoc',
    'nbsphinx',
    'IPython.sphinxext.ipython_console_highlighting',
    'IPython.sphinxext.ipython_directive',
    'sphinx_search.extension',
    'jupyter_sphinx',
    'sphinx_copybutton'
]

source_suffix = '.rst'
templates_path = ['_templates']
exclude_patterns = []

#autosummary_generate = True
numpydoc_show_class_members = True

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
