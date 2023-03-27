# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import date
from pykasso import __version__

year = date.today().strftime("%Y")
print(os.path.abspath("."))
sys.path.insert(0, os.path.abspath("."))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project   = 'pyKasso'
author    = 'F. Miville, C. Fandel, P. Renard'
copyright = '{}, University of Neuchâtel - CHYN'.format(year)
version   = __version__
release   = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
   'sphinx.ext.duration',
   'sphinx.ext.doctest',
   'sphinx.ext.autodoc',
   'sphinx.ext.autosummary',
   'sphinx.ext.napoleon',
   'sphinx.ext.coverage',
   'sphinx.ext.intersphinx',
   'numpydoc',
    # 'sphinx_design',
    # "sphinx.ext.mathjax",
    # "sphinx.ext.ifconfig",
    # "sphinx.ext.viewcode",
    # "IPython.sphinxext.ipython_console_highlighting",
    # "nbsphinx",
    # "sphinxcontrib.bibtex",
    # "sphinx.ext.autosectionlabel",
]

exclude_patterns = []
templates_path = ['_templates']

# add_function_parentheses = False
# add_module_names = False
# toc_object_entries_show_parents = 'all'





# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'

html_static_path = ['_static']

html_theme_options = {
   "logo": {
      "image_light" : "pykasso_logo.png",
      "image_dark"  : "pykasso_logo.png",
   },
   "external_links": [
      {"name": "github-repo", "url": "https://github.com/randlab/pyKasso"},
  ]
}

# -- Autosummary settings -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/autodoc.html#confval-autoclass_content

autoclass_content = "class"

# autodoc_member_order = "bysource"

autodoc_default_options = {
   #  'members': 'var1, var2',
    'member-order': 'bysource',
    'special-members': '__init__',
   #  'undoc-members': True,
   #  'exclude-members': '__weakref__'
}

# autosummary_generate = True

autodoc_typehints = "description"

autodoc_typehints_format = "fully-qualified" #"short"

# autosectionlabel_prefix_document = True

# -- Sets intersphinx Directories ------------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/extensions/intersphinx.html

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/devdocs", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
}

# -- Napoleon settings ----------------------------------------------------------------

# napoleon_include_init_with_doc = True
# napoleon_use_param = True
# napoleon_type_aliases = {
#     "pk": "pykasso",
# }

# -- Numpydoc settings ----------------------------------------------------------------
# https://numpydoc.readthedocs.io/en/v1.5.0/install.html

# numpydoc_xref_param_type = True

# numpydoc_xref_aliases = {
#    'pk' : 'pykasso',
#    'RandomNumberGenerator' : 'numpy.random._generator.Generator',
#    'DataFrame' : 'pandas.core.frame.DataFrame',
   # '':'',
# }