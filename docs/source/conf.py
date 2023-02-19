# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project   = 'pyKasso'
author    = 'F. Miville, C. Fandel, P. Renard'
copyright = '2023, University of Neuchâtel - CHYN'
version   = '0.1.0'
release   = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
   'sphinx.ext.duration',
   'sphinx.ext.doctest',
   'sphinx.ext.autodoc',
   'sphinx.ext.autosummary',
    # "sphinx.ext.napoleon",
    # "sphinx.ext.intersphinx",
    # "sphinx.ext.mathjax",
    # "sphinx.ext.ifconfig",
    # "sphinx.ext.viewcode",
    # "IPython.sphinxext.ipython_console_highlighting",
    # "nbsphinx",
    # "numpydoc",
    # "sphinxcontrib.bibtex",
    # "sphinx_design",
    # "sphinx.ext.autosectionlabel",
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'

# html_theme_options = {}

html_static_path = ['_static']

html_theme_options = {
   "logo": {
      "image_light" : "pykasso_logo.png",
      "image_dark"  : "pykasso_logo.png",
   },
#    "external_links": [
#       {"name": "link-one-name", "url": "https://<link-one>"},
#       {"name": "link-two-name", "url": "https://<link-two>"}
#   ]
}



# -- Autosummary settings -------------------------------------------------

autosummary_generate = True

# autodoc_typehints = "description"
# autodoc_typehints_format = "short"



# autoclass_content = "class"

# autosectionlabel_prefix_document = True