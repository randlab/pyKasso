# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'pyKasso'
copyright = '2023, F. Miville, C. Fandel, P. Renard'
author = 'F. Miville, C. Fandel, P. Renard'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    # "sphinx.ext.autodoc",
    # "sphinx.ext.autosummary",
    # "sphinx.ext.napoleon",
    # "sphinx.ext.doctest",
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
html_static_path = ['_static']

html_theme_options = {
   "logo": {
      "image_light" : "pykasso_logo.png",
      "image_dark"  : "pykasso_logo.png",
   },
   "external_links": [
      {"name": "link-one-name", "url": "https://<link-one>"},
      {"name": "link-two-name", "url": "https://<link-two>"}
  ]
}