"""Sphinx configuration."""
project = "frogr"
author = "M. F. Bieger"
copyright = "2022, M. F. Bieger"
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx_click",
    "myst_parser",
]
autodoc_typehints = "description"
html_theme = "furo"
