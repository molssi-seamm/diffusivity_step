# -*- coding: utf-8 -*-

"""
diffusivity_step
A SEAMM plug-in for Diffusivity
"""

# Bring up the classes so that they appear to be directly in
# the diffusivity_step package.

from .diffusivity import Diffusivity
from .diffusivity_parameters import DiffusivityParameters
from .diffusivity_step import DiffusivityStep
from .tk_diffusivity import TkDiffusivity

from .metadata import metadata

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@molssi.org"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
