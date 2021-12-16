import sys

from tbpore.constants import *  # noqa: F401,F403

"""
Version has unique source in pyproject.toml.
importlib fetches version from distribution metadata files
(in dist-info or egg-info dirs).
From Python 3.8, importlib_metadata is in standard library as importlib.metadata.
"""
if sys.version_info >= (3, 8):
    from importlib import metadata
else:
    import importlib_metadata as metadata

__version__ = metadata.version("tbpore")
