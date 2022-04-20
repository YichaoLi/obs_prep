import os
import warnings
from setuptools import setup, find_packages
from distutils import sysconfig
import numpy as np

REQUIRES = ['numpy', 'scipy', 'matplotlib']

# Don't install requirements if on ReadTheDocs build system.
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    requires = []
else:
    requires = REQUIRES

setup(
    name = 'obs_prep',
    version = '0.1.0',

    packages = find_packages(),
    ext_modules = [],
    install_requires = requires,
    package_data = {},

    # metadata for upload to PyPI
    author = "Yi-Chao LI",
    description = "Utils",
)
