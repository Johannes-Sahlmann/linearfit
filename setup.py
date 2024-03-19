#!/usr/bin/env python

import os
import sys
import re

try:
    import setuptools
except ImportError:
    from distribute_setup import use_setuptools
    use_setuptools()

from setuptools import find_packages

try:
  from setuptools import setup
  setup
except ImportError:
  from distutils.core import setup
  setup

linearfitVersion = '1.0.2'

setup(
    name="linearfit",
    description="python package that implements a general least-squares fit of a linear model using numpy matrix inversion",
    version=linearfitVersion,
    author="Johannes Sahlmann",
    author_email="jsahlmann@stsci.edu",
    url="https://github.com/Johannes-Sahlmann/linearfit",
    license="BSD 3-Clause",
    long_description="python package that implements a general least-squares fit of a linear model using numpy matrix inversion",
    packages = find_packages(),
    scripts=['test_linearfit.py'],  # this will be installed to a bin/ directory
    package_data={'': ['LICENSE', 'AUTHORS.rst', 'HISTORY.rst', 'INSTALL', 'MANIFEST.in','README.md','README.rst']},
    include_package_data=True,
    install_requires=["numpy"],
    classifiers=[
      "Development Status :: 5 - Production/Stable",
      "Intended Audience :: Developers",
      "Intended Audience :: Science/Research",
      "License :: OSI Approved :: BSD License",
      "Operating System :: OS Independent",
      "Programming Language :: Python",
      "Topic :: Scientific/Engineering :: Astronomy",
      ],
    )
