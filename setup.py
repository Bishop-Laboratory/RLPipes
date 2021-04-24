#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name='RSeq',
    version="0.0.1",  # noqa: F821
    scripts=['bin/RSeq'],
    packages=find_packages(),
    include_package_data=True
)
