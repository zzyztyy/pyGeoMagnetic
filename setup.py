#!/usr/bin/python3
# -*- coding: utf-8 -*-

from __future__ import print_function
from setuptools import setup, find_packages


setup(
    name="pyGeoMagApex",
    version="0.1.0",
    author="zzyztyy",
    author_email="2375672032@qq.com",
    description="Coordination transform between GeoGraphic and GeoMagnetic",
    long_description=open("README.md").read(),
    license="MIT",
    url="https://github.com/zzyztyy/pyGeoMagnetic",
    packages=['pyGeoMagApex'],
    install_requires=[
        "NumPy",
        "pyIGRF"
    ]
)
