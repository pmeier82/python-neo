# -*- coding: utf-8 -*-

from setuptools import setup
import os

version = open('./neo/version.py', 'r').read().split('\'')[1]
long_description = open("README.txt").read()

setup(
    name = "neo",
    version = version,
    packages = ['neo', 'neo.core', 'neo.io', 'neo.test', 'neo.test.io'],
    install_requires=[
                    'numpy>=1.3.0',
                    'quantities>=0.9.0',
                    ],
    author = "Neo authors and contributors",
    author_email = "sgarcia at olfac.univ-lyon1.fr",
    description = "Neo is a package for representing electrophysiology data in Python, together with support for reading a wide range of neurophysiology file formats",
    long_description = long_description,
    license = "BSD",
    url='http://neuralensemble.org/neo',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering']
)



