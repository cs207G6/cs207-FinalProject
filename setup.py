"""A setuptools based setup module.

See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

from setuptools import setup, find_packages

from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='chemkin',

    version='1.0.0',

    description='CS 207 Group 6 chemkin project',

    long_description=long_description,

    url='https://github.com/cs207G6/cs207-FinalProject',

    author='CS 207 Fall 2017 Group 6',

    license='MIT',

    keywords='cs207 group6 chemkin',

    packages=find_packages(exclude=['data', 'documentation', 'tests']),

    install_requires=['numpy', 'pandas', 'pytest-runner'],

    tests_require=['coverage', 'pytest', 'pytest-cov'],

)
