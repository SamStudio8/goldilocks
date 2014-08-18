#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools

README = open('README.rst').read()
CHANGELOG = open('CHANGELOG.rst').read()

requirements = [
    "numpy"
]

test_requirements = [
    "tox",
    "pytest",
    "nose",
    "python-coveralls",
]

setuptools.setup(
    name="goldilocks",
    version="0.0.2",
    url="https://github.com/samstudio8/goldilocks",

    description="Locating genomic regions that are \"just right\".",
    long_description=README + '\n\n' + CHANGELOG,

    author="Sam Nicholls",
    author_email="sam@samnicholls.net",

    maintainer="Sam Nicholls",
    maintainer_email="sam@samnicholls.net",

    packages=setuptools.find_packages(),
    include_package_data=True,

    install_requires=requirements,

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    test_suite='tests',
    tests_require=test_requirements
)
