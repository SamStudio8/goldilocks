#!/usr/bin/env python
# -*- coding: utf-8 -*-

import setuptools

README = open('README.rst').read()
CHANGELOG = open('CHANGELOG.rst').read()

requirements = [
    "numpy",
    "matplotlib"
]

test_requirements = [
    "tox",
    "pytest",
    "nose",
    "python-coveralls",
]

setuptools.setup(
    name="goldilocks",
    version="0.1.1",
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

    entry_points = {
        "console_scripts": ["goldilocks=goldilocks.cmd:main"]
    },

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
    ],

    test_suite='tests',
    tests_require=test_requirements
)
