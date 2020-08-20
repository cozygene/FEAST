# ----------------------------------------------------------------------------
# Copyright (c) 2019, FEAST development team.
#
# Distributed under the terms of the Modified cc-by-sa-4.0 License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
# Note that this file doesn't declare any install_requires packages -- this is
# because q2-FEAST assumes that it's being installed into a "normal" QIIME 2
# conda environment, and thus that all of the following
# (not-in-the-standard-python-library) packages will be available:
# - pandas
# - numpy
# - biom
# - qiime2
# - q2_types
# - q2_taxa
# ... However, this file does declare some packages under extra_requires, which
# should only be needed when running tests.

import versioneer
from setuptools import setup, find_packages

setup(
    name="feast",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    author="Liat Shenhav",
    author_email="liashenhav@gmail.com",
    description="Fast Expectation-mAximization microbial Source Tracking (FEAST)",
    license='cc-by-sa-4.0',
    url="https://github.com/cozygene/FEAST",
    # idiom based on https://github.com/altair-viz/altair/blob/master/setup.py
    extras_require={
        "dev": [
            "nose",
            "coverage"
        ]
    },
    entry_points={
        'qiime2.plugins': ['feast=q2_feast.plugin_setup:plugin']
    },
    scripts=['q2_feast/assets/source_tracking.R'],
    package_data={
        "q2_feast": ['citations.bib'],
    },
    zip_safe=False,
)
