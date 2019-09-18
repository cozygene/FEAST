# ----------------------------------------------------------------------------
# Copyright (c) 2019, FEAST development team.
#
# Distributed under the terms of the Modified cc-by-sa-4.0 License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import versioneer
from setuptools import setup, find_packages

setup(
    name="q2-FEAST",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    packages=find_packages(),
    author="Liat Shenhav",
    author_email="liashenhav@gmail.com",
    description="Fast Expectation-mAximization microbial Source Tracking (FEAST)",
    license='cc-by-sa-4.0',
    url="https://www.bioconductor.org/packages/devel/bioc/html/ALDEx2.html",
    entry_points={
        'qiime2.plugins': ['q2-FEAST=q2_FEAST.plugin_setup:plugin']
    },
    scripts=['q2_FEAST/assets/source_tracking.R'],
    package_data={
        "q2_FEAST": ['citations.bib'],
    },
    zip_safe=False,
)