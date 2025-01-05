#! /usr/bin/env python3

import os
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

package = 'zfc-ebar'
version = '0.1.1'


def readme():
    with open(
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'README.md'
            )
    ) as f:
        return f.read()

long_description=readme()

setup(
    name=package,
    version=version,
    description="ZFC-eBAR is a software tool for analyzing CRISPR-based high-throughput screening data, specifically designed to assess the enrichment of epegRNA/sgRNA from read count data.",
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 1 - Planning',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)'
    ],
    url='https://github.com/UronicAcid/ZFC-eBAR',
    author='Wei Tang',
    author_email='tangwei@stu.pku.edu.cn',
    license='GPL',
    packages=['zfc-ebar'],
    install_requires=[
        'numpy', 'scipy', 'pandas',
        'matplotlib', 'sklearn', 'statsmodels'
    ],
    scripts=[
        'bin/zfc-ebar',
    ],
    package_dir={'zfc': 'zfc-ebar'},
    include_package_data=True,
    zip_safe=False
)
