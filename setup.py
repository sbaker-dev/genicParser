# Copyright (C) 2020 Samuel Baker

DESCRIPTION = "Parse Genetic files via bgi like sql indexing"
LONG_DESCRIPTION = """
# PyGenicParser

Parse Genetic files via bgi like sql indexing

Designed to use sql database style of metadata created by bgenix .bgi files for indexing, but generalised to both plink
and bgen file formats. This was a generalisation of the approach used by the [pybgen], with an attempt to include some
of the features of [pysnptools] such as indexing, although this system isn't as clever in terms of reads.

All the source code can be found at the [PyGenicParser git repository](https://github.com/sbaker-dev/pyGenicParser)

"""
LONG_DESCRIPTION_CONTENT_TYPE = "text/markdown"

DISTNAME = 'pyGenicParser'
MAINTAINER = 'Samuel Baker'
MAINTAINER_EMAIL = 'samuelbaker.researcher@gmail.com'
LICENSE = 'MIT'
DOWNLOAD_URL = "https://github.com/sbaker-dev/pyGenicParser"
VERSION = "0.02.2"
PYTHON_REQUIRES = ">=3.7"

INSTALL_REQUIRES = [
    'numpy',
    'zstd',
]

CLASSIFIERS = [
    'Programming Language :: Python :: 3.7',
    'License :: OSI Approved :: MIT License',
]

if __name__ == "__main__":

    from setuptools import setup, find_packages

    import sys

    if sys.version_info[:2] < (3, 7):
        raise RuntimeError("pyGenicParser requires python >= 3.7.")

    setup(
        name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
        license=LICENSE,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        python_requires=PYTHON_REQUIRES,
        install_requires=INSTALL_REQUIRES,
        packages=find_packages(),
        classifiers=CLASSIFIERS
    )
