#!/usr/bin/env python3

'''Check if course required python version and packages are installed.'''

# 2018-01-30 David.Trudgian@UTSouthwestern.edu
# Quick and dirty script to check for course-required
# python version and packages.

import logging
import sys

logger = logging.getLogger("check_versions")
logger.setLevel(logging.INFO)
logger_stream = logging.StreamHandler()
formatter = logging.Formatter("[%(levelname)6s] %(message)s")
logger_stream.setFormatter(formatter)
logger.addHandler(logger_stream)

REQUIRED_PYTHON = '3.6.4'

REQUIRED_PACKAGES = {
    'numpy': '1.14.0',
    'scipy': '1.0.0',
    'pandas': '0.22.0',
    'matplotlib': '2.1.2',
    'seaborn': '0.8.1',
    'bokeh': '0.12.13',
    'spyder': '3.2.6',
}


def check_package_version(package_name, required_version):
    '''Check package is available and matches desired version.'''

    logger.info("Checking for %s %s", package_name, required_version)

    try:
        mod = __import__(package_name)
    except ImportError:
        logger.error("Package %s not available", package_name)
        return False

    if mod.__version__ != required_version:
        logger.error("%s should be %s, found %s",
                     package_name, required_version, mod.__version__)
        return False

    return True


def main():

    environment_ok = True

    logger.info(
        "Hello - we're checking if your system is ready for the Python 1 Nanocourse")

    # Check Python Version
    if sys.version[:5] != REQUIRED_PYTHON:
        logger.error("Python should be %s, found %s",
                     REQUIRED_PYTHON, sys.version.splitlines()[0])
    else:
        logger.info("Python version OK!")

    # Check Package Versions
    for package_name, required_version in REQUIRED_PACKAGES.items():
        if not check_package_version(package_name, required_version):
            environment_ok = False

    if environment_ok:
        print("\nWoo! - Ready to go, see you at the nanocourse :-)\n")
    else:
        print("\nUh-oh! - Your environment isn't quite right")
        print("Please check you followed the setup instructions and contact the")
        print("course co-ordinator if you continue to have problems.\n")


if __name__ == "__main__":
    main()
