# -*- coding: utf-8 -*-

from os.path import abspath, dirname, isfile
from sys import path
from setuptools import setup, find_packages
from biggr import __version__

# To temporarily modify sys.path
SETUP_DIR = abspath(dirname(__file__))

requirement_path = f"{SETUP_DIR}/requirements.txt"
install_requires = []
if isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()

setup(
    name="biggr",
    version=__version__,
    description="""Easy access to the BiGGr API.""",
    url="https://github.com/pascalaldo/biggr",
    author="Pascal A. Pieters",
    author_email="pascalaldo@gmail.com",
    license="GPL",
    classifiers=[
        "License :: OSI Approved :: GPL License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
    ],
    keywords="systems biology, genome-scale model",
    packages=find_packages(),
    install_requires=install_requires,
)
