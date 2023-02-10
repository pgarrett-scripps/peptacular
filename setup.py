from pathlib import Path

from setuptools import setup

setup(
    name='peptacular',
    version='0.0.3',
    packages=['peptacular'],
    url='',
    license='',
    author='Patrick Garrett',
    author_email='pgarrett@scripps.edu',
    description='Utility package for handling peptide and protein sequences/files',
    install_requires=['regex==2022.10.31'],
    python_requires='>=3.6'
)
