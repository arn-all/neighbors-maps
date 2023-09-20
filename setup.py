from setuptools import setup

setup(
    name='neighbors_map',
    version='0.1.0',
    description='Python module for computing the Neighbors Maps of atoms.',
    packages=['neighbors_map'],
    install_requires=[
        'numpy',
        'scipy',
        'ase',
        'matplotlib'
    ],
)