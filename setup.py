#!/usr/bin/env python
# Install script for some random utils.
from setuptools import setup

setup(
    name = 'glm-helperfuncs',
    version = '0.1.0',
    author = 'Max Melin',
    author_email = 'mmelin@g.ucla.edu',
    description = 'Helper functions for training glmhmm.',
    packages = ['glm_hmm_utils'],
    install_requires = [
        'numpy',
        'scipy',
        'matplotlib',
        'pandas',
        'statsmodels',
        'scikit-learn',
        'ipykernel',
        'jupyter',
        'multiprocess',
        'mat73',
        'NeuroDataTypes @ git+https://github.com/mdmelin/NeuroDataTypes.git@state_manuscript',
        'ssm @ git+https://github.com/lindermanlab/ssm.git@master',
    ],
)