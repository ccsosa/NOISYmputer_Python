#!/usr/bin/env python3

"""
Parts of this file were taken from the pyzmq project
(https://github.com/zeromq/pyzmq) which have been permitted for use under the
BSD license. Parts are from lxml (https://github.com/lxml/lxml)
"""

import numpy as np
import pandas as pd
import multiprocessing
import os
import sys
import cython
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from tqdm import tqdm

#versioning
import versioneer

cmdclass = versioneer.get_cmdclass()


# note: sync with pyproject.toml
#https://uoftcoders.github.io/studyGroup/lessons/python/packages/lesson/
#https://github.com/montoyamoraga/tutorial-python-module
#https://stackoverrun.com/es/q/4141323
setup(
    # Needed to silence warnings (and to be a worthwhile package)
    name='NOISYmputer',
    url='https://github.com/ccsosa/NOISYmputer_Python',
    author='Mathias Lorieux',
    author_email='jm.lorieux@gmail.com ',
    # Needed to actually package something
    packages=['NOISYmputer'],
    # Needed for dependencies
    install_requires=['numpy','pandas','os','sys','cython'],
    # *strongly* suggested for sharing
    version='0.1',
    # The license can be anything you like
    license='MIT',
    description='An example of a python package from pre-existing code',
    # We will also need a readme eventually (there will be a warning)
    # long_description=open('README.txt').read(),
)