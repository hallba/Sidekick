#!/usr/bin/python

'''Atomistic Restarts for Sidekick'''

from __future__ import with_statement

import os,sys

import HAConf
#import AnalyseTrajectory
import AnalyseSimulation
import GromacsInterface
import GenerateDimerCGSystem
import shutil,time
import BatchTools
from random import randint
from signal import signal, SIGTERM