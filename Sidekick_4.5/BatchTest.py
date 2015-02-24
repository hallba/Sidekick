#!/usr/bin/python2.5

import os,sys,time

#change to path to give system installed versions of software preference in importing
sys.path = sys.path[1:] + [sys.path[0]]

import HAConf
#import AnalyseTrajectory
import AnalyseSimulation
import GromacsInterface
import GenerateCGSystem
import shutil,time
import BatchTools
from random import randint

#Get options

import SidekickOptParse
#options = SidekickOptParse.get_options(("-s",))
#parser = SidekickOptParse.Helix_Options(("-s",))
parser = SidekickOptParse.Helix_Options(("-s",))
options = parser.options

#Read options

initial_sequence=options.sequence
seed = options.seed
MDsteps = options.simlength*1000./0.02
model_type = options.model_type

#Super option "preformed insertion" resets several other options

if options.preformed_insertion:
	options.bias = True
	options.angle = 90
	options.position = 4.

if options.batch:
	batch_entry = BatchTools.helix_batch_enter(options,networked_tmp_location=False)
	import atexit
	from signal import signal, SIGTERM
	atexit.register(BatchTools.cleanup,batch_entry.temp_folder_name,batch_entry.original_location)
	signal(SIGTERM, lambda signum, stack_frame: exit(1))

time.sleep(20)

#Finally, if in batch mode, create a results directory and move all files there
if options.batch:
	batch_exit = BatchTools.helix_batch_exit(options,temp_folder_location=batch_entry.temp_folder_name,original_folder_location=batch_entry.original_location)
	
