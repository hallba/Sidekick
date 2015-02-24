#!/usr/bin/python

'''Restart Manager for Sidekick'''

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

class empty_object:
	pass

def dimer_restart(sequenceA,sequenceB,batch_entry,seed):
	print "Dimer start"
	system = empty_object()
	options = empty_object()
	'''rebuild system/options from topology file'''
	
	system.charge = 0
	system.lipid_type="DPPC" #No options at present
	system.number_lipid_types = 1
	system.bias = True
	system.topology_name = "ioned_topol"
	
	
	with open("ioned_topol.top") as input_topol:
		for line in input_topol:
			if line[:2] == "NA":
				system.charge -= int(line.split()[-1])
			elif line[:2] == "CL":
				system.charge += int(line.split()[-1])
			
			
	if not os.system("grep Aniparallel Sidekick_log.txt"):
		system.parallel = False
	else:
		system.parallel = True
	
	system.initial_sequence_A = sequenceA
	system.initial_sequence_B = sequenceB
	system.initial_sequence = sequenceA
	
	system.topology_file = "ioned_topol"
	
	options.seed = seed
	options.destination = "rescue"
	options.special = None
	options.temperature = None
	
	options.sequence1, options.sequence2 =system.initial_sequence_A, system.initial_sequence_B
	
	options.model_type = "Unknown"

	'''restart simulation'''
	
	print "Restart MD"
	print os.getcwd()
	
	try:
		if not os.path.exists("vector.dat"):
			#If we've started the analysis we shouldn't repeat the last portion of the simulation- just analyse and exit
			GromacsInterface.restartMD()
		
		#topology_file="ioned_topol",input_structure="em",output_name="t_0",seed=seed,steps=MDsteps,pcoupling=pcoupling,wallclock=options.wallclock,epsilon=System.epsilon,temperature=options.temperature
	except:
			print "Something went wrong"
			#if options.batch:
			#BatchTools.unclean_exit_dimer(batch_entry.temp_folder_name)
			sys.exit()
	
	if HAConf.debug_status:
		try:
			Analysis = AnalyseSimulation.AnalyseHelixDimer(system)
		except:
			pass
	else:
		Analysis = AnalyseSimulation.AnalyseHelixDimer(system)
	
	#Clean up
	
	os.system("rm *.trr")
	
	#Finally, if in batch mode, create a results directory and move all files there
	#if options.batch:
	batch_exit = BatchTools.helix_dimer_batch_exit(options,temp_folder_location=batch_entry.temp_folder_name,original_folder_location=batch_entry.original_location)


'''Input directory'''

'''
Sidekick_temp_998188701_glados-node-08_DILVVLLSVMGAILLIGLAALLIWKLLITIH_________-IPIWWVLVSVLGGLLLLTILVLAMWKV______________000065465
Sidekick_temp_998283769_glados-node-11_GKALWLALALALALALALWLAKA_000046241
'''

class fake_batch:
	def __init__(self,tmpdir,resultsdir):
		self.temp_folder_name = tmpdir
		self.original_location = resultsdir

print "Repeating job " + sys.argv[1]

job_directory = HAConf.results_location + "/tmp/" + sys.argv[1]
contents = sys.argv[1].replace("_"," ").split()

if not os.path.exists(job_directory):
	print "Temporary files not present. Completed"
	sys.exit()

seed = int(contents[-1])

if contents[-2][0] == "-":
	'''dimer'''
	sequenceB = contents[-2][1:]
	sequenceA = contents[-3]
	seed = int(contents[-1])
	
	batch_entry = fake_batch(job_directory,os.getcwd())
	os.chdir(job_directory)
	
	import atexit
	atexit.register(BatchTools.cleanup,batch_entry.temp_folder_name,batch_entry.original_location)
	signal(SIGTERM, lambda signum, stack_frame: exit(1))
	
	dimer_restart(sequenceA,sequenceB,batch_entry,seed)
else:
	'''monomer'''
	pass
