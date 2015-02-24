#!/usr/bin/python2.5

#Kato

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

#Get options

import SidekickOptParse
#options = SidekickOptParse.get_options(("-s",))
#parser = SidekickOptParse.Helix_Options(("-s",))
parser = SidekickOptParse.Helix_Dimer_Options(("--s1","--s2"))
options = parser.options

#Read options

#HelixA and HelixB are always stored in alphabetical order
sequences=[options.sequence1,options.sequence2]
sequences.sort()
if sequences[0] != options.sequence1:
	sequences = [options.sequence2,options.sequence1]
	[options.sequence2,options.sequence1] = [options.sequence1,options.sequence2]

seed = options.seed
MDsteps = options.simlength*1000./0.02
model_type = options.model_type

#Run checks for sensible options
if model_type not in GenerateDimerCGSystem.installed_models:
	print "Model not available at present. Please restart with one of the following options for \"-t\""
	for model in GenerateDimerCGSystem.installed_models.keys():
		print "\t"+model
	exit()
		
#Select a model

try:
	BuildSystem = GenerateDimerCGSystem.installed_models[model_type]
except:
	print "Broken option- model installed but not in lookup list. Exiting"
	sys.exit()

if options.batch:
	batch_entry = BatchTools.helix_dimer_batch_enter(options)
	import atexit
	atexit.register(BatchTools.cleanup,batch_entry.temp_folder_name,batch_entry.original_location)
	signal(SIGTERM, lambda signum, stack_frame: exit(1))

System = BuildSystem(options.sequence1,options.sequence2,parallel=options.parallel,special=options.special,rotatebyangle=options.rotatebyangle)

GromacsInterface.SDMinimise(topology_file="ioned_topol",input_structure="ioned.pdb",output_name="em")

pcoupling = "semiisotropic"

try:
	GromacsInterface.RobustMolecularDynamics(topology_file="ioned_topol",input_structure="em",output_name="t_0",seed=seed,steps=MDsteps,pcoupling=pcoupling,wallclock=options.wallclock*2,epsilon=System.epsilon,temperature=options.temperature)
except:
	if options.batch:
		#BatchTools.unclean_exit_dimer(batch_entry.temp_folder_name)
		sys.exit()

if HAConf.debug_status:
	try:
		#AnalyseTrajectory.visualise(initial_sequence,System.charge,bias=options.bias,number_lipid_types=number_lipid_types)
		Analysis = AnalyseSimulation.AnalyseHelixDimer(System)
	except:
		pass
else:
	#AnalyseTrajectory.visualise(initial_sequence,System.charge,bias=options.bias,number_lipid_types=number_lipid_types)
	Analysis = AnalyseSimulation.AnalyseHelixDimer(System)

#Clean up

os.system("rm *.trr")

#Finally, if in batch mode, create a results directory and move all files there
if options.batch:
	batch_exit = BatchTools.helix_dimer_batch_exit(options,temp_folder_location=batch_entry.temp_folder_name,original_folder_location=batch_entry.original_location)