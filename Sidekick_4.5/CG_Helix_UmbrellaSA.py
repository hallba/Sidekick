#!/usr/bin/python2.5

#Robin

import os,sys

#change to path to give system installed versions of software preference in importing
sys.path = sys.path[1:] + [sys.path[0]]

import HAConf
#import AnalyseTrajectory
import AnalyseSimulation
import GromacsInterface
import GenerateCGSystem
import shutil,time
import BatchTools
import PMFTools
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

#Super option "preformed insertion" ignored

if options.preformed_insertion:
	print "Cannot run preformed bilayer for this type of simulation. Exiting"
	sys.exit()

if options.batch:
	batch_entry = BatchTools.helix_batch_enter(options,networked_tmp_location=False)
	import atexit
	from signal import signal, SIGTERM
	atexit.register(BatchTools.cleanup,batch_entry.temp_folder_name,batch_entry.original_location)
	signal(SIGTERM, lambda signum, stack_frame: exit(1))

#Do I want to specify a seed for the lipid positions?
if options.random_lipids:
	lipid_seed = seed
else:
	lipid_seed=None

#Build system, EM (by SD), MD

#Run checks for sensible options
if model_type not in GenerateCGSystem.installed_models:
	print "Model not available at present. Please restart with one of the following options for \"-t\""
	for model in GenerateCGSystem.installed_models.keys():
		print "\t"+model
	exit()

if options.system_size not in GenerateCGSystem.defined_system_sizes.keys():
	print "Size not available at present. Please restart with one of the following options for \"--system-size\""
	for model in GenerateCGSystem.defined_system_sizes.keys():
		print "\t"+model, "~" + str(GenerateCGSystem.defined_system_sizes[model]["lipid"]*64), "lipids"
	exit()
		
#Select a model

try:
	BuildSystem = GenerateCGSystem.installed_models[model_type]
except:
	print "Broken option- model installed but not in lookup list. Exiting"
	sys.exit()

System = BuildSystem(initial_sequence,bias=options.bias,seed=lipid_seed,angle=options.angle,position=options.position,system_size=options.system_size,preformed_insertion=options.preformed_insertion,special=options.special,lipid_type=options.headgroup_mix)

GromacsInterface.SDMinimise(topology_file="ioned_topol",input_structure="ioned.pdb",output_name="em")

print "System set up and minimised"

PMFTools.UmbrellaGenerate(name="umbrella",window=options.position)

if options.bias:
	pcoupling = "semiisotropic"
else:
	pcoupling = "anisotropic"

try:
	GromacsInterface.RobustMolecularDynamics(topology_file="ioned_topol",input_structure="em",output_name="t_0",seed=seed,steps=MDsteps,pcoupling=pcoupling,wallclock=options.wallclock,epsilon=System.epsilon,temperature=options.temperature,apolar=System.apolar,pme=System.pme,umbrella_name="umbrella")
except:
	if options.batch:
		sys.exit()

#analyse- this is optional if in debug mode

if HAConf.debug_status:
	try:
		#AnalyseTrajectory.visualise(initial_sequence,System.charge,bias=options.bias,number_lipid_types=number_lipid_types)
		Analysis = AnalyseSimulation.AnalyseHelix(System)
	except:
		pass
else:
	#AnalyseTrajectory.visualise(initial_sequence,System.charge,bias=options.bias,number_lipid_types=number_lipid_types)
	Analysis = AnalyseSimulation.AnalyseHelix(System)

#Cleanup code
##Bzips the tachyon input files
##Removes all trr files (if everything went well, we shouldn't need them...)

os.system("bzip2 *0.dat")
os.system("rm *0.dat")
os.system("rm *.trr")

#Finally, if in batch mode, create a results directory and move all files there
if options.batch:
	batch_exit = BatchTools.helix_batch_exit(options,temp_folder_location=batch_entry.temp_folder_name,original_folder_location=batch_entry.original_location)
	
