#!/usr/bin/python

import os,sys

from HAConf import configuration
import HAConf
import AnalyseTrajectory
import GromacsInterface
import Generate_PMFCG_System

#Get options

import SidekickOptParse
#options = SidekickOptParse.get_options(("-s",))
parser = SidekickOptParse.PMF_Helix_Options(("-s","-w"))
options = parser.options

#Read options

initial_sequence=options.sequence
seed = options.seed
window = options.window
MDsteps = options.simlength*1000./0.02
model_type = options.model_type

if options.batch:
	#Batch mode creates a directory tree for the simulation to run in
	if options.destination:
		os.chdir(HAConf.results_location)
		try:
			os.mkdir(options.destination)
		except:
			pass
		os.chdir(options.destination)
	#Create a "PMF" directory,if non exists"
	try:
		os.mkdir("PMF")
	except:
		pass
	os.chdir("PMF")	
	#Subtle change to the directory names- include info on biasing in the name to avoid confusion
	try:
		os.mkdir(model_type)
	except:
		pass
	os.chdir(model_type)
	
	try:
		os.mkdir(initial_sequence)
	except:
		pass
	os.chdir(initial_sequence)
	
	os.mkdir("w"+str(window))
	os.chdir("w"+str(window))
	
	#write logs out to the new directory, rather than stdout
	stdout = "Sidekick_log.txt"
	stderr = "Sidekick_err.txt"
	for f in sys.stdout, sys.stderr: f.flush()
	so = file(stdout, 'a+')
	se = file(stderr, 'a+',0)
	os.dup2(so.fileno(),sys.stdout.fileno())
	os.dup2(se.fileno(),sys.stderr.fileno())
	

#Build system, EM (by SD), MD

if model_type not in HAConf.installed_models:
	print "Model not available at present. Please restart with one of the following options for \"-t\""
	for model in HAConf.installed_models:
		print "\t"+model
	exit()
		
if model_type == "MARTINI2.1.1":
	System = Generate_PMFCG_System.MARTINI_system_2_1_1(initial_sequence,window=options.window)
elif model_type == "Bond":
	System = Generate_PMFCG_System.BOND_system(initial_sequence,window=options.window)
elif model_type == "Bond0.9.5":
	System = Generate_PMFCG_System.BOND_system_0_9_5(initial_sequence,window=options.window)
else:
	System = Generate_PMFCG_System.MARTINI_system(initial_sequence,window=options.window)

GromacsInterface.SDMinimise(topology_file="ioned_topol",input_structure="ioned.pdb",output_name="em")

GromacsInterface.PRMolecularDynamics(topology_file="ioned_topol",input_structure="em",output_name="pr")

GromacsInterface.PMFMolecularDynamics(topology_file="ioned_topol",input_structure="pr",output_name="t_0")

#analyse

AnalyseTrajectory.visualise(initial_sequence,System.charge)

#Cleanup code
##Bzips the tachyon input files
##Removes all trr files (if everything went well, we shouldn't need them...)

os.system("bzip2 *0.dat")
os.system("rm *0.dat")
os.system("rm *.trr")
