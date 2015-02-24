#!/usr/bin/python

import sys,os
import HAConf
import xgrid_tools
import SidekickOptParse

#This script calculates a PMF from a fasta sequence using the Sidekick tools, by creating a batchfile which runs CG_PMG_window.py

parser = SidekickOptParse.Multi_PMF_Helix_Options(("-s","-j","-w"))
options = parser.options

from xgrid_tools import *
starting_point = os.getcwd()

def write_batch_task(batch_file,mutant,model_type,window=0,bias=True,simlength=200,jobname=""):
			mutant_name = mutant + "_" + model_type + "_" + str(window)
			print >> batch_file, "\t"*3 + mutant_name + " = {"
			print >> batch_file, "\t"*4 + "arguments = ("
			print >> batch_file, "\t"*5 + '"-s",'
			print >> batch_file, "\t"*5 + '"' + mutant + '",'
			print >> batch_file, "\t"*5 + '"-t",'
			print >> batch_file, "\t"*5 + '"' + model_type + '",'
			print >> batch_file, "\t"*5 + '"-w",'
			print >> batch_file, "\t"*5 + '"' + str(window) + '",'
			#print >> batch_file, "\t"*5 + '"-l",'
			#print >> batch_file, "\t"*5 + '"' + str(simlength) + '",'
			print >> batch_file, "\t"*5 + '"-j",'
			print >> batch_file, "\t"*5 + '"' + jobname + '",'
			if not bias:
				print >> batch_file, "\t"*5 + '"-u",'
			print >> batch_file, "\t"*5 + '"-b"'
			print >> batch_file, "\t"*4 + ");"
			print >> batch_file, "\t"*4 + 'command = "' + HAConf.programs['PMF_CG_Helix'][:-1] + '";'
			#print >> batch_file, "\t"*4 + 'command = "/Users/Shared/Sidekick/CG_Helix.py";'
			print >> batch_file, "\t"*3 + "};"


if options.batch:
	print "This program always runs batch jobs by default. The batch option (-b) is unnecessary"
#if options.duplicate != 1:
	#print "This program will only run a single seed for each window at present, to minimise load"
	#options.duplicate = 1
	
#Determine if this is a scan, or input from a file

if options.sequence_file == None and options.sequence == None:
	print "This program requires input in terms of a file of fasta sequences (-f) or a sequence (-s), or mutate (-m)"
	print "Try -h as an option for more info"
	exit()
elif options.sequence_file != None and options.sequence != None:
	print "This only works on a sequence, or a file of sequences"
	print "Try -h as an option for more info"
	exit()

if options.batch:
	print "This program always runs batch jobs by default. The batch option (-b) is unnecessary"

if options.window:
	print "Running " + str(int(options.window)) + " 1A windows around the bilayer center"


#Models or model scan?

if options.model_type not in HAConf.installed_models:
	if options.model_type == "search":
		models = HAConf.installed_models
	else:
		print "Choose an installed model from the list below or 'search' to run all models"
		for item in HAConf.installed_models:
			print '\t' +  item
		print "Try -h as an option for more info"
		exit()
else:
	models = [options.model_type,]
	
if options.sequence_file:
	sequences = [line[:-1] for line in open(options.sequence_file,'r')]
elif options.mutant:
	sequences = []
	for i in range(len(options.sequence)):
		scan_residues = options.mutant
		if i+len(scan_residues) > len(options.sequence): continue
		mutant = options.sequence[:i] + scan_residues + options.sequence[(i+len(scan_residues)):]
		if mutant not in sequences:
			sequences.append(mutant)
else:
	sequences = [options.sequence,]
	
window_frame = [(window - int(options.window)) for window in range(int(options.window)*2)]

#Generate batch file

batch_file = open(options.destination + ".batch","w")

##Write Header

print >> batch_file, """{
	jobSpecification = {
		name = \"""" + options.destination + """\";
		taskSpecifications = {"""

for model_type in models:
	for mutant in sequences:
		for window in window_frame:
			write_batch_task(batch_file,mutant,model_type,window=window,bias=options.bias,simlength=options.simlength,jobname=options.destination)
	
print >> batch_file, "\t"*2 + "};"
print >> batch_file, "\t" + "};"
print >> batch_file, "}"

batch_file.close()

#Submit batch job
job_list = []
job_list.append(batch_submit(options.destination+".batch"))