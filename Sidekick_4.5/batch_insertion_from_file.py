#!/usr/bin/python

import sys,os
import HAConf
from random import randint

from xgrid_tools import *
starting_point = os.getcwd()

#Get options

import SidekickOptParse
options = SidekickOptParse.get_options(("-f",))

#Read options

input_file = options.sequence_file
number_of_repeats = options.duplicates

job_list =[]
seed_packet = [randint(0,100000) for i in range(number_of_repeats)]

#Get sequences from file

sequences = [line[:-1] for line in open(input_file,'r')]

#Generate batch file

batch_file = open("insertion.batch","w")

##Write Header

print >> batch_file, """{
	jobSpecification = {
		name = "Insertion Efficiency";
		taskSpecifications = {"""

for mutant in sequences:
	for seed in seed_packet:
			mutant_name = mutant + "_" + str(seed)
			print >> batch_file, "\t"*3 + mutant_name + " = {"
			print >> batch_file, "\t"*4 + "arguments = ("
			print >> batch_file, "\t"*5 + '"-s",'
			print >> batch_file, "\t"*5 + '"' + mutant + '",'
			print >> batch_file, "\t"*5 + '"-r",'
			print >> batch_file, "\t"*5 + '"' + str(seed) + '",'
			print >> batch_file, "\t"*5 + '"-l",'
			print >> batch_file, "\t"*5 + '"' + str(50) + '",'
			print >> batch_file, "\t"*5 + '"-b"'
			print >> batch_file, "\t"*4 + ");"
			print >> batch_file, "\t"*4 + 'command = "' + HAConf.programs['CG_Helix'][:-1] + '";'
			#print >> batch_file, "\t"*4 + 'command = "/Users/Shared/Sidekick/CG_Helix.py";'
			print >> batch_file, "\t"*3 + "};"

print >> batch_file, "\t"*2 + "};"
print >> batch_file, "\t" + "};"
print >> batch_file, "}"

batch_file.close()

#Submit batch job

job_list.append(batch_submit("insertion.batch"))

#Monitor for job completion

collect_results_daemon(job_list,starting_point+"/daemon.log",starting_point+"/restart.pickle")