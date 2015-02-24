#/usr/bin/python

import sys,os
import HAConf

from xgrid_tools import *

initial_sequence=sys.argv[1]
scan_residues=sys.argv[2]
starting_point = os.getcwd()

submission_list = []
job_list =[]

for i in range(len(initial_sequence)):
	if i+len(scan_residues) > len(initial_sequence): continue
	mutant = initial_sequence[:i] + scan_residues + initial_sequence[(i+len(scan_residues)):]
	#mutant = mutant[:len(initial_sequence)]
	print mutant
	try:
		os.mkdir(mutant)
	except:
		continue
	os.chdir(mutant)
	
	#Now create run directories and run
	os.mkdir("Rotation-Translation")
	os.chdir("Rotation-Translation")
	job_list.append(job_submit(HAConf.programs['hippo_tr'] + " " + mutant + " -dAgg "))
	os.chdir("..")
	
	os.mkdir("MARTINI")
	os.chdir("MARTINI")
	job_list.append(job_submit(HAConf.programs['MARTINI'] + " " + mutant + " " + str(5) + " -dAgg "))
	os.chdir("..")
	
	#Finally, return to the starting directory
	os.chdir(starting_point)
	
collect_results_daemon(job_list,starting_point+"/daemon.log",starting_point+"/restart.pickle")

