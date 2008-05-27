#/usr/bin/python

import sys,os
import HAConf

initial_sequence=sys.argv[1]
scan_residues=sys.argv[2]
starting_point = os.getcwd()

for i in range(len(initial_sequence)):
	if i+len(scan_residues) > len(initial_sequence): continue
	mutant = initial_sequence[:i] + scan_residues + initial_sequence[(i+len(scan_residues)):]
	#mutant = mutant[:len(initial_sequence)]
	print mutant
	os.mkdir(mutant)
	os.chdir(mutant)
	
	#Now create run directories and run
	os.mkdir("Rotation-Translation")
	os.chdir("Rotation-Translation")
	os.system(HAConf.queue_initiator + HAConf.programs['hippo_tr'] + " " + mutant)
	os.chdir("..")
	
	os.mkdir("MARTINI")
	os.chdir("MARTINI")
	os.system(HAConf.queue_initiator + HAConf.programs['MARTINI'] + " " + mutant)
	os.chdir("..")
	
	#Finally, return to the starting directory
	os.chdir(starting_point)