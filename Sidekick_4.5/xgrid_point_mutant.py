#/usr/bin/python

import sys,os
import HAConf
from xgrid_tools import *

amino_converter = {	'ALA': 'A',			'ARG': 'R',			'ASN': 'N',
			'ASP': 'D',			'CYS': 'C',			'GLU': 'E',
			'GLN': 'Q',			'GLY': 'G',			'HIS': 'H',
			'ILE': 'I',			'LEU': 'L',			'LYS': 'K',
			'MET': 'M',			'PHE': 'F',			'PRO': 'P',
			'SER': 'S',			'THR': 'T',			'TRP': 'W',
			'TYR': 'Y',			'VAL': 'V'
			}

initial_sequence=sys.argv[1]
position=int(sys.argv[2])
starting_point = os.getcwd()
seed = 5

job_list = []

for value in amino_converter.keys():
	mutant = initial_sequence[:(position-1)] + amino_converter[value] + initial_sequence[position:]
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
	job_list.append(job_submit(HAConf.programs['hippo_tr'] + " " + mutant + " -dAgg &"))
	os.chdir("..")
	
	os.mkdir("MARTINI")
	os.chdir("MARTINI")
	job_list.append(job_submit(HAConf.programs['MARTINI'] + " " + mutant + " " + str(seed) + " -dAgg &" ))
	os.chdir("..")
	
	#Finally, return to the starting directory
	os.chdir(starting_point)

collect_results_daemon(job_list,starting_point+"/daemon.log")
