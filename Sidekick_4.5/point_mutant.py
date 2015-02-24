#/usr/bin/python

import sys,os
import HAConf

amino_converter = {	'ALA': 'A',			'ARG': 'R',			'ASN': 'N',
			'ASP': 'D',			'CYS': 'C',			'GLU': 'E',
			'GLN': 'Q',			'GLY': 'G',			'HIS': 'H',
			'ILE': 'I',			'LEU': 'L',			'LYS': 'K',
			'MET': 'M',			'PHE': 'F',			'PRO': 'P',
			'SER': 'S',			'THR': 'T',			'TRP': 'W',
			'TYR': 'Y',			'VAL': 'V'
			}

initial_sequence=sys.argv[1]
position=sys.argv[2]
starting_point = os.getcwd()
seed = 4


for value in amino_converter.keys():
	mutant = initial_sequence[:(position-1)] + amino_converter[value] + initial_sequence[position:]
	#mutant = mutant[:len(initial_sequence)]
	print mutant
	os.mkdir(mutant)
	os.chdir(mutant)
	
	#Now create run directories and run
	os.mkdir("Rotation-Translation")
	os.chdir("Rotation-Translation")
	os.system(HAConf.queue_initiator + HAConf.programs['hippo_tr'] + " " + mutant + " " + seed)
	os.chdir("..")
	
	os.mkdir("MARTINI")
	os.chdir("MARTINI")
	os.system(HAConf.queue_initiator + HAConf.programs['MARTINI'] + " " + mutant + " " + seed)
	os.chdir("..")
	
	#Finally, return to the starting directory
	os.chdir(starting_point)