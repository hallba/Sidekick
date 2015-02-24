#!/usr/bin/python
#MDP generator

from __future__ import with_statement

import shutil
#These scripts take an mdp file (gromacs input), and change one (or more ) of the run parameters, producing a new file

def replace_seed(filename,output_filename,seed=5):
	#seed has a default here just as a safety precaution- the input must have a seed
	#guaranteed to be random; selected by roll of fair die ;) http://xkcd.com/221/
	output_file = open(output_filename,"w")
	for line in open(filename):
		if line[:8] != "gen_seed": print >>output_file, line,
		else: print >> output_file, "gen_seed                 = " + str(seed)
	output_file.close()
	
def replace_steps(filename,output_filename,steps=100):
	output_file = open(output_filename,"w")
	for line in open(filename):
		if line[:6] != "nsteps": print >>output_file, line,
		else: print >> output_file, "nsteps                   = " + str(steps)
	output_file.close()

def replace_temperature(filename,output_filename,temp=None):
	if temp == None:
		print "No changes to be made for temperature"
		shutil.copy(filename,output_filename)
	else:
		output_file = open(output_filename,"w")
		for line in open(filename):
			if line[:8] == "gen_temp": print >> output_file, "gen_temp                   = " + str(temp)
			elif line[:5] == "ref_t": print >> output_file,  "ref_t                      = " + str(temp) + " " + str(temp)
			else: print >>output_file, line,
		output_file.close()
	
def replace_epsilon(filename,output_filename,epsilon=20):
	output_file = open(output_filename,"w")
	for line in open(filename):
		if line[:9] != "epsilon_r": print >>output_file, line,
		else: print >> output_file, "epsilon_r                = " + str(epsilon)
	output_file.close()

def replace_pcoupling(filename,output_filename,coupling="semiisotropic"):
	output_file = open(output_filename,"w")
	coupling_types =	{	"semiisotropic"	:	"""Pcoupltype               = semiisotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
tau_p                    = 10 10
compressibility          = 3e-5  3e-5  
ref_p                    = 1.0 1.0""",
							"anisotropic"	:	"""Pcoupltype               = anisotropic
; Time constant (ps), compressibility (1/bar) and reference P (bar) = 
tau_p                    = 10 10 10
compressibility          = 3e-5  3e-5 3e-5 0 0 0
ref_p                    = 1.0 1.0 1.0 1.0 1.0 1.0"""
						}
	
	for line in open(filename):
		if line[:5] == "tau_p" or  line[:5] == "ref_p" or line[:5] == "compr": pass 
		elif line[:10] != "Pcoupltype": print >>output_file, line,
		else: print >> output_file, coupling_types[coupling]
	output_file.close()

def pmartini_update(filename):
	print "apolar-"+filename
	shutil.move(filename,"apolar-"+filename)
	with open(filename,"w") as oup:
		with open("apolar-"+filename) as inp:
			for line in inp:
				if "epsilon_r" in line:
					print >> oup, "epsilon_r                = 2.5"
				elif "coulombtype" in line:
					print >> oup, "coulombtype              = PME"
				elif "rcoulomb" in line and "switch" not in line:
					print >> oup, "rcoulomb                 = 1.3"
				else:
					print >>oup, line,

def pme_update(filename):
	print "apolar-"+filename
	shutil.move(filename,"apolar-"+filename)
	with open(filename,"w") as oup:
		with open("apolar-"+filename) as inp:
			for line in inp:
				if "coulombtype" in line:
					print >> oup, "coulombtype              = PME"
				elif "rcoulomb" in line and "switch" not in line:
					print >> oup, "rcoulomb                 = 1.3"
				else:
					print >>oup, line,