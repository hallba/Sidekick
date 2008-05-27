#!/usr/bin/python
#MDP generator

#These scripts take an mdp file (gromacs input), and change one (or more ) of the run parameters, producing a new file

def replace_seed(filename,output_filename,seed=5):
	#seed has a default here just as a safety precaution- the input must have a seed
	#guaranteed to be random; selected by roll of fair die ;) http://xkcd.com/221/
	output_file = open(output_filename,"w")
	for line in open(filename):
		if line[:8] != "gen_seed": print >>output_file, line,
		else: print >> output_file, "gen_seed                 = " + str(seed),
	output_file.close()