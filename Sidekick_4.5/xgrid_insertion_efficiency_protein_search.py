#!/usr/bin/python

import sys,os
import HAConf
from random import randint

from xgrid_tools import *

starting_point = os.getcwd()
input_file = sys.argv[1]
number_of_repeats = 10

job_list =[]

seed_packet = [randint(0,100000) for i in range(number_of_repeats)]

#Tree structure is different here
#Type/Sequence/Repeats

os.mkdir("Hippo_Rotation-Translation")
os.chdir("Hippo_Rotation-Translation")

#sequences = [line[:-1] for line in open(input_file,'r')]

#Read the sequence of the complete protein
fasta_fh = open(input_file,'r')
fasta_sequence = ""

for line in fasta_fh:
	if line[0] == ">" : continue
	fasta_sequence += line[:-1]

fasta_fh.close()

residue_efficiencies = [[0.0,0] for residue in fasta_sequence]
window_size = 30
sequences = []

for i in range(len(fasta_sequence) - window_size + 1):
	sequences.append(fasta_sequence[i:i+window_size]

for mutant in sequences:
        print mutant
        try:
                os.mkdir(mutant)
        except:
                continue
        os.chdir(mutant)
        
        #Now create run directories and run
        job_list.append(job_submit(HAConf.programs['hippo_tr'] + " " + mutant + " -dAgg "))
        os.chdir("..")
        
os.chdir(starting_point)

os.mkdir("MARTINI")
os.chdir("MARTINI")

for mutant in sequences:
        print mutant
        try:
                os.mkdir(mutant)
        except:
                continue
        os.chdir(mutant)
        
        for random_seed in seed_packet:
        	os.mkdir(str(random_seed))
        	os.chdir(str(random_seed))
        	job_list.append(job_submit(HAConf.programs['Insertion_Efficiency'] + " " + mutant + " " + str(random_seed) + " -dAgg "))
        	os.chdir("..")
        
        #Finally, return to the starting directory
        os.chdir("..")

os.chdir(starting_point)

#Hang around and collect the results as a daemon

collect_results_daemon(job_list,starting_point+"/daemon.log",starting_point+"/restart.pickle")

#Scrape off the calculated insertion efficiencies from each MARTINI run and produce a net efficiency for each sequence

martini_location = starting_point + "/MARTINI"
summary_table = open(starting_point + "/MARTINI_Insertion_Summary.txt", "w")

position = 0

for mutant in sequences:
	sequence_location = martini_location + "/" + mutant
	efficiency = 0
	for random_seed in seed_packet:
		efficiency_file = sequence_location + "/" + str(seed) + "/efficiency.dat"
		efficiency += float(efficiency_file.readline())
		efficiency_file.close()
	efficiency /= len(seed_packet)
	for i in range(position,position+window_size):
		residue_efficiencies[i][0] += efficiency
		residue_efficiencies[i][1] += 1
	position += 1
	print >> summary_table, mutant, efficiency

summary_table.close()

residue_efficiency_file = open("PerResidueEfficiency.txt","w")
for value_pair in residue_efficiencies:
	per_res_value = value_pair[0]/value_pair[1]
	print >> residue_efficiency_file, per_res_value

residue_efficiency_file.close()
