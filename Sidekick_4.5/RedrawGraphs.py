#!/usr/bin/python

import AnalyseTrajectory

sequence_file = open("sequence.seq","r")
sequence_file_contents = sequence_file.readlines()
sequence = sequence_file_contents[1]

charge = 0
for line in open("ionisation_state.log"):
	if line[:36] == "  System has non-zero total charge: ": charge = float(line[36:]) 

AnalyseTrajectory.visualise(sequence,charge)
