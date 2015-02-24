#!/usr/bin/python

import os

starting_point = os.getcwd()

for line in os.walk(os.getcwd()) :  
	if line[0][-7:] == "MARTINI":
		os.chdir(line[0])
		os.system("/nfsmount/Sidekick/RedrawGraphs.py -dAgg")
		os.chdir(starting_point)
