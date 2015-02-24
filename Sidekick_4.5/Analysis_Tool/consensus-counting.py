#!/usr/bin/python

import os

#This scrapes the "efficiency" data to calculate the percentage insertion over the entire trajectory

def get_efficiency(run_directory):
	'''This returns a limited set of results corresponding to what is found
	[flag, percentage]
	0 = file not found or otherwise abnormal-ignore
	1 = not inserted for > 50% of the time
	2 = inserted for > 50 % of the time
	'''
	try:
		run_output = open(run_directory+"/efficiency.dat")
	except:
		#print "Cannot read file", run_directory+"/efficiency.dat"
		return [0, None]
	for line in run_output:
		if line[:5] == "Numbe":
			contents = line.split()
			if int(contents[-1]) != 1:
				#print "Wrong number of bilayers..."
				return [0, None]
		elif line[:5] == "Unusu":
			#print "Inappropriate mixing"
			return [0, None]
		elif line[:5] == "Helix":
			#print line
			contents = line.split()
			#print contents[-2]
			if float(contents[-2]) > 50:
				return [2, float(contents[-2]) ]
			else:
				return [1, float(contents[-2]) ]
		else:
			print "Error- no case is satisfied"
			exit()

#Where am I? I'll start by assuming I'm in a "results" directory- that is /data_location/job_name 

starting_point = os.getcwd()
results = {}

for forcefield in os.listdir(starting_point):
	for sequence in os.listdir(starting_point + "/" + forcefield):
		#inserted = 0.
		#successful_runs = 0.
		#percentage_insertion = 0.
		#now check for lipids and contents
		for lipid in os.listdir(starting_point + "/" + forcefield + "/" + sequence):
			for mix in os.listdir(starting_point + "/" + forcefield + "/" + sequence + "/" + lipid):
				inserted = 0.
				successful_runs = 0.
				percentage_insertion = 0.
				for run in os.listdir(starting_point + "/" + forcefield + "/" + sequence + "/" + lipid + "/" + mix):
					
			
					run_directory = starting_point + "/" + forcefield + "/" + sequence + "/" + lipid + "/" + mix + "/" + run 
					[efficiency, percent] = get_efficiency(run_directory)
					if efficiency == 0:
						continue
					elif efficiency == 2:
						inserted += 1
					percentage_insertion += percent
					successful_runs += 1
				if successful_runs > 0:
					results[forcefield+"_"+sequence + "_" + lipid + "_" + mix] = (percentage_insertion/successful_runs, inserted/successful_runs, inserted, successful_runs)
				else:
					results[forcefield+"_"+sequence + "_" + lipid + "_" + mix] = ("NA", inserted, successful_runs)
keyring = results.keys()
keyring.sort()
for item in keyring:
	print item, results[item]
