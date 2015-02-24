#!/usr/bin/python

import sys,os
import HAConf, GenerateDimerCGSystem
from random import randint, shuffle

from xgrid_tools import *
starting_point = os.getcwd()

from BatchFileTools import open_batch_file, close_batch_file
from BatchFileTools import write_dimer_batch_task as write_batch_task

#Get options

import SidekickOptParse
#options = SidekickOptParse.get_options(("-j",))
parser = SidekickOptParse.Multi_Helix_Dimer_Options(("-j",))
options = parser.options

if options.sequence_file == None and (options.sequence1 == None or options.sequence2 == None):
	print "This program requires input in terms of a file of fasta sequences (-f) or sequences (--s1,--s2) to repeat (-d), or mutate (-m)"
	print "Try -h as an option for more info"
	exit()
elif options.sequence_file != None and options.sequence1 != None and options.sequence2 != None:
	print "This only works on a sequence, or a file of sequences"
	print "Try -h as an option for more info"
	exit()


if options.batch:
	print "This program always runs batch jobs by default. The batch option (-b) is unnecessary"




#Models or model scan?

if options.model_type not in GenerateDimerCGSystem.installed_models:
	if options.model_type == "search":
		models = GenerateDimerCGSystem.installed_models.keys()
	else:
		print "Choose an installed model from the list below or 'search' to run all models"
		for item in GenerateDimer.installed_models:
			print '\t' +  item
		print "Try -h as an option for more info"
		exit()
else:
	models = [options.model_type,]

#If the input is sensible, generate a list of seeds, a list of sequences, and generate a batch file in the target directory

number_of_repeats = options.duplicates

job_list =[]
#seed_packet = [randint(0,99999) for i in range(number_of_repeats)]
seed_packet = []
old_seeds = []#Germinated! Ho ho ho

try:
	#Find out what seeds have been already run for jobs with this job name
	job_location = HAConf.results_location + "/" + options.destination
	
	for i in range(4):
		job_location += "/" + os.listdir(job_location)[0]
	old_seeds = [int(i.split("-")[0]) for i in os.listdir(job_location) ]
except:
	print "No old seeds to analyse"

for i in range(number_of_repeats):
	tmp_rand = randint(0,99999)
	while tmp_rand in seed_packet or tmp_rand in old_seeds:
		tmp_rand = randint(0,99999)
	seed_packet.append(tmp_rand)

if options.vary_angle:
	seed_combo = [[ seed_packet[i] , 15 *(i %24) ] for i in range(number_of_repeats)]
else:
	seed_combo = [[seed,0] for seed in seed_packet  ]

if options.sequence_file:
	sequences = [line.split() for line in open(options.sequence_file,'r')]
else:
	sequences = [[options.sequence1,options.sequence2],]
	
#Generate batch file

batch_file = open(options.destination + ".batch","w")

##Write Header

open_batch_file(batch_file,options.destination)

for model_type in models:
	for pair in sequences:
		for [seed,angle] in seed_combo:
			write_batch_task(	batch_file,
								pair[0],
								pair[1],
								model_type,
								seed=seed,
								simlength=options.simlength,
								jobname=options.destination,
								stripe=options.stripe,
								special=options.special,
								wallclock=options.wallclock,
								temperature=options.temperature,
								parallel=options.parallel,
								rotatebyangle=angle)
								
close_batch_file(batch_file)

batch_file.close()

#Submit batch job

job_list.append(batch_submit(options.destination+".batch"))
