#!/usr/bin/python

import sys,os,math
import HAConf, GenerateCGSystem
from random import randint, shuffle

from xgrid_tools import *
starting_point = os.getcwd()

def write_batch_task(batch_file,mutant,model_type,seed=5,angle=0,bias=True,simlength=200,jobname="",stripe=False,random_lipids=False,position=0,system_size="XS",special=None,wallclock=24,lipid_types="DPPC",temperature=None):
			seed_name = "%09d" % (seed)
			position_name = "%09.3f" % (position)

			if not stripe:
				mutant_name = position_name + "_" + seed_name + "_" + model_type + "_" + lipid_types + "_" + mutant
			else:
				mutant_name = seed_name + "_" + mutant + "_" + position_name + "_" + model_type + "_" + lipid_types
			print >> batch_file, "\t"*3 + mutant_name + " = {"
			print >> batch_file, "\t"*4 + "arguments = ("
			print >> batch_file, "\t"*5 + '"-s",'
			print >> batch_file, "\t"*5 + '"' + mutant + '",'
			print >> batch_file, "\t"*5 + '"-t",'
			print >> batch_file, "\t"*5 + '"' + model_type + '",'
			print >> batch_file, "\t"*5 + '"--lipid_headgroup_mix",'
			print >> batch_file, "\t"*5 + '"' + lipid_types + '",'
			print >> batch_file, "\t"*5 + '"-r",'
			print >> batch_file, "\t"*5 + '"' + str(seed) + '",'
			print >> batch_file, "\t"*5 + '"-a",'
			print >> batch_file, "\t"*5 + '"' + str(angle) + '",'
			print >> batch_file, "\t"*5 + '"-p",'
			print >> batch_file, "\t"*5 + '"' + str(position) + '",'
			print >> batch_file, "\t"*5 + '"--wallclock",'
			print >> batch_file, "\t"*5 + '"' + str(wallclock) + '",'
			print >> batch_file, "\t"*5 + '"-l",'
			print >> batch_file, "\t"*5 + '"' + str(simlength) + '",'
			print >> batch_file, "\t"*5 + '"-j",'
			print >> batch_file, "\t"*5 + '"' + jobname + '",'
			if not bias:
				print >> batch_file, "\t"*5 + '"-u",'
			if random_lipids:
				print >> batch_file, "\t"*5 + '"--randomize_lipids",'
			print >> batch_file, "\t"*5 + '"--system_size",'
			print >> batch_file, "\t"*5 + '"' + system_size + '",'
			if special:
				print >> batch_file, "\t"*5 + '"--special",'
				print >> batch_file, "\t"*5 + '"' + special + '",'
			if temperature:
				print >> batch_file, "\t"*5 + '"--change_temperature",'
				print >> batch_file, "\t"*5 + '"' + str(temperature) + '",'
			print >> batch_file, "\t"*5 + '"-b"'
			print >> batch_file, "\t"*4 + ");"
			print >> batch_file, "\t"*4 + 'command = "' + HAConf.programs['CG_Helix_UmbrellaSA'][:-1] + '";'
			#print >> batch_file, "\t"*4 + 'command = "/Users/Shared/Sidekick/CG_Helix.py";'
			print >> batch_file, "\t"*3 + "};"

#Get options

import SidekickOptParse
#options = SidekickOptParse.get_options(("-j",))
parser = SidekickOptParse.Multi_UmbrellaSA_Helix_Options(("-j",))
options = parser.options

#Determine if this is a scan, or input from a file

if options.sequence_file == None and options.sequence == None:
	print "This program requires input in terms of a file of fasta sequences (-f) or a sequence (-s) to repeat (-d), or mutate (-m)"
	print "Try -h as an option for more info"
	exit()
elif options.sequence_file != None and options.sequence != None:
	print "This only works on a sequence, or a file of sequences"
	print "Try -h as an option for more info"
	exit()

if options.batch:
	print "This program always runs batch jobs by default. The batch option (-b) is unnecessary"

#Models or model scan?

if options.model_type not in GenerateCGSystem.installed_models:
	if options.model_type == "search":
		models = GenerateCGSystem.installed_models
	else:
		print "Choose an installed model from the list below or 'search' to run all models"
		for item in GenerateCGSystem.installed_models:
			print '\t' +  item
		print "Try -h as an option for more info"
		exit()
else:
	models = [options.model_type,]

if options.system_size not in GenerateCGSystem.defined_system_sizes.keys():
	print "Size not available at present. Please restart with one of the following options for \"--system-size\""
	for model in GenerateCGSystem.defined_system_sizes.keys():
		print "\t"+model, "~" + str(GenerateCGSystem.defined_system_sizes[model]["lipid"]*64), "lipids"
	print "Try -h as an option for more info"
	exit()
	
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

angle_ranges = [i%13*15 for i in range(number_of_repeats)]

#windows need to be expressed in nm explicitly

number_of_windows = int(math.ceil(options.window_max/options.window_size)) + 1
windows = [item * options.window_size for item in range(number_of_windows)]

shuffle(angle_ranges)
seed_angle = zip(seed_packet,angle_ranges)
seed_angle_combos = []
for item in windows:
	for jitem in seed_angle:
		seed_angle_combos.append(jitem+(item,))

if options.binary_lipid_scan:
	lipid_options = options.binary_lipid_scan.split("/")
	lipid_types = lipid_options[0]
	lipid_steps = [int(item) for item in lipid_options[1].split(":")]
	if lipid_steps[-1] > 100:
		print "Cannot get more than 100 %. Exiting"
		sys.exit()
	lipid_mixtures = [ lipid_types + "/" + str(item)+":"+str(100 - item) for item in range(lipid_steps[0],lipid_steps[2]+lipid_steps[1],lipid_steps[1]) ]

if options.sequence_file:
	sequences = [line[:-1] for line in open(options.sequence_file,'r')]
elif options.mutant:
	sequences = []
	for i in range(len(options.sequence)):
		scan_residues = options.mutant
		if i+len(scan_residues) > len(options.sequence): continue
		mutant = options.sequence[:i] + scan_residues + options.sequence[(i+len(scan_residues)):]
		if mutant not in sequences:
			sequences.append(mutant)
else:
	sequences = [options.sequence,]

#Rudimentary test for sequence sanity...

vetted_sequences= []

for item in sequences:
	if len(item) > 0 :
		vetted_sequences.append(item)

sequences = vetted_sequences

#Everything should be in place now- we have x sequences, being repeated with y seeds. Write a batch file to the current directory. If you try to create the job directory it gets confused about permissions
'''os.chdir(HAConf.results_location)
try:
	os.mkdir(options.destination)
except:
	pass
os.chdir(options.destination)'''

#Generate batch file

batch_file = open(options.destination + ".batch","w")

##Write Header

print >> batch_file, """{
	jobSpecification = {
		name = \"""" + options.destination + """\";
		taskSpecifications = {"""

#if options.stripe:
#	for model_type in models:
#		for seed in seed_packet:
#			for mutant in sequences:
				#write_batch_task(batch_file,mutant,model_type,seed=seed,bias=options.bias,simlength=options.simlength,jobname=options.destination,stripe=True,system_size=options.system_size,random_lipids=options.random_lipids)
#else:
for model_type in models:
		for mutant in sequences:
			#for seed in seed_packet:
			if options.binary_lipid_scan:
				for mix in lipid_mixtures:
					for combo in seed_angle_combos:
						seed = combo[0]
						if options.vary_angle:
							angle = combo[1]
						else:
							angle = 0
						
						position = combo[2]

						write_batch_task(batch_file,mutant,model_type,seed=seed,bias=options.bias,
							random_lipids=options.random_lipids,simlength=options.simlength,jobname=options.destination,
							angle=angle,position=position,system_size=options.system_size,stripe=options.stripe,
							special=options.special,wallclock=options.wallclock,lipid_types=mix,temperature=options.temperature)
				
			else:
				for combo in seed_angle_combos:
					seed = combo[0]
					if options.vary_angle:
						angle = combo[1]
					else:
						angle = 0

					position = combo[2]

					write_batch_task(batch_file,mutant,model_type,seed=seed,bias=options.bias,
						random_lipids=options.random_lipids,simlength=options.simlength,jobname=options.destination,
						angle=angle,position=position,system_size=options.system_size,stripe=options.stripe,
						special=options.special,wallclock=options.wallclock,lipid_types=options.headgroup_mix,temperature=options.temperature)

print >> batch_file, "\t"*2 + "};"
print >> batch_file, "\t" + "};"
print >> batch_file, "}"

batch_file.close()

#Submit batch job

job_list.append(batch_submit(options.destination+".batch"))