#!/usr/bin/python

import sys,os
import HAConf
from random import randint
import time

starting_point = os.getcwd()

if "-d" in sys.argv:
	print "Daemonizing..."
	time.sleep(5)
	#To allow for very long gaps between submissions turn into a daemon to start with
	import daemonize
	daemonize.daemonize('/dev/null',starting_point + "/submission.log",starting_point + "/submission.log")
	print "Submission Daemon started with pid %d" % os.getpid()
        print "Started %s" % time.ctime()
	os.chdir(starting_point)

from xgrid_tools import *

input_file = sys.argv[1]
number_of_repeats = 50

job_list =[]

seed_packet = [randint(0,100000) for i in range(number_of_repeats)]

#Tree structure is different here
#Type/Sequence/Repeats

os.mkdir("Rotation-Translation")
os.chdir("Rotation-Translation")

sequences = [line[:-1] for line in open(starting_point + "/" + input_file,'r')]

for mutant in sequences:
        print mutant
	#mutant_directory = starting_point + "/" + mutant
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

for mutant in sequences:
	sequence_location = martini_location + "/" + mutant
	efficiency = 0
	for random_seed in seed_packet:
		efficiency_file = sequence_location + "/" + str(seed) + "/efficiency.dat"
		efficiency += float(efficiency_file.readline())
		efficiency_file.close()
	efficiency /= len(seed_packet)
	print >> summary_table, mutant, efficiency

summary_table.close()
