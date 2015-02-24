#!/usr/bin/python

import sys,os
import HAConf

from xgrid_tools import *

starting_point = os.getcwd()
input_file = sys.argv[1]

job_list =[]

for line in open(input_file,'r'):
        mutant = line[:-1]
        print mutant
        try:
                os.mkdir(mutant)
        except:
                continue
        os.chdir(mutant)
        
        #Now create run directories and run
        os.mkdir("Rotation-Translation")
        os.chdir("Rotation-Translation")
        job_list.append(job_submit(HAConf.programs['hippo_tr'] + " " + mutant + " -dAgg "))
        os.chdir("..")
        
        os.mkdir("MARTINI")
        os.chdir("MARTINI")
        job_list.append(job_submit(HAConf.programs['MARTINI'] + " " + mutant + " " + str(5) + " -dAgg "))
        os.chdir("..")
        
        #Finally, return to the starting directory
        os.chdir(starting_point)
        
collect_results_daemon(job_list,starting_point+"/daemon.log",starting_point+"/restart.pickle")

