#!/usr/bin/python

import sys,os
import HAConf

from xgrid_tools import *

import SidekickOptParse
options = SidekickOptParse.get_options(("-f",))

starting_point = os.getcwd()
input_file = options.sequence_file

job_list =[]

for line in open(input_file,'r'):
        mutant = line[:-1]
        print mutant
        #job_list.append(job_submit(HAConf.programs['CG_Helix'] + " -t MARTINI -s " + mutant + " -b "))
        #job_list.append(job_submit(HAConf.programs['CG_Helix'] + " -t MARTINI2.1.1 -s " + mutant + " -b "))
        job_list.append(job_submit(HAConf.programs['CG_Helix'] + " -t Bond -s " + mutant + " -b "))
        job_list.append(job_submit(HAConf.programs['CG_Helix'] + " -t Bond0.9.5 -s " + mutant + " -b "))

        
collect_results_daemon(job_list,starting_point+"/daemon.log",starting_point+"/restart.pickle")

