#!/usr/bin/python

import os, sys, pickle
import xgrid_tools

starting_point = os.getcwd()

#Determine pid of the old job and kill

log_file = open("daemon.log", "r")
start_line = log_file.readline()

pid = int(start_line[24:])

##Detect restart pid
for line in log_file:
	if line[:6] == "Daemon": pid = int(start_line[24:])
log_file.close()

print "Old monitor pid", pid

##Try to kill old daemon (I need to check that the pid is reported correctly)
try:
	os.kill(pid,15)
except:
	pass

#Now load restart pickle and spin off a fresh daemon, again writing to daemon.log/restart.pickle
job_list = pickle.load(open("restart.pickle",'r'))
print "Pickle restored\nJob restarting"
xgrid_tools.collect_results_daemon(job_list,starting_point+"/daemon.log",starting_point+"/restart.pickle")
