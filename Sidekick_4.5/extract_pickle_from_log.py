#!/usr/bin/python

import pickle

old_list = []
job_list = []

for line in open("daemon.log","r"):
  if line[0] == "#": old_list = job_list; job_list = []; continue
  contents = line.split()
  if len(contents) != 3: continue
  if contents[1] == "Finished": continue
  job_list.append((int(contents[0]), contents[2]))

if len(job_list) < len(old_list): job_list = old_list

restart_file = open("restart.pickle","w")
pickle.dump(job_list,restart_file)
restart_file.close()
