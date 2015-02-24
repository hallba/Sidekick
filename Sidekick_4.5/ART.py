#!/usr/bin/python

import os,sys,re

'''
This ART script is designed to try and set up the client for use by the Sidekick grid
It first attempts to mount the nfs exported directory, before reporting a score of 0, 1, or 2:
	0	No mounted nfs drive or no available cpp
	1	Workstation (no "glados" in hostname)
	2	Server
'''

def report_and_exit(score):
	print score
	sys.exit()

score = 0

location = "/nfsmount"

#mount_command = "mkdir " + location +" ; mount_nfs glados.oerc.ox.ac.uk:/Shared\ Items/Public/ " + location + " "
#os.system(mount_command)

if (os.path.exists(location + "/Sidekick") + os.path.exists("/usr/lib/cpp")) == 2:
	score += 1
else:
	report_and_exit(score)
	
hostname = os.uname()[1]

# and hostname[-9:] == ".ox.ac.uk"

if hostname[:6] == "glados":
	score += 1
	
report_and_exit(score)