#!/usr/bin/python

'''ReQueue'''

import sys,os,shutil,plistlib
import HAConf, GenerateDimerCGSystem
from random import randint, shuffle

from xgrid_tools import *

broken_jobs = sys.argv[1:]

start = []
taskspec = {}

for job in broken_jobs:
	task_name = "rescue_" + job
	taskspec[task_name] = { 	"command"		:		HAConf.sidekick_location+"/RestartManager.py",
								"arguments"		:		[	job,	]
	}
	
start.append(	{	"name"					:	"rescue",
					"taskSpecifications"	:	taskspec
					})
	
filename = os.getcwd() + "/rescue.batch"
plistlib.writePlist(start, filename)
os.system("chmod a+w "+filename)
