#!/usr/bin/python

import os, HAConf, shutil

working_dir = HAConf.web_data_location #"/tmp/sidekick-web"
archive_dir = HAConf.web_data_location + "/archive"

def job_submit(batch_file):
	base_command = "xgrid -h "+HAConf.xgrid["controller"]
	if "password" in HAConf.xgrid:
		base_command += " -p "+ HAConf.xgrid["password"]
	base_command += " -job batch "+ working_dir + "/" + batch_file
	os.system(base_command)
	shutil.move(working_dir + "/" + batch_file,archive_dir + "/" + batch_file)

if not os.path.exists(working_dir):
	os.mkdir(working_dir)
os.system("chmod a+rwx " + working_dir)
	
if not os.path.exists(archive_dir):
	os.mkdir(archive_dir)
os.system("chmod a+rwx " + archive_dir)

new_batch = os.listdir(working_dir)

for item in new_batch:
	if item != "archive":
		job_submit(item)