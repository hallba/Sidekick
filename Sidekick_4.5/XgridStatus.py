#!/usr/bin/python

from __future__ import with_statement
import HAConf
import subprocess,os,time

def job_attributes(jid):
	time.sleep(10)
	base_command = "xgrid -h "+HAConf.xgrid["controller"]
	if "password" in HAConf.xgrid:
		base_command += " -p "+ HAConf.xgrid["password"]
	base_command += " -job attributes -id " + str(jid)
	#print base_command
	result = subprocess.Popen(base_command, shell=True, stdout=subprocess.PIPE)
	job_description = {}
	for line in result.stdout:
		if "{" in line or "}" in line or "(" in line or ")" in line:
			continue
		contents = line.split()
		job_attribute = contents[0]
		if '"' in line:
			attribute_contents = line.split('"')[1]
		else:
			attribute_contents = contents[-1][:-1]
		job_description[job_attribute] = attribute_contents
	return job_description

def grid_power():
	time.sleep(10)
	base_command = "xgrid -h "+HAConf.xgrid["controller"]
	if "password" in HAConf.xgrid:
		base_command += " -p "+ HAConf.xgrid["password"]
	base_command += " -grid attributes -gid 0"
	result = subprocess.Popen(base_command, shell=True, stdout=subprocess.PIPE)
	for line in result.stdout:
		if "gridMegahertz" in line:
			#print line
			#print line.split()
			#print line.split()[-1]
			#print line.split()[-1][:-1]
			return float(line.split()[-1][:-1])
	return 0
'''	
xgrid -h localhost -p thecakeisalie -grid attributes -gid 0
{
    gridAttributes =     {
        gridMegahertz = 496755;
        isDefault = YES;
        name = Xgrid;
    };
}
'''

def available_data():
	available=[]
	for item in os.listdir(HAConf.visualisation_location):
		#print item
		if os.path.exists(HAConf.visualisation_location + "/" + item+"/index.html"):
			available.append(item)
	return available

def get_disk_space():
	p = subprocess.Popen(['df','-h',HAConf.results_location],stdout=subprocess.PIPE)
	line= p.stdout.next()
	space= int(round(float(p.stdout.next().split()[4][:-1])))
	return space

def comma_clip(stuff):
	if stuff[-1] == ",":
		stuff = stuff[:-1]
	return stuff

def job_list():
	base_command = "xgrid -h "+HAConf.xgrid["controller"]
	if "password" in HAConf.xgrid:
		base_command += " -p "+ HAConf.xgrid["password"]
	base_command += " -job list"
	#print base_command
	result = subprocess.Popen(base_command, shell=True, stdout=subprocess.PIPE)
	jobs =[]
	for line in result.stdout:
		if "{" in line or "}" in line or "(" in line or ")" in line:
			continue
		if line[-1] == "\n":
			line = line[:-1]
		if line[-1] == ",":
			line = line[:-1]
		jobs.append(int(line))
	return jobs
	
def update_job_status():
	available_power = grid_power()
	pending = []
	running = []
	for job in job_list():
		attributes = job_attributes(job)
		if "/" in attributes["name"]:
			continue
		try:
			#print HAConf.visualisation_location + "/" + attributes["name"]
			os.mkdir(HAConf.visualisation_location + "/" + attributes["name"])
		except:
			pass
		status_file = HAConf.visualisation_location + "/" + attributes["name"] + "/status.html"
		with open(status_file,"w") as status_output:
			print >> status_output, """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>
<head>
<title>Sidekick @ University of Oxford</title>
<meta name="viewport" content="width=device-width; initial-scale=1.0; maximum-scale=1.0;">
<link rel="apple-touch-icon-precomposed" href="../../images/sidekick_ruby.png"/>
<link rel="stylesheet" type="text/css"
              href="../../Sidekick.css"
              title="standard" >
<script type="text/javascript">
<!--
function timedRefresh(timeoutPeriod) {
	setTimeout("location.reload(true);",timeoutPeriod);
}
//   -->
</script>
</head>

<div align="center">
<h1>""" + attributes["name"] + """</h1>
<p>
<table class="datasheet" border="3" cellpadding="10" cellspacing="0" width="100%" id="AutoNumber1">
  <tr class="even">

"""
			summary_location = HAConf.public_location + attributes["name"] + '/'
			status_location = summary_location + "status.html"
			gallery_location = summary_location + "gallery.html"
			grid_status_location = summary_location + "grid_status.html"
			#print >> status_output, "<th>Attribute</th><th>Status</th>"
			print >> status_output, "<td>Name</td><td>" + attributes["name"] + "</td>"
			print >> status_output, "</tr><tr class=\"odd\">"
			print >> status_output, "<td>Status</td><td>" + attributes["jobStatus"] + "</td>"
			print >> status_output, "</tr><tr class=\"even\">"
			print >> status_output, "<td>Percent Complete</td><td>" + attributes["percentDone"] + " %</td>"
			print >> status_output, "</tr><tr class=\"odd\">"
			if available_power == 0:
				print >> status_output, "<td>Active CPU Power</td><td>0 % (Grid off)</td>"
			else:
				print >> status_output, "<td>Active CPU Power</td><td>" + str(round((float(attributes["activeCPUPower"])/available_power*100),1)) + " %</td>"
			print >> status_output, "</tr></table></p><div align='center'>" + '<p><a href="' + status_location + '">Status</a> <a href="' + summary_location + '">Summary</a> <a href="' + gallery_location + '">Gallery</a> <a href="' + grid_status_location + '">Grid Status</a></p>' + "</body></html>"
			if attributes["jobStatus"] == "Pending" and len(running) > 0: pending.append(attributes["name"])
			if attributes["jobStatus"] == "Running": running.append( (attributes["name"],str(round((float(attributes["activeCPUPower"])/available_power*100))),str(round(float(attributes["percentDone"]),1))) )
	available = available_data()
	diskspace = get_disk_space()
	with open(HAConf.visualisation_location + "/status.xml","w") as grid_status_file:
		print >> grid_status_file, "<node>"
		print >> grid_status_file, "<diskspace>"+str(diskspace)+"</diskspace>"
		#for item in running:
		#	print >> grid_status_file, "<tr class=\"" + tableclass + "\"><td>" + item[0] + "</td><td>" + item[1] + " %</td><td>" + item[2] + " %</td></tr>"
		running_text = "<working>"
		usage_text = "<usage>"
		completion_text = "<completion>"
		for item in running:
			running_text += item[0] + ","
			completion_text += item[2] + ","
			usage_text += item[1] + ","
		completion_text = comma_clip(completion_text) + "</completion>"
		usage_text = comma_clip(usage_text) +"</usage>"
		running_text = comma_clip(running_text)+ "</working>"
		print >> grid_status_file, running_text
		print >> grid_status_file, completion_text
		print >> grid_status_file, usage_text
		pendingtext = "<pending>"
		for item in pending:
			pendingtext += item + ","
		pendingtext = comma_clip(pendingtext) + "</pending>"
		print >> grid_status_file, pendingtext
		availtext = "<available>"
		for item in available:
			availtext += item + ","
		availtext = comma_clip(availtext) + "</available>"
		print >> grid_status_file, availtext
		print >> grid_status_file, "</node>"
	#print pending, running
if __name__ == "__main__":
	update_job_status()

'''
xgrid -h glados.bioch.ox.ac.uk -p thecakeisalie -job list
xgrid -h glados.bioch.ox.ac.uk -p thecakeisalie -job attributes -id 216
'''
