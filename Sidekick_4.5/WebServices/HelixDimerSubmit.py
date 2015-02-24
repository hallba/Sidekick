#!/usr/bin/python

#HelixDimerSubmit.py

from cgi import *
import sys, re, os, plistlib
from random import randint

#set the stage

configuration = plistlib.readPlist("info.plist")

try:
	working_dir = configuration["web_data_location"]
except:
	working_dir = "/tmp/sidekick-web"

if not os.path.exists(working_dir):
	os.mkdir(working_dir)

os.system("chmod a+rwx " + working_dir)
	
sys.stderr=sys.stdout


def input_clean(sequences):
	clean_sequences = []
	for item in sequences:
		if len(item) <= 2: continue
		if item[-1] == "\r": item = item[:-1]
		if item[-1] == "\n": item = item[:-1]
		clean_sequences.append(item)
	return clean_sequences

def grab_data():

    global form

    details = {	"job_name"		:	"",
    			"sequenceInput"	:	"",
    			"length"		:	"",
    			"stripe"		:	"",
    			"forcefield"	:	"",
    			"duplicates"	:	"",
    			"special"		:	"",
    			"orientation"	:	""
    			}
    
    form =  FieldStorage()
    print "<b>Input values:</b>"
    print "<table>"
    print "<th>Input</th><th>Value</th>"
    for item in details.keys():
    	details[item] = form.getvalue(item)
    	print "<tr>","<td>",item,"</td>","<td>", details[item],"</td>","</tr>"
    print "</table>"
    return details

def generate_plist(input_options):
	start = []
	seed_packet = []
	#print "duplicates", input_options["duplicates"], int(input_options["duplicates"])
	for i in range(int(input_options["duplicates"])):
		tmp_rand = randint(0,99999)
		while tmp_rand in seed_packet:
			tmp_rand = randint(0,99999)
		seed_packet.append(tmp_rand)

	sequences = [line for line in input_options["sequenceInput"].split("\n")]
	sequences = input_clean(sequences)
	
	def generate_job_name(sequence,seed,forcefield):
		seed_name = "%09d" % (seed)
		if input_options["stripe"] == "False":
			[seed_name,sequence] = [sequence,seed_name]
			
		name = "%s_%s_%s" % (seed_name, forcefield, sequence)
		
		return name
	
	taskspec = {}
	
	for sequence in sequences:
		for seed in seed_packet:
			task_name = generate_job_name(sequence,seed,input_options["forcefield"])
			sequence_pair = sequence.split()
			taskspec[task_name] = { 	"command"	:	configuration['sidekick_location']+"/CG_Helix_Dimer.py",
										"arguments"	:	[	"--s1",	sequence_pair[0],
															"--s2",	sequence_pair[1],
															"-t",	input_options["forcefield"],
															"-r",	str(seed),
															"-l",	input_options["length"],
															"-j",	input_options["job_name"],
															]
										}
			if len(input_options["special"]) > 0:
				taskspec[task_name]["arguments"] += ["--special",input_options["special"]]
			if input_options["orientation"] == "A":
				taskspec[task_name]["arguments"].append("-a")
			taskspec[task_name]["arguments"].append("-b")
	
	start.append(	{	"name"					:	input_options["job_name"],
					"taskSpecifications"	:	taskspec
					})
			
			
	
	filename = working_dir + "/"+input_options["job_name"]+".batch"
	plistlib.writePlist(start, filename)
	os.system("chmod a+w "+filename)




if __name__ == "__main__":
	print 'Content-type: text/html\n\n'
	print '<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">'
	print '\n'
	print """<html>"""
	print """<head>"""
	print """<title>Sidekick for Single Helices</title>"""
	print """<meta name="viewport" content="width=device-width; initial-scale=1.0; maximum-scale=1.0;">"""
	print """<link rel="apple-touch-icon-precomposed" href="images/sidekick_dimer.png"/>"""
	print """</head>"""
	print """<h1>Checking your submitted data... </h1>"""


	input_options = grab_data()
	generate_plist(input_options)
	#print_confirmation(input_options)
	print "<p>Your job results will appear at</p>"
	resloc="http://glados.bioch.ox.ac.uk/Results/"+ input_options["job_name"] +"/"
	print '<p><b><a href="' + resloc + '">' + resloc + "</a></b></p>"
	#print resloc
	print "<p>Get status updates at</p>"
	resloc += "status.html"
	print '<p><b><a href="' + resloc + '">' + resloc + "</a></b></p>"
	print """</html>"""

