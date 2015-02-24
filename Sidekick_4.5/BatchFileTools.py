#!/usr/bin/python

import HAConf
'''
import platform
#Leopard has a different format for batch files
if platform.release().split(".")[0] == "11":
		def open_batch_file(batch_file,jobname):
			print >> batch_file, """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple Computer//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<array>
	<dict>
		<key>name</key>
			<string>""" + jobname + """</string>
             <key>taskSpecifications</key>
             <dict>"""
		
		def close_batch_file(batch_file):
			print >> batch_file, "\t"*2 + "</dict>"
			print >> batch_file, "\t" + "</dict>"
			print >> batch_file, "</array>"
			print >> batch_file, "</plist>"

		def write_dimer_batch_task(batch_file,mutant1,mutant2,model_type,seed=5,simlength=200,jobname="",stripe=False,special=None,wallclock=320,temperature=None):
					seed_name = "%09d" % (seed)
					if stripe:
						mutant_name = seed_name + "_" + model_type + "_" + mutant1 + "_" + mutant2
					else:
						mutant_name = mutant1 + "_" + mutant2 + "_" + model_type + "_" + seed_name
					print >> batch_file, "\t"*3 + "<key>" + mutant_name + "</key>"
					print >> batch_file, "\t"*3 + "<dict>"
					print >> batch_file, "\t"*4 + "<key>command</key>"
					print >> batch_file, "\t"*4 + "<string>" + HAConf.programs['CG_Helix_Dimer'][:-1] + "</string>"
					print >> batch_file, "\t"*4 + "<key>arguments</key>"
					print >> batch_file, "\t"*4 + "<array>"
					print >> batch_file, "\t"*5 + '<string>--s1</string>'
					print >> batch_file, "\t"*5 + '<string>' + mutant1 + '</string>'
					print >> batch_file, "\t"*5 + '<string>--s2</string>'
					print >> batch_file, "\t"*5 + '<string>' + mutant2 + '</string>'
					print >> batch_file, "\t"*5 + '<string>-t</string>'
					print >> batch_file, "\t"*5 + '<string>' + model_type + '</string>'
					print >> batch_file, "\t"*5 + '<string>-r</string>'
					print >> batch_file, "\t"*5 + '<string>' + str(seed) + '</string>'
					print >> batch_file, "\t"*5 + '<string>--wallclock</string>'
					print >> batch_file, "\t"*5 + '<string>' + str(wallclock) + '</string>'
					print >> batch_file, "\t"*5 + '<string>-l</string>'
					print >> batch_file, "\t"*5 + '<string>' + str(simlength) + '</string>'
					print >> batch_file, "\t"*5 + '<string>-j</string>'
					print >> batch_file, "\t"*5 + '<string>' + jobname + '</string>'
					if special:
						print >> batch_file, "\t"*5 + '<string>--special</string>'
						print >> batch_file, "\t"*5 + '<string>' + special + '</string>'
					if temperature:
						print >> batch_file, "\t"*5 + '<string>--change_temperature</string>'
						print >> batch_file, "\t"*5 + '<string>' + str(temperature) + '</string>'
					print >> batch_file, "\t"*5 + '<string>-b</string>'
					print >> batch_file, "\t"*4 + "</array>"
					#print >> batch_file, "\t"*4 + 'command = "/Users/Shared/Sidekick/CG_Helix.py";'
					print >> batch_file, "\t"*3 + "</dict>"

			
else:

'''
def open_batch_file(batch_file,jobname):
			print >> batch_file, """{
	jobSpecification = {
		name = \"""" + jobname + """\";
		taskSpecifications = {"""
		
def close_batch_file(batch_file):
			print >> batch_file, "\t"*2 + "};"
			print >> batch_file, "\t" + "};"
			print >> batch_file, "}"
			
def write_batch_task(batch_file,mutant,model_type,seed=5,angle=0,bias=True,simlength=200,jobname="",stripe=False,random_lipids=False,position=0,system_size="XS",special=None,wallclock=24,lipid_types="DPPC",temperature=None):
					seed_name = "%09d" % (seed)
					if stripe:
						mutant_name = seed_name + "_" + model_type + "_" + lipid_types + "_" + mutant
					else:
						mutant_name = mutant + "_" + model_type + "_" + seed_name + "_" + lipid_types
					print >> batch_file, "\t"*3 + mutant_name + " = {"
					print >> batch_file, "\t"*4 + "arguments = ("
					print >> batch_file, "\t"*5 + '"-s",'
					print >> batch_file, "\t"*5 + '"' + mutant + '",'
					print >> batch_file, "\t"*5 + '"-t",'
					print >> batch_file, "\t"*5 + '"' + model_type + '",'
					print >> batch_file, "\t"*5 + '"--lipid_headgroup_mix",'
					print >> batch_file, "\t"*5 + '"' + lipid_types + '",'
					print >> batch_file, "\t"*5 + '"-r",'
					print >> batch_file, "\t"*5 + '"' + str(seed) + '",'
					print >> batch_file, "\t"*5 + '"-a",'
					print >> batch_file, "\t"*5 + '"' + str(angle) + '",'
					print >> batch_file, "\t"*5 + '"-p",'
					print >> batch_file, "\t"*5 + '"' + str(position) + '",'
					print >> batch_file, "\t"*5 + '"--wallclock",'
					print >> batch_file, "\t"*5 + '"' + str(wallclock) + '",'
					print >> batch_file, "\t"*5 + '"-l",'
					print >> batch_file, "\t"*5 + '"' + str(simlength) + '",'
					print >> batch_file, "\t"*5 + '"-j",'
					print >> batch_file, "\t"*5 + '"' + jobname + '",'
					if not bias:
						print >> batch_file, "\t"*5 + '"-u",'
					if random_lipids:
						print >> batch_file, "\t"*5 + '"--randomize_lipids",'
					print >> batch_file, "\t"*5 + '"--system_size",'
					print >> batch_file, "\t"*5 + '"' + system_size + '",'
					if special:
						print >> batch_file, "\t"*5 + '"--special",'
						print >> batch_file, "\t"*5 + '"' + special + '",'
					if temperature:
						print >> batch_file, "\t"*5 + '"--change_temperature",'
						print >> batch_file, "\t"*5 + '"' + str(temperature) + '",'
					print >> batch_file, "\t"*5 + '"-b"'
					print >> batch_file, "\t"*4 + ");"
					print >> batch_file, "\t"*4 + 'command = "' + HAConf.programs['CG_Helix'][:-1] + '";'
					#print >> batch_file, "\t"*4 + 'command = "/Users/Shared/Sidekick/CG_Helix.py";'
					print >> batch_file, "\t"*3 + "};"
		
		
def write_dimer_batch_task(batch_file,mutant1,mutant2,model_type,seed=5,simlength=200,jobname="",stripe=False,special=None,wallclock=320,temperature=None,parallel=True,rotatebyangle=0):
					seed_name = "%09d" % (seed)
					if stripe:
						mutant_name = seed_name + "_" + model_type + "_" + mutant1 + (40-len(mutant1))*"x" + "_" + mutant2 + (40-len(mutant2))*"x"
					else:
						mutant_name = mutant1 + (40-len(mutant1))*"x" + "_" + mutant2 + (40-len(mutant2))*"x"+ "_" + model_type + "_" + seed_name
					print >> batch_file, "\t"*3 + mutant_name + " = {"
					print >> batch_file, "\t"*4 + "arguments = ("
					print >> batch_file, "\t"*5 + '"--s1",'
					print >> batch_file, "\t"*5 + '"' + mutant1 + '",'
					print >> batch_file, "\t"*5 + '"--s2",'
					print >> batch_file, "\t"*5 + '"' + mutant2 + '",'
					print >> batch_file, "\t"*5 + '"-t",'
					print >> batch_file, "\t"*5 + '"' + model_type + '",'
					print >> batch_file, "\t"*5 + '"-r",'
					print >> batch_file, "\t"*5 + '"' + str(seed) + '",'
					print >> batch_file, "\t"*5 + '"--wallclock",'
					print >> batch_file, "\t"*5 + '"' + str(wallclock) + '",'
					print >> batch_file, "\t"*5 + '"-l",'
					print >> batch_file, "\t"*5 + '"' + str(simlength) + '",'
					print >> batch_file, "\t"*5 + '"-j",'
					print >> batch_file, "\t"*5 + '"' + jobname + '",'
					if special:
						print >> batch_file, "\t"*5 + '"--special",'
						print >> batch_file, "\t"*5 + '"' + special + '",'
					if temperature:
						print >> batch_file, "\t"*5 + '"--change_temperature",'
						print >> batch_file, "\t"*5 + '"' + str(temperature) + '",'
					if not parallel:
						print >> batch_file, "\t"*5 + '"-a",'
					print >> batch_file, "\t"*5 + '"--rotatebyangle",'
					print >> batch_file, "\t"*5 + '"' + str(rotatebyangle) + '",'
					print >> batch_file, "\t"*5 + '"-b"'
					print >> batch_file, "\t"*4 + ");"
					print >> batch_file, "\t"*4 + 'command = "' + HAConf.programs['CG_Helix_Dimer'][:-1] + '";'
					#print >> batch_file, "\t"*4 + 'command = "/Users/Shared/Sidekick/CG_Helix.py";'
					print >> batch_file, "\t"*3 + "};"
