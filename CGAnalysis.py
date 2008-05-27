#!/usr/bin/python

import os
import HAConf
import AutomatedPlot

def run(command_string, keystrokes=''): 
	"""run a UNIX command but use os.popen, rather than os.system so that key strokes can be passed to the process. e.g. if you run make_ndx in gromacs you will be asked for the group number. With run you can include it in the call. Dont forget to put '\\n' after your keystrokes if they are needed"""
	command=os.popen(command_string,'w')
	command.flush() 
	command.write(keystrokes) 
	command.flush() 
	command.close()

gromacs = HAConf.configuration["gromacs_location"]

length = int(open("ssdump.ssd").readline())

make_ndx_command = 'echo "r 1-7\n r ' + str(length-6) + '-' + str(length) +'\n q\n" | '  + gromacs + 'make_ndx -f t_0 -o system.ndx'

g_bundle_command = "echo '16\n17\n'| g_bundle -f t_0 -s t_0 -na 1 -z -n system.ndx "


os.system(make_ndx_command)
#print make_ndx_command
os.system(g_bundle_command)

AutomatedPlot.tilt_data_plot("bun_tilt.xvg","tilt_V_time.png","sequence")
AutomatedPlot.hist_tilt_data_plot("bun_tilt.xvg","hist_tilt.png","sequence")