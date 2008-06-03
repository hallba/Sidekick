#!/usr/bin/python

import os

debug_status = True

configuration = {	'gromacs_location': '/usr/local/gromacs-3.3.3/bin/',
					'mdp_files':		'/Users/benjamin/Work/Runfiles/',
					'vmd_state_files':	'/Users/benjamin/Work/Sidekick/VMD/',
					'hippo_files':		'/Users/benjamin/Work/Hippo_runfiles/',
					'martini_systems':	'/Users/benjamin/Work/MARTINI/systems/',
					'martini_setup':	'/Users/benjamin/Work/MARTINI/Setup/',
					'lipids':			'/Users/benjamin/Work/MARTINI/Setup/AddLipid+Solvent/'
			}

programs = {		'pymol':			'/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL ',
					'MARTINI':			'/Users/benjamin/Work/Tar/CGAuto.py ',
					'hippo_tr':			'/Users/benjamin/Work/Tar/HippoAuto.py ',
					'vmd':				'/Applications/VMD\ 1.8.6.app/Contents/MacOS/startup.command '

			}

queue_initiator = ' '


for item in configuration.keys():
	if not os.path.isdir(configuration[item]): print configuration[item], "missing. Exiting"; exit() 

if __name__=="__main__":
	print "This is a configuration file. All values are stored in a dictionary and are printed below"
	for item in configuration.keys():
		print item+"\t"+configuration[item]
	
