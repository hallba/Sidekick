#!/usr/bin/python

from __future__ import with_statement
import os, base64, pickle, plistlib, sys

debug_status = False
thread = 2
try:
	configuration = plistlib.readPlist(sys.path[0]+"/conf.plist")
	sidekick_location = configuration['sidekick_location']
	results_location = configuration['results_location']
	#results_location = "/Volumes/GladOS Data"
	web_data_location = configuration['web_data_location']
	visualisation_location = configuration['visualisation_location']
	public_location = configuration['public_location']
	python_location = configuration['python']

	with open(sidekick_location + "/pw.pickle") as inp:
		password = base64.b64decode(pickle.load(inp))
	
	xgrid = { 			'password': 		password,
						'controller':		configuration['controller']
				}
except:
	#fallback defaults
	sidekick_location = "/Users/Shared/Sidekick"
	results_location = "/Users/Shared/Data"
	#results_location = "/Volumes/GladOS Data"
	web_data_location = "/tmp/Web_Sidekick"
	visualisation_location = "/Library/WebServer/Documents/Results"
	public_location = "http://glados.bioch.ox.ac.uk/Results/"
	python_location = "/usr/bin/python"

	with open(sidekick_location + "/pw.pickle") as inp:
		password = base64.b64decode(pickle.load(inp))
	
	xgrid = { 			'password': 		password,
						'controller':		'localhost'
				}

#installed_models = ["MARTINI", "MARTINI_1.1.1", "MARTINI_1.1.2", "MARTINI_1.1.2.b", "MARTINI_1.1.2.b.ENM", "Bond", "Bond0.9.5","BondDev","BondUnchargedDE"]

#from GenerateCGSystem import installed_models

configuration = {	'gromacs_location': sidekick_location + '/gromacs/bin/',
					'mdp_files':		sidekick_location +'/MDP/',
					#'vmd_state_files':	sidekick_location +'/VMD/',
					#'hippo_files':		sidekick_location +'/hippo_runfiles/',
					'martini_systems':	sidekick_location +'/_MARTINI/',
					'martini_setup':	sidekick_location +'/MARTINI/',
					'DPlipids':			sidekick_location +'/MARTINI/AddLipid+Solvent/',
					'DLlipids':			sidekick_location +'/MARTINI/DL_Lipid+Solvent/',
					'POlipids':			sidekick_location +'/MARTINI/PO_Lipid+Solvent/',
					'DOlipids':			sidekick_location +'/MARTINI/DO_Lipid+Solvent/',
					'misc_scripts':		sidekick_location +'/Scripts/',
					'bond_setup':		sidekick_location +'/BOND/',
					'pmf':				sidekick_location +'/PMF/',
					'python':		python_location
			}

programs = {		#'pymol':			sidekick_location +'/pymol/pymol ',
					'MARTINI':			sidekick_location +'/Bilayer_Insertion.py  ',
					#'hippo_tr':			sidekick_location +'/HippoAuto.py ',
					#'vmd':				sidekick_location +'/VMD_install/vmd ',
					'Insertion_Efficiency':		sidekick_location +'/Insertion_Efficiency.py ',
					'CG_Helix':			sidekick_location + '/CG_Helix.py ',
					'CG_Helix_Dimer':			sidekick_location + '/CG_Helix_Dimer.py ',
					'PMF_CG_Helix':			sidekick_location + '/CG_Helix_PMF.py ',
					'CG_Helix_UmbrellaSA':			sidekick_location + '/CG_Helix_UmbrellaSA.py '
			}

if len(xgrid['password']) > 0:
	queue_initiator = 'xgrid -job run -h ' + xgrid['controller'] + " -p " + xgrid['password'] + " -so xgrid_log.txt -se xgrid_err.txt -out . " 
else:
	queue_initiator = 'xgrid -job run -h ' + xgrid['controller'] + " -so xgrid_log.txt -se xgrid_err.txt -out ." 

if False:
	for item in configuration.keys():
		if not os.path.isdir(configuration[item]): print configuration[item], "missing. Exiting"; exit() 

if __name__=="__main__":
	print "This is a configuration file. All values are stored in a dictionary and are printed below"
	for item in configuration.keys():
		print item+"\t"+configuration[item]
