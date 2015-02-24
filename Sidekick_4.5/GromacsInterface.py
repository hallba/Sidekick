#!/usr/bin/python
from __future__ import with_statement
import sys,os,pdbio

from HAConf import configuration,programs,debug_status
from LipidBox import add_solvent_dppc
from MDPGenerator import replace_seed, replace_steps, replace_pcoupling, replace_epsilon, replace_temperature, pmartini_update, pme_update
import AutomatedPlot
import process_manager
import HAConf
gromacs = configuration['gromacs_location']
scripts = configuration['misc_scripts']
#pymol = programs['pymol']
#vmd = programs['vmd']



def UserDefSimulation(topology_file,input_structure,output_name,mdpfile,wallclock=False,wallclock_lifetime=24,wallclock_poll=0.1,wallclock_time_unit="hours",umbrella=None,maxwarn=0):
	grompp_modifier = "" # for dumping extra options
	mdrun_modifier = ""
	if umbrella:
		#grompp_modifier += " -n %s.ndx" % umbrella # Don't need this + it can break the grompp command if not all needed gromacs files are there
		mdrun_modifier += " -pi %(name)s.ppa -po %(name)s.gmx.ppa -pd %(name)s.pdo -pn %(name)s.ndx" % {"name" : umbrella}
	#os.system(gromacs+"grompp -f " + mdpfile + " -c " + input_structure + " -p " + topology_file + " -o " + output_name + " >& grompp_" + output_name+ ".log" + grompp_modifier)
	grompp(options={"f":mdpfile,"c":input_structure,"p":topology_file,"o":output_name,"maxwarn":maxwarn},logfile="grompp_" + output_name+ ".log",mod=grompp_modifier)
	#os.system(gromacs+"mdrun -deffnm " + output_name)
	if wallclock:
		if debug_status:
			process_manager.wallclock(gromacs+"mdrun " + "-nt " + str(HAConf.thread) + " -deffnm " + output_name + mdrun_modifier,poll=0.5,lifetime=5,time_unit="minutes")
		else:
			process_manager.wallclock(gromacs+"mdrun " + "-nt " + str(HAConf.thread) + " -deffnm " + output_name + mdrun_modifier,poll=wallclock_poll,lifetime=wallclock_lifetime,time_unit=wallclock_time_unit)
	else:
		os.system(gromacs+"mdrun " + "-nt " + str(HAConf.thread) + " -deffnm " + output_name + mdrun_modifier)

def SDMinimise(topology_file="ioned_topol",input_structure="ioned.pdb",output_name="em",maxwarn=0):
	UserDefSimulation(topology_file,input_structure,output_name,configuration['mdp_files']+"cg-em_extended",maxwarn=maxwarn)

def MolecularDynamics(topology_file="ioned_topol",input_structure="em",output_name="t_0",seed=1993,steps=100,pcoupling="semiisotropic"):
	if debug_status:
		replace_steps(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","runfile.mdp")
	elif pcoupling == "semiisotropic":
		replace_seed(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","seeded.mdp",seed)
		replace_steps("seeded.mdp","runfile.mdp",steps)
	elif pcoupling == "anistropic":
		replace_seed(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","seeded.mdp",seed)
		replace_pcoupling("seeded.mdp","aniso_pcouple.mdp",coupling="anisotropic")
		replace_steps("aniso_pcouple.mdp","runfile.mdp",steps)
	#os.system(gromacs+"grompp -f runfile -c " + input_structure + " -p " + topology_file + " -o " + output_name + " >& grompp_" + output_name + ".log")
	#os.system(gromacs+"mdrun -deffnm " + output_name)
	UserDefSimulation(topology_file,input_structure,output_name,"runfile.mdp",wallclock=True)

def PMFMolecularDynamics(topology_file="ioned_topol",input_structure="pr",output_name="t_0",seed=1993,steps=2500000):
	if debug_status:
		replace_steps(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","runfile.mdp")
	else:
		replace_seed(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","seeded.mdp",seed)
		replace_steps("seeded.mdp","runfile.mdp",steps)
	os.system(gromacs+"grompp -f runfile -n pull_system.ndx -c " + input_structure + " -p " + topology_file + " -o " + output_name + " >& grompp_" + output_name + ".log")
	os.system(gromacs+"mdrun -pi umbrella.ppa -po umbrella.gmx.ppa -pd pull.pdo -pn pull_system.ndx -deffnm " + output_name)
	#UserDefSimulation(topology_file,input_structure,output_name,"runfile.mdp")

def PRMolecularDynamics(topology_file="ioned_topol",input_structure="em",output_name="t_0",seed=1993,steps=100):
	if debug_status:
		replace_steps(configuration['mdp_files'] + "cg-pr.mdp","runfile.mdp")
	else:
		replace_seed(configuration['mdp_files'] + "cg-pr.mdp","runfile.mdp",seed)
	#os.system(gromacs+"grompp -f runfile -c " + input_structure + " -p " + topology_file + " -o " + output_name + " >& grompp_" + output_name + ".log")
	#os.system(gromacs+"mdrun -deffnm " + output_name)
	UserDefSimulation(topology_file,input_structure,output_name,"runfile.mdp")

def RobustMolecularDynamics(topology_file="ioned_topol",input_structure="em",output_name="t_0",seed=1993,steps=100,pcoupling="semiisotropic",wallclock=24,epsilon=20,temperature=None,apolar=True,pme=False,umbrella_name=None):
	if debug_status:
		replace_steps(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","runfile.mdp")
	elif pcoupling == "semiisotropic":
		replace_seed(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","seeded.mdp",seed)
		replace_epsilon("seeded.mdp","epsilon-correct.mdp",epsilon)
		replace_temperature("epsilon-correct.mdp","retemp.mdp",temperature)
		replace_steps("retemp.mdp","runfile.mdp",steps)
		#replace_steps("seeded.mdp","runfile.mdp",steps)
	elif pcoupling == "anisotropic":
		replace_seed(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","seeded.mdp",seed)
		replace_pcoupling("seeded.mdp","aniso_pcouple.mdp",coupling="anisotropic")
		replace_epsilon("aniso_pcouple.mdp","epsilon-correct.mdp",epsilon)
		replace_temperature("epsilon-correct.mdp","retemp.mdp",temperature)
		replace_steps("retemp.mdp","runfile.mdp",steps)
		#replace_steps("aniso_pcouple.mdp","runfile.mdp",steps)
		#replace_steps("epsilon-correct.mdp","runfile.mdp",steps)
	else:
		print "Error- pcoupling type not recognised"
		sys.exit()
	if not apolar:
		print "Polar martini parameters in use: changing the appropriate parameters"
		pmartini_update("runfile.mdp")
		#This is a special change- so that all CG systems can have a common wallclock we specifically increase it for PMARTINI
		wallclock *= 4
	if pme:
		print "Turning on PME. Unknown effects. Here be dragons"
		pme_update("runfile.mdp")
		#This is a special change- so that all CG systems can have a common wallclock we specifically increase it for PMARTINI
		wallclock *= 4
		
	#os.system(gromacs+"grompp -f runfile -c " + input_structure + " -p " + topology_file + " -o " + output_name + " >& grompp_" + output_name + ".log")
	#os.system(gromacs+"mdrun -deffnm " + output_name)
	UserDefSimulation(topology_file,input_structure,output_name,"runfile",wallclock=wallclock,wallclock_lifetime=wallclock,umbrella=umbrella_name)
	
	default_timestep = 0.02
	reducer=1
	
	reduce_values = ("dt")
	increase_values = ("nstxout","nstvout","nstlog","nstenergy","nstxtcout")
	
	if output_name + ".gro" not in os.listdir("."):
			reducer *= 0.5
			#if reducer >= 0.25: 
			timestep = reducer * default_timestep
			restart_name = "restart_" + str(timestep)
			#Determine last frame
			logfile = open(output_name+".log")
			logcontents = logfile.readlines()
			for i in range(len(logcontents)):
				line = logcontents[i]
				contents = line.split()
				if len(contents) == 3:
					if contents[0] == "Step":
						if len (logcontents[i+1].split()) == 3:
							lasttime = float(logcontents[i+1].split()[1])
			logfile.close()
			#final time is the final 2-thousandth frame 
			finaltime = int(lasttime/2000)*2000 
			if finaltime < 0: finaltime = 0
			os.system("echo '0\n' | " + gromacs + "trjconv -s " + output_name + ".tpr -f " + output_name + ".trr -dump " + str(finaltime) + " -o " + restart_name + "_ini.gro")
			
			#read the runfile and write out updated lines (with smaller timestep etc) to restart.mdp
			restart_mdp = open(restart_name + ".mdp","w")
			old_mdp = open("runfile.mdp","r")
			for line in old_mdp:
				contents = line.split()
				if len(contents) == 0 or len(contents) == 1:
					print >> restart_mdp, line
				elif contents[0] in reduce_values:
					for_restart = contents[0] + " = " + str(float(contents[2]) * reducer)
					print >> restart_mdp, for_restart
				elif contents[0] in increase_values:
					for_restart = contents[0] + " = " + str(float(contents[2]) / reducer)
					print >> restart_mdp, for_restart
				elif contents[0] == "tinit":
					print >> restart_mdp, "tinit                    = " + str(finaltime)
				elif contents[0] == "nsteps":
					#The number of steps should be the new total number of steps minus the number that would have been done
					previous_steps = float(finaltime)/timestep
					new_total_steps = steps / reducer -  previous_steps
					print >> restart_mdp, "nsteps                   = " + str(new_total_steps)
				elif contents[0] == "gen_vel":
					print >> restart_mdp, "gen_vel = no"
				else:
					print >> restart_mdp, line
			restart_mdp.close()
			old_mdp.close()
			#Run
			UserDefSimulation(topology_file,restart_name+"_ini.gro",restart_name,restart_name+".mdp",umbrella=umbrella_name)
			#Remove broken frames from the trajectory, after the time point
			os.system(gromacs + "trjconv -f " + output_name + ".xtc -e "+ str(finaltime) + " -o " + output_name + ".xtc")
			#concatenated output_name.xtc and restart.xtc to output_name.xtc
			os.system(gromacs + "trjcat -f " + output_name + ".xtc " + restart_name + ".xtc -o " + output_name + ".xtc") 
			#move restart.gro to output_name.gro
			os.rename(restart_name + ".gro", output_name + ".gro")

def restartMD(topology_file="ioned_topol",input_structure="em",output_name="t_0",seed=1993,steps=100,pcoupling="semiisotropic",wallclock=1000,epsilon=20,temperature=None,apolar=True):
		'''Have we tried this before? If so cat everything and impersonate t_0'''
		if os.path.exists("crash_rerun.log"):
			os.system("cat " + output_name + ".log crash_rerun.log > cat.log")
			shutil.move("cat.log",output_name+".log")
			os.system(gromacs + "trjcat -f " + output_name + ".xtc crash_rerun.xtc -o " + output_name + ".xtc")
			
		
		
		
		'''determine the endpoint of the old file'''
		#Determine last frame
		logfile = open(output_name+".log")
		logcontents = logfile.readlines()
		for i in range(len(logcontents)):
				line = logcontents[i]
				contents = line.split()
				if len(contents) == 3:
					if contents[0] == "Step":
						if len (logcontents[i+1].split()) == 3:
							lasttime = float(logcontents[i+1].split()[1])
		logfile.close()
		
		#final time is the final 2-thousandth frame 
		finaltime = int(lasttime/2000)*2000 
		if finaltime < 0: finaltime = 0
		os.system("echo '0\n' | " + gromacs + "trjconv -s " + output_name + ".tpr -f " + output_name + ".trr -dump " + str(finaltime) + " -o " +  "crash.gro")
		timestep = 0.02
		#read the runfile and write out updated lines (with smaller timestep etc) to restart.mdp
		restart_mdp = open("crash_restart.mdp","w")
		old_mdp = open("runfile.mdp","r")
		for line in old_mdp:
				contents = line.split()
				if len(contents) == 0 or len(contents) == 1:
					print >> restart_mdp, line
				elif contents[0] == "tinit":
					print >> restart_mdp, "tinit                    = " + str(finaltime)
				elif contents[0] == "nsteps":
					#The number of steps should be the new total number of steps minus the number that would have been done
					previous_steps = float(finaltime)/timestep
					new_total_steps = float(contents[2]) -  previous_steps
					print >> restart_mdp, "nsteps                   = " + str(new_total_steps)
				elif contents[0] == "gen_vel":
					print >> restart_mdp, "gen_vel = no"
				else:
					print >> restart_mdp, line
		restart_mdp.close()
		old_mdp.close()
		
		UserDefSimulation(topology_file,"crash.gro","crash_rerun","crash_restart",wallclock=True,wallclock_lifetime=wallclock)

		os.system(gromacs + "trjconv -f " + "t_0.xtc -e "+ str(finaltime) + " -o " + output_name + ".xtc")
		#concatenated output_name.xtc and restart.xtc to output_name.xtc
		os.system(gromacs + "trjcat -f t_0.xtc crash_rerun.xtc -o t_0.xtc") 
		#move restart.gro to output_name.gro
		os.rename("crash_rerun.gro", output_name + ".gro")

def pdb2gmx(options={"f":"initial_helix.pdb", "p":"initial_helix", "o":"initial_helix", "ignh":False}, selections=["13","3"],logfile="/dev/null"):
	if options["ignh"]:
		options["ignh"] = "-ignh"
	else:
		options["ignh"] = ""
        [pdb2gmx_stdin, pdb2gmx_stout_sterr] = os.popen4(HAConf.configuration['gromacs_location'] + "pdb2gmx -f %(f)s -p %(p)s -o %(o)s %(ignh)s" % options)
        #Select gmx53a6
	for type in selections:
        	print >> pdb2gmx_stdin, type
        pdb2gmx_stdin.flush()
        with open(logfile,"w") as log:
		for line in pdb2gmx_stout_sterr: print >> log, line,

def grompp(options={"f":"quick_sd.mdp","c":"initial_helix","p":"initial_helix","o":"helix_em","maxwarn":0},logfile="/dev/null",mod=""):
	[grompp_stdin, grompp_stout_sterr] = os.popen4(HAConf.configuration['gromacs_location'] + ("grompp -f %(f)s -c %(c)s -p %(p)s -o %(o)s -maxwarn %(maxwarn)d" % options) + mod)
        with open(logfile,"w") as log:
		for line in grompp_stout_sterr: print >> log, line,

def mdrun(options={"deffnm":"helix_em"}, logfile="/dev/null"):
	[mdrun_stdin, mdrun_stout_sterr] = os.popen4(HAConf.configuration['gromacs_location'] + "mdrun -deffnm %(deffnm)s" % options)
        with open(logfile,"w") as log:
		for line in mdrun_stout_sterr: print >> log, line
