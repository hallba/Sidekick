#!/usr/bin/python

from __future__ import with_statement

#Batch tools- reduce the amount of code in top level scripts
import HAConf
import os,sys,time,shutil,pickle,filecmp
from random import randint

def blackbox(options,name="blackbox.pickle"):
	with open(name,"w") as bl_fh:
		pickle.dump(options,bl_fh)
	
class batch_event:
	def __init__(self,options,temp_folder_location=None,original_folder_location=None,networked_tmp_location=True):
		self.temp_folder_name=temp_folder_location
		#print self.temp_folder_name
		#generate a folder location
		initial_sequence=options.sequence
		seed = options.seed
		MDsteps = options.simlength*1000./0.02
		model_type = options.model_type
		self.options = options
		
		self.seed = "%09d" % seed
		self.sequence = initial_sequence
		self.nfs_tmp = networked_tmp_location
	
		#Test if this has been run previously-if so exit (to restart batch jobs)
	
		translation = int(options.position*10)
		if options.bias:
			model_type_name = model_type
		else:
			model_type_name = model_type + "_unbias"
		hg_mix = options.headgroup_mix.split("/")
		#model_type_name += "_" + hg_mix[0] + "-" + hg_mix[1]
		if len(hg_mix) > 1 :
			lipid_type_name = "/" + hg_mix[0] + "/" + hg_mix[1]
		else:
			lipid_type_name = "/" + hg_mix[0] + "/" + "1"
		if options.special:
			model_type_name += "_" + options.special
		if options.temperature:
			temperature = options.temperature
		else:
			temperature = 323
		
		seed_name = "%09d-%03d-%+03d-%03d" % (seed,options.angle,translation,temperature)
		self.directory_tree = HAConf.results_location + "/" + options.destination + "/" + model_type_name + "/" + initial_sequence + "/" + lipid_type_name + "/" + seed_name
		
		if original_folder_location == None:
			self.original_location = os.getcwd()
		else:
			self.original_location = original_folder_location
			
		self.go()
	def tmp_location(self,temp_mod):
		if self.nfs_tmp:
			return HAConf.results_location + "/tmp/" + "Sidekick_temp_%09d" % (randint(0,999999999)) + temp_mod
		else:
			return "./Sidekick_temp_%09d" % (randint(0,999999999)) + temp_mod

class helix_batch_enter(batch_event):
	def go(self):
		if os.path.exists(self.directory_tree):
			print "Already done this one; must be a restart following a crash.\nExiting after short, random sleep to avoid overloading controller"
			time.sleep(30+randint(0,15)*20)
			sys.exit()
		
	
		#New batch mode works in two parts- the first creates a temporary directory and goes there to do work
		##At the end, it moves the files to the final location
		#from random import randint
		
		#New additions- now write to a temporary folder on the shared disk space, before copying the complete system (if everything has worked)
		#This kind of reverses a previous design decision to write to the shared disks only when the simulation is complete. It may be reverted
		#Master temporarylocation is a folder called tmp in the data folder
		
		if not os.path.exists(HAConf.results_location):
			os.mkdir(HAConf.results_location)
			os.system("chmod a+rwx "+HAConf.results_location)
		
		if not os.path.exists(HAConf.results_location + "/tmp"):
			os.mkdir(HAConf.results_location + "/tmp")
			os.system("chmod a+rwx "+HAConf.results_location  + "/tmp")		
		
		temp_mod = "_" + os.uname()[1].split(".")[0] + "_" + self.sequence + "_" + self.seed
		
		self.temp_folder_name = self.tmp_location(temp_mod) 
		
		while(os.path.exists(self.temp_folder_name)):
			self.temp_folder_name = self.tmp_location(temp_mod) 
		
		try:
			os.mkdir(self.temp_folder_name)
		except:
			print "Cannot create temporary folder. Exiting"
			sys.exit()
		os.chdir(self.temp_folder_name)
			
		#write out a pickle containing the original input options
		blackbox(self.options)
		
		#write logs out to the new directory, rather than stdout
		stdout = os.getcwd() + "/" + "Sidekick_log.txt"
		stderr = os.getcwd() + "/" + "Sidekick_err.txt"
		for f in sys.stdout, sys.stderr: f.flush()
		so = file(stdout, 'a+')
		se = file(stderr, 'a+',0)
		os.dup2(so.fileno(),sys.stdout.fileno())
		os.dup2(se.fileno(),sys.stderr.fileno())
		print os.uname()[1]

class helix_batch_exit(batch_event):
	def go(self):
		os.chdir("..")
		try:
			shutil.rmtree(self.original_location + "/.matplotlib")
		except:
			pass
			
		#Ok, try 30 times over the course of 5 hours to move the files. If unsuccessful, to prevent data induced XGrid crashes, delete the files and free the node
		print "Moving", self.temp_folder_name, "to", self.directory_tree
		for i in range(30):
			
			try:
				print "Attempting move..."
				if os.path.exists(self.directory_tree):
					#delete any failed copy events
					#Something complex happened whereby a move was somehow interupted
					print "Previous error created a results directory. Comparing and updating"
					sys.exit()
					#shutil.rmtree(self.directory_tree)
				#shutil.move(self.temp_folder_name,self.directory_tree)
				shutil.copytree(self.temp_folder_name,self.directory_tree)
				print "Files copied. Removing temp folder"
				#shutil.rmtree(self.temp_folder_name)
				os.system("rm -rf " + self.temp_folder_name)
				print "Successful Completion"
				break
			except:
				print "Unexpected error:", sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2], 
				if not os.path.exists(self.temp_folder_name):
					#folder has gone: break and finish up
					print "Directory gone: finishing"
					break
				time.sleep(60)
		try:
			shutil.rmtree(self.temp_folder_name)
		except:
			pass
		
		#Not sure I need to do this- hoping that it will mean that I don't end up copying the files back to the controller after all...
		os.chdir(self.original_location)

		sys.exit()
		
class helix_dimer_batch_event:
	def __init__(self,options,temp_folder_location=None,original_folder_location=None,networked_tmp_location=True):
		self.temp_folder_name=temp_folder_location
		#print self.temp_folder_name
		#generate a folder location
		initial_sequence = ("%-040s-%-040s" % (options.sequence1, options.sequence2)).replace(" ","_")
		seed = options.seed
		#MDsteps = options.simlength*1000./0.02
		model_type = options.model_type
		self.options = options
		
		self.seed = "%09d" % seed
		self.sequence = initial_sequence
		
		self.nfs_tmp = networked_tmp_location
	
		translation = int(0) #options.position*10)
		angle = 0 #options.angle
		
		model_type_name = model_type

		hg_mix = ["DPPC",]# options.headgroup_mix.split("/")
		#model_type_name += "_" + hg_mix[0] + "-" + hg_mix[1]
		if len(hg_mix) > 1 :
			lipid_type_name = "/" + hg_mix[0] + "/" + hg_mix[1]
		else:
			lipid_type_name = "/" + hg_mix[0] + "/" + "1"
		if options.special:
			model_type_name += "_" + options.special
		if options.temperature:
			temperature = options.temperature
		else:
			temperature = 323
		
		seed_name = "%09d-%03d-%+03d-%03d" % (seed,angle,translation,temperature)
		self.directory_tree = HAConf.results_location + "/" + options.destination + "/" + model_type_name + "/" + initial_sequence + "/" + lipid_type_name + "/" + seed_name
		
		if original_folder_location == None:
			self.original_location = os.getcwd()
		else:
			self.original_location = original_folder_location
		
		self.go()

class helix_dimer_batch_enter(helix_dimer_batch_event,helix_batch_enter):
	pass

class helix_dimer_batch_exit(helix_dimer_batch_event,helix_batch_exit):
	pass	

class helix_oligomer_batch_event:
	def __init__(self,options,temp_folder_location=None,original_folder_location=None,networked_tmp_location=True):
		self.temp_folder_name=temp_folder_location
		#print self.temp_folder_name
		#generate a folder location
		initial_sequence = "" #("%-040s-%-040s" % (options.sequence1, options.sequence2)).replace(" ","_")
		for sequence in options.sequences.split(","):
			initial_sequence += ("%-040s" % sequence).replace(" ","_")
		seed = options.seed
		#MDsteps = options.simlength*1000./0.02
		model_type = options.model_type
		self.options = options
		
		self.seed = "%09d" % seed
		self.sequence = initial_sequence
		
		self.nfs_tmp = networked_tmp_location
	
		translation = int(0) #options.position*10)
		angle = 0 #options.angle
		
		model_type_name = model_type

		hg_mix = ["DPPC",]# options.headgroup_mix.split("/")
		#model_type_name += "_" + hg_mix[0] + "-" + hg_mix[1]
		if len(hg_mix) > 1 :
			lipid_type_name = "/" + hg_mix[0] + "/" + hg_mix[1]
		else:
			lipid_type_name = "/" + hg_mix[0] + "/" + "1"
		if options.special:
			model_type_name += "_" + options.special
		if options.temperature:
			temperature = options.temperature
		else:
			temperature = 323
		
		seed_name = "%09d-%03d-%+03d-%03d" % (seed,angle,translation,temperature)
		self.directory_tree = HAConf.results_location + "/" + options.destination + "/" + model_type_name + "/" + initial_sequence + "/" + lipid_type_name + "/" + seed_name
		
		if original_folder_location == None:
			self.original_location = os.getcwd()
		else:
			self.original_location = original_folder_location
		
		self.go()

class helix_oligomer_batch_enter(helix_oligomer_batch_event,helix_batch_enter):
	pass

class helix_oligomer_batch_exit(helix_oligomer_batch_event,helix_batch_exit):
	pass

class unclean_exit(helix_batch_exit):
	def go(self):
		os.chdir("..")
		try:
			shutil.rmtree(self.original_location + "/.matplotlib")
		except:
			pass
		if not os.path.exists(HAConf.results_location + "/failed"):
			os.mkdir(HAConf.results_location + "/failed")
		try:
			shutil.move(self.temp_folder_name,HAConf.results_location + "/failed/")
		except:
			pass
			#shutil.rmtree(temp_folder_name)
		finally:
			pass
		sys.exit()

class unclean_exit_dimer(helix_dimer_batch_event,unclean_exit):
	pass
	
def cleanup(temporary_location,original_location):
	'''simple function to ensure the deletion of temporary writing location and .matplotlib'''
	try:
		shutil.rmtree(original_location + "/.matplotlib")
	except:
		pass
	if not os.path.exists(HAConf.results_location + "/failed"):
		os.mkdir(HAConf.results_location + "/failed")
		os.system("chmod a+rwx "+HAConf.results_location  + "/failed")
	try:
		try:
			shutil.move(temporary_location,HAConf.results_location + "/failed/"+temporary_location)
			os.rmdir(temporary_location)
		except:
			shutil.rmtree(temporary_location)
			os.rmdir(temporary_location)
		finally:
			os.rmdir(temporary_location)
	except:
		pass