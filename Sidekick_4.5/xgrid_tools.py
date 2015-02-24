#!/usr/bin/python

import os,time,pickle
import HAConf  

controller = HAConf.xgrid['controller']
password = HAConf.xgrid['password']

class job_save:
	def __init__(self,filename):
		self.restart_filename = filename
	def update(self,job_list):
		restart_file = open(self.restart_filename,'w')
		pickle.dump(job_list,restart_file)
		restart_file.close()

def job_submit(command,target_directory=None):
	
	if target_directory == None:
		target_directory =  os.getcwd()

	if password == "":
		command = 'xgrid -job submit -h ' + controller  + ' ' + command
	else:
		command = 'xgrid -job submit -h ' + controller + " -p " + password + ' ' + command
	
	results = os.popen4(command)
	for line in results[1]:
		if line[4:17] == "jobIdentifier": job_id = int(line.split()[2][:-1])
	time.sleep(120)
	return (job_id, target_directory)
	
def batch_submit(command,target_directory=None):
	
	if target_directory == None:
		target_directory =  os.getcwd()

	if password == "":
		command = 'xgrid -job batch -h ' + controller  + ' ' + command
	else:
		command = 'xgrid -job batch -h ' + controller + " -p " + password + ' ' + command
	
	results = os.popen4(command)
	for line in results[1]:
		#print line
		if line[4:17] == "jobIdentifier": job_id = int(line.split()[2][:-1])
	return (job_id, target_directory)

class xgrid_monitor():
	def __init__(self,controller,password,job_id):
		self.controller = controller
		self.job_id = job_id
		def_command = "xgrid -job attributes -h " + controller + " -id " + str(job_id) + " -p " + password
		results = os.popen4(def_command)
		self.stdout = results[1]
		self.verbose = False
		self.process()
	def process(self):
		#by default, print every line with an attribute
		for line in self.stdout:
			if line[0] == "{" or line[0] == "}": continue
			print line.split()[0] + "\t" + line.split()[2][:-1]

class get_status(xgrid_monitor):
	def process(self):
		for line in self.stdout:
			if self.verbose: print line
			if line.split()[0] != "jobStatus": continue
			status = line.split()[2][:-1]
		#print status
		try:
			self.status = status
		except:
			self.status = "Unknown"
			self.verbose = True
	def status_report(self): return self.status
	
class get_status_nd(xgrid_monitor):
        def __init__(self,controller,password,job_id,log):
                self.controller = controller
                self.job_id = job_id
                def_command = "xgrid -job attributes -h " + controller + " -id " + str(job_id) + " -p " + password
                results = os.popen4(def_command)
                self.stdout = results[1]
                self.verbose = False
		self.log = log
                self.process()
        def status_report(self): print >> self.log, self.job_id, self.status; return self.status

def collect_results(job_list,log):
	#job_list is a list of tuples, with [0] as job_id and [1] as the target directory
	global controller
	global password
	log_location = open(log,"w")
	print >> log_location, "Started %s" % time.ctime()
	while len(job_list) > 0:
                time.sleep(3000)
                print >> log_location, "################################"
                print >> log_location, "Tick %s" % time.ctime()
                #print job_list
                #use a working list to prevent problems associated with in loop editing of the job_list
                working_list = list(job_list)
                for job in working_list:
                        time.sleep(7)
			job_id = job[0]
                        target_directory = job[1]
                        status = get_status_nd(controller, password, job_id, log_location).status_report()
                        #print job_id, status
                        if status == "Finished":
                                retrieve_command = "xgrid -job results -h " + controller + " -id " + str(job_id) + " -out " + target_directory + " -so " + target_directory + "/xgrid_log.txt -se " + target_directory + "/xgrid_err.txt -p " + password
                                os.system(retrieve_command)
                                job_list.pop(job_list.index(job))
                        elif status == "Failed":
                                job_list.pop(job_list.index(job))
                        #print "Loop", job_list
        print >> log_location, "Finished %s" % time.ctime()
	exit()

def collect_results_daemon(job_list,log_file,restart_file):
	#job_list is a list of tuples, with [0] as job_id and [1] as the target directory
	import daemonize
	global controller
	global password
	daemonize.daemonize('/dev/null',log_file,log_file)
	restart_pickle = job_save(restart_file)
	restart_pickle.update(job_list)
	print "Daemon started with pid %d" % os.getpid()
	print "Started %s" % time.ctime()
	while len(job_list) > 0:
		time.sleep(1800)
		print "################################"
		print "Tick %s" % time.ctime()
		#print job_list
		#use a working list to prevent problems associated with in loop editing of the job_list
		working_list = list(job_list)
		for job in working_list:
			job_id = job[0]
			target_directory = job[1]
			status = get_status(controller, password, job_id).status_report()
			print job_id, status, target_directory
			if status == "Finished":
				retrieve_command = "xgrid -job results -h " + controller + " -id " + str(job_id) + " -out " + target_directory + " -so " + target_directory + "/xgrid_log.txt -se " + target_directory + "/xgrid_err.txt -p " + password  
				os.system(retrieve_command)
				job_list.pop(job_list.index(job))
			elif status == "Failed":
				job_list.pop(job_list.index(job))
			#print "Loop", job_list
		restart_pickle.update(job_list)
	print "Finished %s" % time.ctime()
	exit()
