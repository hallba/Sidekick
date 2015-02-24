#!/usr/bin/env python
#ActivityRelationship.py
#This replaces consensus.py in that it takes generates a full set of consensus data from all existing data and then can attempt to relate it back to activities

import SidekickOptParse
import HAConf
from os.path import exists
import sys, math, os

def read_sidekick_listfile(filename):
	input_file = open(filename,'r')
	data = []
	for line in input_file:
		if line[0] == "@" or line[0] == "#": continue
		[time, current_data] = line.split()
		data.append(float(current_data))
	input_file.close()
	return data
		
class Datastore:
	'''A slightly more refined class for keeping data and associated analyses on said data'''
	def __init__(self,content="Distance",units="Angstroms"):
		self.surface = []
		self.center = []
		self.content = content
		self.units = units
		self.mean_surface = "NA"
		self.mean_center = "NA"
		self.sd_surface = "NA"
		self.sd_center  = "NA"
	def statistics(self):
		self.means()
		self.standard_deviations()
	def means(self):
		if len(self.surface) > 0:
			self.mean_surface = mean(self.surface)
		if len(self.center) > 0:
			self.mean_center = mean(self.center)
	def standard_deviations(self):
		if len(self.surface) > 0:
			self.sd_surface = (mean([(datapoint - self.mean_surface)**2 for datapoint in self.surface]))**0.5
		if len(self.center) > 0:
			self.sd_center = (mean([(datapoint - self.mean_center)**2 for datapoint in self.center]))**0.5
		
class AngularDatastore(Datastore):
	def means(self):
		if len(self.surface) > 0:
			self.mean_surface = angle_mean_degree(self.surface)
		if len(self.center) > 0:
			self.mean_center = angle_mean_degree(self.center)
	def standard_deviations(self):
		#angular deviations formally. From Batchelet 1981
		def CalculateVariance(data):
			R = CalculateMeanVectorLength(data)
			return 2*(1 - R)
		def CalculateMeanVectorLength(data):
			#Kept as a seperate function in case I ever want to use it again
			return ((sum([math.cos(degree_to_radian(angle)) for angle in data])**2 + sum([math.sin(degree_to_radian(angle)) for angle in data])**2)**0.5)/len(data)
		#r2 = sum([math.cos(degree_to_radian(angle)) for angle in self.surface])**2 + sum([math.sin(degree_to_radian(angle)) for angle in self.surface])**2
		#R = (r2 ** 0.5)/len(self.surface)
		if len(self.surface) > 0:
			self.sd_surface = CalculateVariance(self.surface) ** 1/2
		#r2 = sum([math.cos(degree_to_radian(angle)) for angle in self.center])**2 + sum([math.sin(degree_to_radian(angle)) for angle in self.center])**2
		#R = (r2 ** 0.5)/len(self.surface)
		if len(self.center) > 0:
			self.sd_center = CalculateVariance(self.center) ** 1/2

def mean(number_list):
	return sum(number_list)/len(number_list)
	
def degree_to_radian(angle):
	return angle*math.pi/180
	
def radian_to_degree(angle):
	return angle*180/math.pi

def angle_mean_degree(angle_list):
	#convert to radians
	sin_mean = mean([math.sin(degree_to_radian(angle)) for angle in angle_list])
	cos_mean = mean([math.cos(degree_to_radian(angle)) for angle in angle_list])
	mean_angle = math.atan2(sin_mean,cos_mean)
	return radian_to_degree(mean_angle)

def consensus_search(starting_point):
	#Main loop
	for forcefield in os.listdir(starting_point):
		for sequence in os.listdir(starting_point + "/" + forcefield):
			positions = datastore()
			rotations = datastore()
			tilts = datastore()
			kinks = datastore()
			
			for run in os.listdir(starting_point + "/" + forcefield + "/" + sequence):
				run_directory = starting_point + "/" + forcefield + "/" + sequence + "/" + run
				#Read in the list-style data files
				position_data = read_sidekick_listfile(run_directory+"/bilayer_position.dat")
				rotation_data = read_sidekick_listfile(run_directory+"/rot.dat")
				tilt_data = read_sidekick_listfile(run_directory+"/bun_tilt.xvg")
				kink_data = read_sidekick_listfile(run_directory+"/bun_kink.xvg")
				#
				for timepoint in zip(position_data,rotation_data,tilt_data,kink_data):
					if timepoint[0] >= 10 or timepoint[0] <= -10:
						#surface_data
						if timepoint[0] > 0:
							positions.surface.append(timepoint[0])
						else:
							positions.surface.append(-timepoint[0])
						rotations.surface.append(timepoint[1])
						if tilt <= 90:
							tilts.surface.append(timepoint[2])
						else:
							tilts.surface.append(180 - timepoint[2])
						kink_data.surface.append(timepoint[3])
					else:
						#center_data
						positions.center.append(timepoint[0])
						rotations.center.append(timepoint[1])
						if tilt <= 90:
							tilts.center.append(timepoint[2])
						else:
							tilts.center.append(180 - timepoint[2])
						kink_data.center.append(timepoint[3])
						
if __name__ == "__main__":
	parser = SidekickOptParse.Analysis_Options(("-j",))
	
	#Initialise
	
	starting_point = HAConf.results_location + "/" + parser.options.destination
	
	#Is the job name sensible?
	if not exists(starting_point):
		print "Job name does not correspond with an existing job directory. Exiting"
		print "Failure to read", starting_point
		sys.exit()

	if parser.options.summary_file != None:
		write_to_file = True
		csv_output = open(parser.options.summary_file,"w")
	else:
		write_to_file = False
	
	
	results = {}
	metric_list = ["Position","Tilt","Rotation","Kink"]
	
	for forcefield in os.listdir(starting_point):
		results[forcefield] = {}
		for sequence in os.listdir(starting_point + "/" + forcefield):
			#Now start collecting data for each sequence
			results[forcefield][sequence] = {	"Position"	:	Datastore(content="Distance From Bilayer Center",units="Angstroms"),
												"Tilt"		:	AngularDatastore(content="Tilt",units="Degrees"),
												"Rotation"	:	AngularDatastore(content="Rotation",units="Degrees"),
												"Kink"		:	AngularDatastore(content="Kink",units="Degrees")
			}
			for run in os.listdir(starting_point + "/" + forcefield + "/" + sequence):
				run_directory = starting_point + "/" + forcefield + "/" + sequence + "/" + run
				position_data = read_sidekick_listfile(run_directory+"/bilayer_position.dat")
				rotation_data = read_sidekick_listfile(run_directory+"/rot.dat")
				tilt_data = read_sidekick_listfile(run_directory+"/bun_tilt.xvg")
				kink_data = read_sidekick_listfile(run_directory+"/bun_kink.xvg")
				#
				for timepoint in zip(position_data,rotation_data,tilt_data,kink_data):
					if timepoint[0] >= 10 or timepoint[0] <= -10:
						#surface_data- negative vs positive membrane position is unimportant
						if timepoint[0] > 0:
							results[forcefield][sequence]["Position"].surface.append(timepoint[0])
						else:
							results[forcefield][sequence]["Position"].surface.append(-timepoint[0])
						results[forcefield][sequence]["Rotation"].surface.append(timepoint[1])
						if timepoint[2] <= 90:
							results[forcefield][sequence]["Tilt"].surface.append(timepoint[2])
						else:
							results[forcefield][sequence]["Tilt"].surface.append(180 - timepoint[2])
						results[forcefield][sequence]["Kink"].surface.append(timepoint[3])
					else:
						#center_data
						results[forcefield][sequence]["Position"].center.append(timepoint[0])
						results[forcefield][sequence]["Rotation"].center.append(timepoint[1])
						if timepoint[2] <= 90:
							results[forcefield][sequence]["Tilt"].center.append(timepoint[2])
						else:
							results[forcefield][sequence]["Tilt"].center.append(180 - timepoint[2])
						results[forcefield][sequence]["Kink"].center.append(timepoint[3])
			#Get the stats from the run data
			for metric in results[forcefield][sequence].keys():
				results[forcefield][sequence][metric].statistics()

	#Finally, print all the results out in a nice looking format
	if write_to_file:
		data_line = "Forcefield,Sequence,Position,"
		for metric in metric_list:
			data_line += metric + "," + metric + "sd,"
		print >> csv_output, data_line
	for forcefield in results.keys():
		print "-"*60
		print "-"*60
		print "Forcefield: " + forcefield
		for sequence in results[forcefield].keys():
			print "-"*60
			print "\tSequence: " + sequence
			print "\tCenter"
			for metric in results[forcefield][sequence].keys():
					print "\t" + metric + "\t" + str(results[forcefield][sequence][metric].mean_center) + " sd " + str(results[forcefield][sequence][metric].sd_center)
			print "\tSurface"
			for metric in results[forcefield][sequence].keys():
					print "\t" + metric + "\t" + str(results[forcefield][sequence][metric].mean_surface) + " sd " + str(results[forcefield][sequence][metric].sd_surface)
			if write_to_file:
				data_line = "%s,%s,Surface," % (forcefield,sequence)
				for metric in metric_list:
					data_line += str(results[forcefield][sequence][metric].mean_surface) + ","
					data_line += str(results[forcefield][sequence][metric].sd_surface) + ","
				print >> csv_output, data_line
				data_line = "%s,%s,Center," % (forcefield,sequence)
				for metric in metric_list:
					data_line += str(results[forcefield][sequence][metric].mean_center) + ","
					data_line += str(results[forcefield][sequence][metric].sd_center) + ","
				print >> csv_output, data_line
				csv_output.flush()
	csv_output.close()