#!/usr/bin/python

from __future__ import with_statement
import os,pickle,math,numpy,sys
import HAConf, AutomatedPlot

job_pickle = "/tmp/sidekick.pickle"
top_dir = HAConf.results_location
vis_loc = HAConf.visualisation_location
pub_loc = HAConf.public_location

block_number = 5

class dataset:
	def __init__(self,sequence,thickness=41):
		self.position = []
		self.tilt = []
		self.rotation = []
		self.sequence = sequence
		self.thickness = thickness
		self.quadrupolar_splitting = [[] for item in sequence]
	def analyse(self,xml_name):
		try:
			#for metric in (self.position, self.tilt, self.rotation):
			position_hist, position_bins = numpy.histogram(self.position,40,(-10,30))			
			tilt_hist, tilt_bins = numpy.histogram(self.tilt,45,(0,90))
			rotation_hist, rotation_bins = numpy.histogram(self.rotation,24,(-180,180))
			#qs_hist = [numpy.histogram(item,50,(0,50)) for item in self.quadrupolar_splitting]
			qs_averages= [sum(item)/len(item) for item in self.quadrupolar_splitting]
			#print len(qs_hist), len(qs_hist[0]), len(qs_hist[0][0])
			#print qs_hist[0][0]
			#print qs_hist[0][1]
			
			#2D
			pvt_hist, pvt_pbins, pvt_tbins = numpy.histogram2d(self.position,self.tilt,(40,45),((-10,30),(0,90)))
			pvr_hist, pvr_pbins, pvr_rbins = numpy.histogram2d(self.position,self.rotation,(40,24),((-10,30),(-180,180)))
			tvr_hist, tvr_tbins, tvr_rbins = numpy.histogram2d(self.tilt,self.rotation,(45,24),((0,90),(-180,180)))
	
			#3D
			self.descriptors = numpy.array(zip(self.position,self.tilt,self.rotation))
			pvtvr_hist, pvtvr_edges = numpy.histogramdd(self.descriptors,bins=(40,45,24),range=((-10,30),(0,90),(-180,180) ) )
	
			print position_hist, len(position_hist)
			print tilt_hist, len(tilt_hist)
			print rotation_hist, len(rotation_hist)
			with open(xml_name,"w") as oup:
				print >> oup, "<dataset>"
				print >> oup, "<sequence>" + self.sequence + "</sequence>"
				print >> oup, "<thickness>" + str(self.thickness) + "</thickness>"
				self.print_data(oup,"position",position_hist,position_bins)
				self.print_data(oup,"tilt",tilt_hist,tilt_bins)
				self.print_data(oup,"rotation",rotation_hist,rotation_bins)
				#for resid in range(len(qs_hist)):
					#self.print_data(oup,"quadrupolar_splitting_r_"+str(resid+1),qs_hist[resid][0],qs_hist[resid][1])
				self.print_list_data(oup,"quadrupolar_splitting_mean",qs_averages)
				self.print_2d_data(oup,"pvt",pvt_hist,pvt_pbins,pvt_tbins)
				self.print_2d_data(oup,"pvr",pvr_hist,pvr_pbins,pvr_rbins)
				self.print_2d_data(oup,"tvr",tvr_hist,tvr_tbins,tvr_rbins)
				self.print_3d_data(oup,"pvtvr",pvtvr_hist,pvtvr_edges)
				print >> oup, "</dataset>"
		except:
			print "Unexpected error:", sys.exc_info()[0]
	def print_data(self,oup,name,hist,bins):
                        print >> oup, "<"+name+">",
                        datastring = ""
                        for item in hist: datastring += str(item) + ","
                        print >> oup, datastring[:-1],
                        print >> oup, "</"+name+">"
                        print >> oup, "<"+name+"_bin>",
                        datastring = ""
                        for item in bins: datastring += str(item) + ","
                        print >> oup, datastring[:-1],
                        print >> oup, "</"+name+"_bin>"
	def print_list_data(self,oup,name,hist):
                        print >> oup, "<"+name+">",
                        datastring = ""
                        for item in hist: datastring += str(item) + ","
                        print >> oup, datastring[:-1],
                        print >> oup, "</"+name+">"
	def print_2d_data(self,oup,name,hist,xbins,ybins):
			print >> oup, "<"+name+">",
			datastring = ""
			for row in hist:
				for item in row:
					datastring += str(item) + ","
				datastring = datastring[:-1] + "/"
			print >> oup, datastring[:-1],
			print >> oup, "</"+name+">"
			print >> oup, "<"+name+"_xbin>",
			datastring = ""
			for item in xbins: datastring += str(item) + ","
			print >> oup, datastring[:-1],
			print >> oup, "</"+name+"_xbin>"
			print >> oup, "<"+name+"_ybin>",
			datastring = ""
			for item in ybins: datastring += str(item) + ","
			print >> oup, datastring[:-1],
			print >> oup, "</"+name+"_ybin>"
	def print_3d_data(self,oup,name,hist,edges):
			print >> oup, "<"+name+">",
			datastring = ""
			for row in hist:
				for column in row:
					for item in column:
						datastring += str(item) + ","
					datastring = datastring[:-1] + "/"
				datastring = datastring[:-1] + ":"
			print >> oup, datastring[:-1],
			print >> oup, "</"+name+">"
 			print >> oup, "<"+name+"_xbin>",
			datastring = ""
			for item in edges[0][:-1]: datastring += str(item) + ","
			print >> oup, datastring[:-1],
			print >> oup, "</"+name+"_xbin>"
			print >> oup, "<"+name+"_ybin>",
			datastring = ""
			for item in edges[1][:-1]: datastring += str(item) + ","
			print >> oup, datastring[:-1],
			print >> oup, "</"+name+"_ybin>"
			print >> oup, "<"+name+"_zbin>",
			datastring = ""
			for item in edges[2][:-1]: datastring += str(item) + ","
			print >> oup, datastring[:-1],
			print >> oup, "</"+name+"_zbin>"

def stats(some_list):
	if len(some_list) == 0:
		return [0,0,0]
	list_mean = mean(some_list)
	list_sd = sample_sd(some_list,list_mean)
	list_abdev = mean_abs_dev(some_list,list_mean)
	return [list_mean, list_sd, list_abdev]

def mean(some_list):
	return sum(some_list)/len(some_list)
def sample_sd(some_list,mean):
	return (sum([(item - mean)**2 for item in some_list])/(len(some_list)-1))**0.5
def mean_abs_dev(some_list,mean):
	return sum([math.fabs(item - mean) for item in some_list])/len(some_list)

def get_efficiency(run_directory):
	'''This returns a limited set of results corresponding to what is found
	[flag, percentage]
	0 = file not found or otherwise abnormal-ignore
	1 = not inserted for > 50% of the time
	2 = inserted for > 50 % of the time
	'''
	try:
		run_output = open(run_directory+"/efficiency.dat")
	except:
		#print "Cannot read file", run_directory+"/efficiency.dat"
		return [0, None]
	for line in run_output:
		if line[:5] == "Numbe":
			contents = line.split()
			if int(contents[-1]) != 1:
				#print "Wrong number of bilayers..."
				return [0, None]
		elif line[:5] == "Unusu":
			#print "Inappropriate mixing"
			return [0, None]
		elif line[:5] == "Helix":
			#print line
			contents = line.split()
			#print contents[-2]
			if float(contents[-2]) > 50:
				return [2, float(contents[-2]) ]
			else:
				return [1, float(contents[-2]) ]
		else:
			print "Error- no case is satisfied"
			exit()

def subdir_contents(original_directory,depth):
	'''Returns all the directories in subdirectories at a given depth. 0 returns the directories contents'''
	curr_dir = [original_directory,]
	for level in range(depth):
		#curr_dir = [fullpath_listdirindir(item) for item in curr_dir]
		new_dir = []
		for item in curr_dir:
			#print "Item, fullpath", item, fullpath_listdirindir(item)
			new_dir += fullpath_listdirindir(item)
		curr_dir = new_dir
		#print curr_dir
	return curr_dir
		
def fullpath_listdir(directory):
	return [directory + "/" + item for item in os.listdir(directory)]

def fullpath_listdirindir(directory):
	current = []
	for item in os.listdir(directory):
		if os.path.isdir(directory + "/" + item):
			current.append(directory + "/" +item)
	return current

def backup_file(filename):
	if os.path.exists(filename):
		#target_name = "#" + filename
		last_name = filename.split("/")[-1]
		target_name = filename[:-len(last_name)] + "#" + last_name
		#print target_name
		failure = True
		#if not os.path.exists(target_name):
		os.rename(filename,target_name)
		failure = False
		#else:
		#	for i in range(20):
		#		alt_target_name = target_name + "." + str(i)
		#		if os.path.exists(alt_target_name):
		#			continue
		#		else:
		#			os.rename(filename,alt_target_name)
		#			failure = False
		#			break
		if failure:
			print "Too many backups. Clean up and try again"
			exit()

#Generate visualisation files

lipid_thickness = {	"DP" : 41,
					"DL" : 33,
					"DO" : 45,
					"PO" : 43 }

def deg2rad(angle):
	return angle/180*math.pi

def quadrupolar_calculation(tilt,rotation,resnumber=0):
	#transform rotation into common coordinates with roger; this is +200 (for G0) -180 (for different frame of reference)
	rotation = (20 + rotation )%360
	
	K = 49 #Hz
	Epar = 58.9 #degrees, taken from ideality
	Eper = -53.2 #degrees, from ideality
	delta = (rotation + Eper +200 - 100*(resnumber-1))%360 #delta wrt residue 0 (the first residue), so + 200 (rotation derived from 3)
	#print delta
	tilt = deg2rad(tilt)
	rotation = deg2rad(rotation)
	delta = deg2rad(delta)
	Epar = deg2rad(Epar)
	Eper = deg2rad(Eper)
	qs = K*0.5*0.75*(3*(math.cos(Epar * (math.cos(tilt)-math.sin(tilt)*math.cos(delta)*math.tan(Epar))**2))**2-1)
	return abs(qs)

def visgen(forcefield,lipid,sequence,mix,position_plots,rot_plots,tilt_plots,results,mix_dir):				 			
							position_file = job_plots + "/" + forcefield + "-" + sequence + "-" + lipid + "-" + mix + "-" + "bilayer_position.dat"
							tilt_file = job_plots + "/" + forcefield + "-" + sequence + "-" + lipid + "-" + mix + "-" + "tilt.dat"
							rot_file = job_plots + "/" + forcefield + "-" + sequence + "-" + lipid + "-" + mix + "-" + "rot.dat"
							kink_file = job_plots + "/" + forcefield + "-" + sequence + "-" + lipid + "-" + mix + "-" + "kink.dat"
							internal_tilt_file = job_plots + "/" + forcefield + "-" + sequence + "-" + lipid + "-" + mix + "-" + "internal_tilt.dat"
							detail_directory = job_plots + "/" + forcefield + "-" + sequence + "-" + lipid + "-" + mix
						
							for item in (position_file,tilt_file,rot_file,kink_file,internal_tilt_file):
								backup_file(item)
								os.system("touch " + item)	
							position_plots.append(position_file)
							rot_plots.append(rot_file)
							tilt_plots.append(tilt_file)
							tilt_plots.append(kink_file)
							tilt_plots.append(internal_tilt_file)

							inserted = 0.
							successful_runs = 0.
							percentage_insertion = 0.
							frame_wise_pi = []

				 			#And calculate a summary
				 			for seed_stamp in os.listdir(mix_dir):
				 					current_dir = mix_dir + "/" + seed_stamp
				 					
									[efficiency,percent] = get_efficiency(current_dir)
									print job, forcefield, sequence, lipid, mix, seed_stamp, efficiency
									if efficiency == 0:
										continue
									elif efficiency == 2:
										inserted += 1
									
									#Cat files for "good" simulations
									#print "cat " + current_dir + "/bilayer_position.dat | grep -v \#  | grep -v @ >> " + position_file
									os.system("cat " + current_dir + "/bilayer_position.dat | grep -v \#  | grep -v @ >> " + position_file)
									os.system("cat " + current_dir + "/bun_tilt.xvg | grep -v \#  | grep -v @ >> " + tilt_file)
									os.system("cat " + current_dir + "/bun_kink.xvg | grep -v \#  | grep -v @ >> " + kink_file)
									os.system("cat " + current_dir + "/rot.dat | grep -v \#  | grep -v @ >> " + rot_file)
									os.system("cat " + current_dir + "/tilt.dat | grep -v \#  | grep -v @ >> " + internal_tilt_file)
									
									successful_runs += 1
									percentage_insertion += percent
									#Assuming default temp of 323...
									#Skip simulations which are 100% or 0%...
									#if percent == 100 or percent == 0:
									#	continue
									#print percent
									#print -8.314 * 323 * math.log(percent/(100-percent))
									#frame_wise_pi.append( -8.314 * 323 * math.log(percent/(100-percent) ) )
									frame_wise_pi.append(percent)

							#Generate single file dataset and xml
							try:
								os.system("mkdir " + detail_directory)
							except:
								pass							
							print "join " + position_file + " " + internal_tilt_file + " > " + detail_directory +"/pos_tilt.dat"
							os.system("join " + position_file + " " + internal_tilt_file + " > " + detail_directory +"/pos_tilt.dat")
							os.system("join " + detail_directory + "/pos_tilt.dat " + rot_file + " > " + detail_directory +"/pos_tilt_rot.dat")
							
							thickness = lipid_thickness[lipid[:2]]
							
							interfacial = dataset(sequence,thickness)
							transmembrane = dataset(sequence,thickness)
							complete = dataset(sequence,thickness)
									
									
							with open(detail_directory +"/pos_tilt_rot.dat") as inp:
									for line in inp:
											contents = line.split()
											#print contents
											position = float(contents[1])
											tilt = float(contents[2])
											rotation = float(contents[3])

											#process 
											if tilt > 90:
												tilt = 180 - tilt
											if position < -10:
												position = -position
											#classify
											if tilt < 65 and position < 10:
												#transmembrane
												set = transmembrane
											else:
												set = interfacial
											for item in (set,complete):
												item.position.append(position)
												item.tilt.append(tilt)
												item.rotation.append(rotation)
												qs_11 = quadrupolar_calculation(tilt,rotation,11)
												for residue in range(len(item.quadrupolar_splitting)):
													qs = quadrupolar_calculation(tilt,rotation,residue)
													item.quadrupolar_splitting[residue].append(qs)
											#print "DEBUG-T/R/qs11", tilt, rotation, qs_11
									
							#generate xml files
							complete.analyse(detail_directory + "/complete.xml")
							transmembrane.analyse(detail_directory + "/transmembrane.xml")
							interfacial.analyse(detail_directory + "/interfacial.xml")
							os.system("cp " + HAConf.sidekick_location + "/webpages/base.html " + detail_directory + "/base.html")

							print results.keys()
							#block analysis of the energies; split into 10 chunks and calculate average and SD, then SE (SD/sqrt(n))
							if len(frame_wise_pi) < 30:
								#No point in trying this with too few simulations performed
								SE_Energy = "TLD"
								Mean_Energy = "TLD"
							else:
								blocks = [ [] for i in range(block_number) ]
								cblock_number = block_number
								for i in range(len(frame_wise_pi)): blocks[i%block_number].append(frame_wise_pi[i])
								#energy_statistic = stats([mean(item) for item in blocks])
								tmp_pmean =  [mean(item) for item in blocks]
								for_del = []
								for i in range(len(tmp_pmean)):
																if tmp_pmean[i] >= 100:
																	for_del.append(i)
																	cblock_number -= 1
								for_del.sort(reverse=True)
								for item in for_del:
									print item, tmp_pmean[item]
									print len(blocks), len(tmp_pmean)
									del blocks[item]
								tmp_emean =  [-8.314 * 323 * math.log( mean(item) / (100-mean(item)) ) / 1000 * 0.239  for item in blocks]
								
								print tmp_pmean
								print tmp_emean
								
								energy_statistic = stats([-8.314 * 323 * math.log( mean(item) / (100-mean(item)) ) / 1000 * 0.239  for item in blocks])
								
								Mean_Energy = energy_statistic[0]
								SE_Energy = energy_statistic[1]/(cblock_number**0.5)
								#print frame_wise_pi
								#print blocks
								#print [mean(item) for item in blocks]
								print "Energy", Mean_Energy, SE_Energy
							if successful_runs > 0:
									results[forcefield][sequence][lipid][mix] = (inserted*100/successful_runs, percentage_insertion/successful_runs, inserted, successful_runs,Mean_Energy,SE_Energy)
							else:
									results[forcefield][sequence][lipid][mix] = ("NA", inserted, successful_runs)

def webpage_header(page_handler):
	print >> page_handler, """<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">

<html>
<head>
<title>Sidekick @ University of Oxford</title>
<script src="/Scripts/dragtable.js"></script>
<script LANGUAGE="JavaScript" src="/Scripts/SidekickVisualise.js"></script>
<meta name="viewport" content="width=device-width; initial-scale=1.0; maximum-scale=1.0;">
<link rel="apple-touch-icon-precomposed" href="../../images/sidekick_ruby.png"/>
<link rel="stylesheet" type="text/css"
              href="../../Sidekick.css"
              title="standard" >
</head>
"""

###################################################
									
#Startup- Check for information gathered previously

if os.path.exists(job_pickle):
	with (open(job_pickle,'r')) as input_pickle:
		job_atime = pickle.load(input_pickle)
else:
	job_atime = {}

#Now- search each job directory, characterise it as new, old or dimer (or not interesting) and perform relevant analysis

for job in os.listdir(top_dir):
	job_dir = top_dir + "/" + job
	if job[0] == "." or not os.path.isdir(job_dir) or job == "tmp":
		continue
	with open(job_pickle,"w") as pickle_location:
		pickle.dump(job_atime,pickle_location)

	if os.path.isdir(job_dir) != True or job == ".TemporaryItems" or job == "tmp" or job == "failed":
		continue

	
	#Test job type
	print "Examining", job
	if os.path.isdir(fullpath_listdir((subdir_contents(job_dir,3)[0]))[0]):
		job_type = "new"
		times = [os.stat(item).st_mtime for item in subdir_contents(job_dir,5)]
	else:
		job_type = "old"
		times = [os.stat(item).st_mtime for item in subdir_contents(job_dir,3)]
	
	print job, "is", job_type, "style"
	
	times.sort()
	most_recent_depo = times[-1]
		
	if job in job_atime:
		if job_atime[job] >= most_recent_depo:
			continue
	
	#else:
	#	job_atime[job] = os.stat(job_dir).st_mtime

	results = {} #dict(dict(dict(dict())))#{ { { {}, }, }, }
	position_plots = []
	tilt_plots = []
	rot_plots = []

	#Files to store plots
	job_plots = vis_loc + "/" + job
	try:
		os.mkdir(job_plots)
	except:
		pass #probably already exists
	position_file = job_plots + "/bilayer_position.dat"
	tilt_file = job_plots + "/tilt.dat"
	rot_file = job_plots + "/rot.dat"
	kink_file = job_plots + "/kink.dat"
	internal_tilt_file = job_plots + "/internal_tilt.dat"

	backup_file(position_file)
	backup_file(tilt_file)
	backup_file(rot_file)
	backup_file(kink_file)
	backup_file(internal_tilt_file)
	
	simulation_type = None
	
	for forcefield in os.listdir(job_dir):
		results[forcefield] = {}
		forcefield_dir = job_dir + "/" + forcefield
		for sequence in os.listdir(forcefield_dir):
				results[forcefield][sequence] = {}
				sequence_dir = forcefield_dir + "/" + sequence
				#Need to test for old and new styles
				#print "Test style", sequence_dir + "/" + os.listdir(sequence_dir)[0] + "/" + os.listdir(sequence_dir + "/" + os.listdir(sequence_dir)[0])[0]
				
				if "-" in sequence:
					#Dimer
					simulation_type = "Dimer"
					break
				
				if job_type == "old":
							lipid = "DPPC"
							mix = "1"
							
							results[forcefield][sequence][lipid] = {}
							visgen(forcefield,lipid,sequence,mix,position_plots,rot_plots,tilt_plots,results,sequence_dir)
							continue

				for lipid in os.listdir(sequence_dir):
					results[forcefield][sequence][lipid] = {}
				 	lipid_dir = sequence_dir + "/" + lipid
				 	for mix in os.listdir(lipid_dir):
				 			mix_dir = lipid_dir + "/" + mix
				 			visgen(forcefield,lipid,sequence,mix,position_plots,rot_plots,tilt_plots,results,mix_dir)
	#Now we have the results, make something pretty to look at 
	if simulation_type == "Dimer":
		continue
	
	gallery_listing =[]
	
	with open(job_plots  + "/index.html","w") as summary_output:
		webpage_header(summary_output)
		print >> summary_output, """
<div align="center">
<h1>""" + job + """</h1>
</div>
<div align="left">
<p>
<table class="datasheet" border="3" cellpadding="10" cellspacing="0" width="100%" id="AutoNumber1">
  <tr class="even">

"""

		header = "%-40s%-40s%-5s%-10s%-40s%-100s" % ("<th>Forcefield</th>","<th>Sequence</th>","<th>Lipid</th>","<th>Mix</th>","<th>% Run Insertion</th>","<th>% Time Insertion</th><th>Total Runs</th><th>Energy of Insertion (kCal)</th><th>Data Files</th>")
		print >> summary_output, header
		print >> summary_output, "</tr><tr class=\"odd\">"
		rclass = "even"
		keyring = results.keys()
		keyring.sort()
		for ff in keyring:
			skeyring = results[ff].keys()
			skeyring.sort()
			for seq in skeyring:
					lkeyring = results[ff][seq].keys()
					lkeyring.sort()
					for lipid in lkeyring:
							mkeyring = results[ff][seq][lipid].keys()
							mkeyring.sort()
							for mix in mkeyring:
								
								result = "<td>%-40s</td><td>%-40s</td><td>%-20s</td><td>%-20s</td><td>%-40.3f</td><td>%-40.3f</td><td>%-5d</td>" %(ff,seq,lipid,mix,results[ff][seq][lipid][mix][0],results[ff][seq][lipid][mix][1],results[ff][seq][lipid][mix][3])
								if results[ff][seq][lipid][mix][-1] == "TLD":
									result += "<td>TLD</td>"
								else:
									energy_addition = "<td>%-10.3f%s%-10.3f</td>" %(results[ff][seq][lipid][mix][-2],"&plusmn; ",results[ff][seq][lipid][mix][-1])
									result += energy_addition
								loc_mod = "http://glados.bioch.ox.ac.uk/Results/" + job + "/" + ff + "-" + seq + "-" + lipid + "-" + mix + "-"
								result += "<td><a href="+ loc_mod + "bilayer_position.dat.png"  +">Position</a>,<a href="+ loc_mod + "internal_tilt.dat.png"  +">Tilt</a>,<a href="+ loc_mod + "rot.dat.png"  +">Rotation</a></td>"
								
								gallery_listing.append({"forcefield":ff,"sequence":seq,"lipid":lipid,"mix":mix,"position_image":loc_mod+ "bilayer_position.dat.png", "tilt_image": loc_mod+ "internal_tilt.dat.png", "rotation_image":loc_mod+ "rot.dat.png", "identity":ff + "-" + seq + "-" + lipid + "-" + mix, "details":loc_mod[:-1]+"/base.html"})

								print >> summary_output, result
								print >> summary_output, "</tr><tr class=\"" + rclass +"\">"
								if rclass == "even":
									rclass= "odd"
								else:
									rclass="even"
		gallery_location = pub_loc + job + "/gallery.html"
		status_location = pub_loc + job + "/status.html"
		grid_status_location = pub_loc + "grid_status.html"
		print >> summary_output, '</tr></table></p></div><div align="center"><p><a href="' + status_location + '">Status</a> <a href="' + pub_loc + job + '/' + '">Summary</a> <a href="' + gallery_location + '">Gallery</a> <a href="' + grid_status_location + '">Grid Status</a></p></div></body></html>'
	for item in position_plots:
		AutomatedPlot.hist_com_difference_plot_focus(item,item+".png",graph_title=item.split("-")[-4])
	for item in tilt_plots:
		AutomatedPlot.hist_tilt_data_plot(item,item+".png",graph_title=item.split("-")[-4],binwidth=2)
	for item in rot_plots:
		AutomatedPlot.rotation_plot(item,item+".png",graph_title=item.split("-")[-4])
	
	with open(job_plots  + "/gallery.html","w") as gallery_output:
		webpage_header(gallery_output)
		draw_string = '\'draw_cartoons(['
		#onload_string = '<body onload=\'draw_cartoons(['
		for item in gallery_listing:
			draw_string += '"' + item["identity"] +'",'
		draw_string = draw_string[:-1] + '])\''
		onload_string =  '<body onload='+ draw_string + ">"
		print >> gallery_output, onload_string
		print >> gallery_output, """
<div align="center">
<h1>""" + job + """</h1>

<div align="left">
<p>
<table  border="3" cellpadding="10" cellspacing="0" width="100%" id="AutoNumber1"  class="draggable datasheet">
  <tr>

"""
		print >> gallery_output, "<th>Sequence</th>"
		for item in gallery_listing:
			print >> gallery_output, "<th>"+ item["sequence"] + "</th>"
		print >> gallery_output, "</tr><tr><td>Forcefield</td>"
		for item in gallery_listing:
			print >> gallery_output, "<td>"+ item["forcefield"] + "</td>"
		print >> gallery_output, "</tr><tr><td>Lipid</td>"
		for item in gallery_listing:
			print >> gallery_output, "<td>"+ item["lipid"] + "</td>"		
		print >> gallery_output, "</tr><tr><td>Mix</td>"
		for item in gallery_listing:
			print >> gallery_output, "<td>"+ item["mix"] + "</td>"
		print >> gallery_output, "</tr><tr><td>Position</td>"
		for item in gallery_listing:
			print >> gallery_output, '<td><a href="'+ item["position_image"] + '"><img src="' + item["position_image"] + '" WIDTH="270" ></a></td>'			
		print >> gallery_output, "</tr><tr><td>Tilt</td>"
		for item in gallery_listing:
			print >> gallery_output, '<td><a href="'+ item["tilt_image"] + '"><img src="' + item["tilt_image"] + '" WIDTH="270" ></a></td>'
		print >> gallery_output, "</tr><tr><td>Rotation</td>"
		for item in gallery_listing:
			print >> gallery_output, '<td><a href="'+ item["rotation_image"] + '"><img src="' + item["rotation_image"] + '" WIDTH="270" ></a></td>'
		print >> gallery_output, "</tr><tr><td>Cartoon</td>"
		for item in gallery_listing:
			print >> gallery_output, '<td><canvas width="270" height="270" id="'+item["identity"]+'"></canvas></td>'
			#print item["identity"]
		print >> gallery_output, "</tr><tr><td>Detail</td>"
		for item in gallery_listing:
			print >> gallery_output, '<td><a href="'+ item["details"] + '">Details</a></td>'
		print >> gallery_output, "</tr></table>" + '<div align="center"><p><a href="' + status_location + '">Status</a> <a href="' + pub_loc + job+ '/' + '">Summary</a> <a href="' + gallery_location + '">Gallery</a> <a href="' + grid_status_location + '">Grid Status</a></p></div> ' 
		print >> gallery_output, """<div id="settings">
<h2>Settings</h2>
<div id="dataset">
<p>Select dataset:
<select id="datasetSelect" onchange =""",
		print >> gallery_output, draw_string
		print >> gallery_output, """>
<option value="complete.xml" >Complete</option>
<option value="transmembrane.xml" >Transmembrane</option>
<option value="interfacial.xml" >Interfacial</option>
</select>
</div>

<div id="rotateView">
View front or back of helix?
<select id="rotateViewSelect" onchange ="""
		print >> gallery_output, draw_string
		print >> gallery_output, """>
<option value="false" >Front</option>
<option value="true" >Back</option>
</p>
</select>
</div>

</div>
"""
		print >> gallery_output, "</body></html>"
	
	job_atime[job] = most_recent_depo

	#print "Job", job, "complete at atime:", str(os.stat(job_dir).st_atime)

#Finally, store the last modified times

with open(job_pickle,"w") as pickle_location:
                pickle.dump(job_atime,pickle_location)

