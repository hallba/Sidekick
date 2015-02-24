#!/usr/bin/python

#HELANAL algorithm. Rewritten for clarity

from CartesianToolkit import *
#import math, numpy, sys
import math,sys,os
import xtcio,ndxio

try:
	import psyco
except:
	pass

def helanal_trajectory(xtcfile,pdbfile,ndxfile,start,end,begin,finish,matrix_filename,origin_pdbfile,summary_filename,verbose):
	pdbdata = open(pdbfile)
	[indices,pdblength] = get_backbone_indices(pdbdata)
	pdbdata.close()
	trajectory  = xtcio.read_xtc(xtcfile)
	#origin_pdbfile = "origin.pdb"
	#matrix_filename = "bending_matrix.txt"
	backup_file(matrix_filename)
	backup_file(origin_pdbfile)
	backup_file(summary_filename)
	global_height = []
	global_twist = []
	global_rnou = []
	global_bending = []
	global_bending_matrix = []
	global_tilt = []		
	
	if start != None and end != None:
			print "Analysing from residue", start, "to", end
	elif start != None and end == None:
			print "Analysing from residue", start, "to the C termini"
	elif start == None and end != None:
			print "Analysing from the N termini to", end
	while(trajectory.next_frame()):
		if begin != None:
			if trajectory.time < begin:
				continue
		if finish != None:
			if trajectory.time > finish:
				break
		if pdblength != trajectory.natoms:
			print "Different numbers of atoms in the pdb and xtc files (", pdblength , "vs" , trajectory.natoms ,")! Exiting"
			exit()
		trajectory.convert_coordinates()
		ca_positions = index_cartesians(indices,trajectory.cartesian)
		max_length = len(ca_positions)
		if start != None and end != None:
			#print "Analysing from residue", start, "to", end
			ca_positions = ca_positions[(start-1):(end)]
		elif start != None and end == None:
			#print "Analysing from residue", start, "to the C termini"
			ca_positions = ca_positions[(start-1):]
		elif start == None and end != None:
			#print "Analysing from the N termini to", end
			ca_positions = ca_positions[:(end)]
		if len(global_height) == 0:
			measured_length = len(ca_positions)
			print "Analysing", measured_length, "/", max_length, "residues"
		
		[twist,bending_angles,height,rnou,origins,local_helix_axes] = main_loop(ca_positions)
		
		origin_pdb(origins,origin_pdbfile)
		
		#calculate local bending matrix( it is looking at all i, j combinations)
		if len(global_bending_matrix) == 0:
			global_bending_matrix = [ [ [] for item in local_helix_axes] for item in local_helix_axes ]
		for i,row in zip(local_helix_axes,global_bending_matrix):
			for j,col in zip(local_helix_axes,row):
				if i == j:
					angle = 0.
				else:
					angle = math.acos(vecscaler(i,j))*180/math.pi
				col.append(angle)
		
		global_height += height
		global_twist  += twist
		global_rnou   += rnou
		if len(global_bending) == 0:
			global_bending = [ [] for item in bending_angles ]
			global_tilt = [ [] for item in local_helix_axes ]
		for store,tmp in zip(global_bending,bending_angles): store.append(tmp)
		for store,tmp in zip(global_tilt,local_helix_axes): store.append(vecangle(tmp,[0,0,1]))
		#simple ticker
		formated_time = "%20.1f" % trajectory.time
		print '\r',formated_time,' ps',
		sys.stdout.flush()

	print '\nComplete'
	[twist_mean, twist_sd, twist_abdev] = stats(global_twist)
	[height_mean, height_sd, height_abdev] = stats(global_height)
	[rnou_mean, rnou_sd, rnou_abdev] = stats(global_rnou)
	
	bending_statistics = [ stats(item) for item in global_bending]
	tilt_statistics =    [ stats(item) for item in global_tilt]

	bending_statistics_matrix = [[stats(col) for col in row] for row in global_bending_matrix]
	mat_output = open(matrix_filename,'w')
	print >> mat_output, "Mean"
	#[ [ print >> mat_output, "%8.3f\t", % col[0] for col in row ] and print '' for row in bending_statistics_matrix ]
	for row in bending_statistics_matrix:
		for col in row:
			formatted_angle = "%6.1f" % col[0]
			print >> mat_output, formatted_angle,
		print >> mat_output, ''
	print >> mat_output, "\nSD"
	#[ [ print >> mat_output, "%8.3f\t", % col[0] for col in row ] and print '' for row in bending_statistics_matrix ]
	for row in bending_statistics_matrix:
		for col in row:
			formatted_angle = "%6.1f" % col[1]
			print >> mat_output, formatted_angle,
		print >> mat_output, ''
	print >> mat_output, "\nABDEV"
	#[ [ print >> mat_output, "%8.3f\t", % col[0] for col in row ] and print '' for row in bending_statistics_matrix ]
	for row in bending_statistics_matrix:
		for col in row:
			formatted_angle = "%6.1f" % col[2]
			print >> mat_output, formatted_angle,
		print >> mat_output, ''
	
	mat_output.close()
	print "Height:", height_mean, "SD", height_sd, "ABDEV", height_abdev, '(nm)'
	print "Twist:", twist_mean, "SD", twist_sd, "ABDEV", twist_abdev
	print "Residues/turn:", rnou_mean, "SD", rnou_sd, "ABDEV", rnou_abdev
	print "Local bending angles:"
	residue_statistics = zip(*bending_statistics)
	measure_names = ["Mean ","SD   ","ABDEV"]
	print "ResID",
	if start == None:
		for item in range(4,len(residue_statistics[0])+4):
			output = "%8d" % item
			print output,
	else:
		for item in range(start+3,len(residue_statistics[0])+start+3):
			output = "%8d" % item
			print output,
	print ""
	for measure,name in zip(residue_statistics,measure_names):
		print name,
		for residue in measure:
			output = "%8.1f" % residue
			print output,
		print ''
	print "Local tilts:"
	residue_statistics = zip(*tilt_statistics)
	measure_names = ["Mean ","SD   ","ABDEV"]
	
	print "ResID",
	if start == None:
		for item in range(1,len(residue_statistics[0])+1):
			#output = "%8d" % item
			output = "%3d->%3d" % (item, item+3)
			print output,
	else:
		for item in range(start,len(residue_statistics[0])+start):
			#output = "%8d" % item
			output = "%3d->%3d" % (item, item+3)
			print output,
	print ""
	
	for measure,name in zip(residue_statistics,measure_names):
		print name,
		for residue in measure:
			output = "%8.1f" % (tilt_correct(180*residue/math.pi))
			print output,
		print ''

	#Verbose options- prints out a tilt for each 
	if verbose:
		if start == None:
			for item,tilts in zip(range(1,len(residue_statistics[0])+1),global_tilt):
				#output = "%8d" % item
				output = summary_filename + "%3d->%3d.tilt" % (item, item+3)
				output_file = open(output,"w")
				verbose_count= 0
				for tilt_value in tilts:
					print >> output_file, verbose_count, (tilt_correct(180*tilt_value/math.pi))
					verbose_count +=1
				output_file.close()				
		else:
			for item,tilts in zip(range(start,len(residue_statistics[0])+start),global_tilt):
				#output = "%8d" % item
				output = summary_filename + "%3d->%3d.tilt" % (item, item+3)
				output_file = open(output,"w")
				verbose_count = 0
				for tilt_value in tilts:
					print >> output_file, verbose_count, (tilt_correct(180*tilt_value/math.pi))
					verbose_count +=1
				output_file.close()
	


	summary_output = open(summary_filename,'w')
	print >> summary_output, "Height:", height_mean, "SD", height_sd, "ABDEV", height_abdev, '(nm)'
	print >> summary_output, "Twist:", twist_mean, "SD", twist_sd, "ABDEV", twist_abdev
	print >> summary_output, "Residues/turn:", rnou_mean, "SD", rnou_sd, "ABDEV", rnou_abdev
	print >> summary_output, "Local bending angles:"
	residue_statistics = zip(*bending_statistics)
	measure_names = ["Mean ","SD   ","ABDEV"]
	print >> summary_output,  "ResID",
	if start == None:
		for item in range(4,len(residue_statistics[0])+4):
			output = "%8d" % item
			print >> summary_output,  output,
	else:
		for item in range(start+3,len(residue_statistics[0])+start+3):
			output = "%8d" % item
			print >> summary_output, output,
	print >> summary_output, ""

	for measure,name in zip(residue_statistics,measure_names):
		print >> summary_output, name,
		for residue in measure:
			output = "%8.1f" % residue
			print >> summary_output, output,
		print >> summary_output, ''
	print  >> summary_output, "Local tilts:"
	residue_statistics = zip(*tilt_statistics)
	measure_names = ["Mean ","SD   ","ABDEV"]
	
	print >> summary_output, "ResID",
	if start == None:
		for item in range(1,len(residue_statistics[0])+1):
			#output = "%8d" % item
			output = "%3d->%3d" % (item, item+3)
			print >> summary_output, output,
	else:
		for item in range(start,len(residue_statistics[0])+start):
			#output = "%8d" % item
			output = "%3d->%3d" % (item, item+3)
			print >> summary_output, output,
	print >> summary_output, ""
	
	for measure,name in zip(residue_statistics,measure_names):
		print >> summary_output, name,
		for residue in measure:
			output = "%8.1f" % (tilt_correct(180*residue/math.pi))
			print >> summary_output, output,
		print >> summary_output, ''

	summary_output.close()

def tilt_correct(number):
	'''Changes an angle (in degrees) so that it is between 0 and 90'''
	if number < 90:
		return number
	else:
		return 180 - number

def backup_file(filename):
	if os.path.exists(filename):
		target_name = "#" + filename
		failure = True
		if not os.path.exists(target_name):
			os.rename(filename,target_name)
			failure = False
		else:
			for i in range(20):
				alt_target_name = target_name + "." + str(i)
				if os.path.exists(alt_target_name):
					continue
				else:
					os.rename(filename,alt_target_name)
					failure = False
					break
		if failure:
			print "Too many backups. Clean up and try again"
			exit()

def stats(some_list):
	list_mean = mean(some_list)
	list_sd = sample_sd(some_list,list_mean)
	list_abdev = mean_abs_dev(some_list,list_mean)
	return [list_mean, list_sd, list_abdev]

def helanal_main(pdbfile,start,end):
	pdbdata = open(pdbfile)
	positions = get_backbone_positions(pdbdata)
	max_length = len(positions)

	if start != None and end != None:
			print "Analysing from residue", start, "to", end
			positions = positions[(start-1):(end)]
	elif start != None and end == None:
			print "Analysing from residue", start, "to the C termini"
			positions = positions[(start-1):]
	elif start == None and end != None:
			print "Analysing from the N termini to", end
			positions = positions[:(end)]

	measured_length = len(positions)
	print "Analysing", measured_length, "/", max_length, "residues"
	[twist,bending_angles,height,rnou,origins,local_helix_axes] = main_loop(positions)

	#TESTED- origins are correct
	#print current_origin
	#print origins
	
	max_angle = max(bending_angles)
	mean_angle = mean(bending_angles)
	#sd calculated using n-1 to replicate original fortran- assumes a limited sample so uses the sample standard deviation
	sd_angle = sample_sd(bending_angles,mean_angle)
	mean_absolute_deviation_angle = mean_abs_dev(bending_angles,mean_angle)
	#TESTED- stats correct
	#print max_angle, mean_angle, sd_angle, mean_absolute_deviation_angle
	
	#calculate local bending matrix(now it is looking at all i, j combinations)
	for i in local_helix_axes:
		for j in local_helix_axes:
			if i == j:
				angle = 0.
			else:
				angle = math.acos(vecscaler(i,j))*180/math.pi
			string_angle = "%6.0f\t" % angle
			#print string_angle,
		#print ''
		#TESTED- local bending matrix!
	
	#Average helical parameters
	mean_twist = mean(twist)
	sd_twist = sample_sd(twist,mean_twist)
	abdev_twist = mean_abs_dev(twist,mean_twist)
	#TESTED-average twists
	#print mean_twist, sd_twist, abdev_twist
	mean_rnou = mean(rnou)
	sd_rnou = sample_sd(rnou,mean_rnou)
	abdev_rnou = mean_abs_dev(rnou,mean_rnou)
	#TESTED-average residues per turn
	#print mean_rnou, sd_rnou, abdev_rnou
	mean_height = mean(height)
	sd_height = sample_sd(height,mean_height)
	abdev_height = mean_abs_dev(height,mean_height)	
	#TESTED- average rises
	
	print "Height:", mean_height, sd_height, abdev_height
	print "Twist:", mean_twist, sd_twist, abdev_twist
	print "Residues/turn:", mean_rnou, sd_rnou, abdev_rnou
	#print mean_height, sd_height, abdev_height
	print "Local bending angles:"
	for angle in bending_angles:
		output = "%8.1f\t" % angle
		print output,
	print ''
	
def origin_pdb(origins,pdbfile):	#special- print origins to pdb, assumes the need to convert nm->A
	output = open(pdbfile,'a')
	i=1
	for item in origins:
		tmp = "ATOM    %3d  CA  ALA   %3d    %8.3f%8.3f%8.3f  1.00  0.00" % (i,i,item[0]*10,item[1]*10,item[2]*10)
		print >> output, tmp
		i += 1
	print >> output, "TER\nENDMDL"
	output.close()

def main_loop(positions):
	twist = []
	rnou = []
	height = []
	origins = [[0.,0.,0.] for item in positions[:-2]]
	local_helix_axes = []
	for i in range(len(positions)-3):
		vec12 = vecsub(positions[i+1],positions[i])
		vec23 = vecsub(positions[i+2],positions[i+1])
		vec34 = vecsub(positions[i+3],positions[i+2])
		
		dv13 = vecsub(vec12,vec23)
		dv24 = vecsub(vec23,vec34)
		#direction of the local helix axis
		current_uloc = vecnorm(veccross(dv13,dv24))
		local_helix_axes.append(current_uloc)
		
		#TESTED- Axes correct
		#print current_uloc

		dmag = veclength(dv13)
		emag = veclength(dv24)
		
		costheta = vecscaler(dv13,dv24)/(dmag*emag)
		#rnou is the number of residues per turn
		current_twist = math.acos(costheta)
		twist.append(current_twist*180/math.pi)
		rnou.append(2*math.pi/current_twist)
		#radius of local helix cylinder radmag
		
		costheta1 = 1.0 - costheta
		radmag = (dmag*emag)**0.5/(2*costheta1)
		
		#Height of local helix cylinder
		current_height = vecscaler(vec23,current_uloc)
		height.append(current_height)
		#TESTED- Twists etc correct
		#print current_twist*180/math.pi, 2*math.pi/current_twist, height
		
		dv13 = vecnorm(dv13)
		dv24 = vecnorm(dv24)
		
		rad = [radmag * item for item in dv13]
		current_origin = [(item[0] - item[1]) for item in zip(positions[i+1],rad)]
		origins[i] = current_origin
		
		#TESTED- origins are correct
		#print current_origin
		
		rad = [radmag * item for item in dv24]
		current_origin = [(item[0] - item[1]) for item in zip(positions[i+2],rad)]
		origins[i+1] = current_origin
	#local bending angles (eg i > i+3, i+3 > i+6)
	bending_angles = [0 for item in range(len(local_helix_axes)-3)]
	for axis in range(len(local_helix_axes)-3):
		angle = math.acos(vecscaler(local_helix_axes[axis],local_helix_axes[axis+3]))*180/math.pi
		bending_angles[axis] = angle
		#TESTED- angles are correct
		#print angle

	return [twist,bending_angles,height,rnou,origins,local_helix_axes]

def fit(origins):
	#Subroutine to fir plane, circle and line to local helix origins
	#INCOMPLETE
	
	#Not sure exactly what these represent
	[x,y,z] = [0., 0., 0.]
	[x2,y2,z2] = [0., 0., 0.]
	[xy,xz,yz] = [0., 0., 0.]
	
	matp = [ [ 0. for i in range(3)] for i in range(3) ] 
	for item in origins:
		x2 += item[0]**2
		y2 += item[1]**2
		z2 += item[2]**2
		xy += item[0]*item[1]
		xz += item[0]*item[2]
		yz += item[1]*item[2]
		x  += item[0]
		y  += item[1]
		z  += item[2]
	matp[0][0] += x2
	matp[0][1] += xy
	matp[1][0] += xy
	matp[0][2] += xz
	matp[2][0] += xz
	matp[1][1] += y2
	matp[2][2] += z2
	matp[1][2] += yz
	matp[2][1] += yz
	
	matp = numpy.matrix(matp)
	pmat = matp.I
	

def index_cartesians(index_list,frame):
	return [frame[item] for item in index_list]	

def mean(some_list):
	return sum(some_list)/len(some_list)
def sample_sd(some_list,mean):
	return (sum([(item - mean)**2 for item in some_list])/(len(some_list)-1))**0.5
def mean_abs_dev(some_list,mean):
	return sum([math.fabs(item - mean) for item in some_list])/len(some_list)
		
def get_backbone_positions(pdbdata):
	positions = []
	for line in pdbdata:
		atom_name = line[11:16].strip()
		if atom_name =="CA" or (atom_name[:1] == "B" or atom_name[:2] == "0B"): #Bondini, atomistic or MARTINI
			positions.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
	return positions

def get_backbone_indices(pdbdata):
	indices = []
	current_index = 0
	pdblength = 0
	for line in pdbdata:
		if line[:4] != 'ATOM':
			continue
		pdblength += 1
		atom_name = line[11:16].strip()
		if atom_name =="CA" or (atom_name[:1] == "B" or atom_name[:2] == "0B"  or atom_name[:2] == "5B" ): #Bondini, atomistic or MARTINI
			indices.append(current_index)
		current_index +=1
	print "Length", len(indices)
	return indices, pdblength

def getoptions():
	from optparse import OptionParser
	parser = OptionParser()
	
	parser.add_option("-p", "--pdb", dest="pdb_file", default=None,
	         help="PDB file", metavar="FILE")
	parser.add_option("-x", "--xtc", dest="xtc_file", default=None,
	         help="XTC file", metavar="FILE")
	parser.add_option("-n", "--ndx", dest="ndx_file", default=None,
	         help="NDX file", metavar="FILE")
	parser.add_option("-m", "--matrix", dest="mat_file", default="bending_matrix.txt",
	         help="Output file- bending matrix", metavar="FILE")
	parser.add_option("-o", "--origin", dest="ori_file", default="origin.pdb",
	         help="Output file- origin pdb file", metavar="FILE")
	parser.add_option("-r", "--roundup", dest="summary_file", default="summary.txt",
	         help="Output file- all of the basic data", metavar="FILE")
	parser.add_option("-s", "--start",
	         dest="start", default=None,type="int",
	         help="Start residue")
	parser.add_option("-e", "--end",
	         dest="end", default=None,type="int",
	         help="End residue")
	parser.add_option("-b", "--begin",
	         dest="begin", default=None,type="int",
	         help="Begin analysing from time (ps)")
	parser.add_option("-f", "--finish",
	         dest="finish", default=None,type="int",
	         help="Stop analysing after time (ps)")
	parser.add_option("-v","--verbose",
			 dest="verbose", default=False, action="store_true",
			 help="Generate a file for each residue's information")

	(options, args) = parser.parse_args()
	return options
	
	

if __name__ == "__main__":
	import sys
	options = getoptions()
	if options.pdb_file == None:
		print "Needs a PDB file"
		exit()
	if options.xtc_file == None:
		print "No xtc file- working on the pdb alone"
		helanal_main(options.pdb_file,options.start,options.end)
		exit()
	helanal_trajectory(options.xtc_file,options.pdb_file,options.ndx_file,options.start,options.end,options.begin,options.finish,options.mat_file,options.ori_file,options.summary_file,options.verbose)
