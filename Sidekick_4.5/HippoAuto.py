#!/usr/bin/python

import os
#Ensure that there is a writable folder for writing configuration files
try:
 os.mkdir('/tmp/.matplotlib')
except:
 pass
os.environ['MPLCONFIGDIR']= '/tmp/.matplotlib'

import HAConf
import numpy
from pylab import *

pymol = HAConf.programs['pymol']
scripts = HAConf.configuration['misc_scripts']

def translation_rotation(mutant):
	source_files = HAConf.configuration['hippo_files']
	run_command = pymol + """ -cq """ + scripts + "BuildHelix.py -d """
	#run_command = HAConf.programs['pymol'] + """ -cq ~/Scripts/pymol/BuildHelix.py -d """
	run_command = run_command + """'BuildMe("""
	run_command = run_command + '"' + mutant + '"' + """)'"""
	os.system(run_command)
	
	for item in os.listdir(source_files):
			try:
				os.symlink(source_files+"/"+item,"./"+item)
			except OSError, err:
				if err.errno == errno.EEXIST:
					print "Run file", item, "already exists. Using new version"
					try:
						os.remove("./"+item)
						os.symlink(source_files+"/"+item
,"./"+item)
					except: raise
				else: raise

	
	os.system("./hippo_mac")

	plot_matrix(mutant,"energy_surface.png")
	
def plot_matrix(graph_title,output_png):
	#Now, generate a plot of the energy surface
	input_file = open("e_z_theta_matrix_tm.dat",'r')
	ignore = input_file.readline().split()
	anno_x = [float(ignore[0]),float(ignore[-1])]

	components = []

	anno_y_low = ''

	for line in input_file:
		#print line
		tmp = line.split()
		if anno_y_low == '': anno_y_low = float(tmp[0])
		components.append(tmp[1:])

	anno_y_high = float(tmp[0])

	input_file.close()
	float_components = [ [float(col) for col in row] for row in components]

	X,Y = meshgrid(arange(anno_x[0],anno_x[1]+0.001,20),arange(anno_y_low,anno_y_high+0.001,2))
	Z = numpy.array(float_components)
	figure()
	#im = imshow(Z, interpolation='bicubic', cmap=cm.gray,extent=(0,180,-50,50),vmin=-40,vmax=120)
	grid = contourf(X,Y,Z,10,cmap=cm.hot,norm=normalize(-40,120))
	xlabel("Angle wrt bilayer (degrees)")
	ylabel("Depth in membrane (A)")
	title(graph_title)
	#hot()
	CB = colorbar(grid, shrink=0.8, extend='both',spacing='proportional')
	#CBI = colorbar(im, orientation='horizontal', shrink=0.8, ticks=range(-40,120,20))
	savefig(output_png)
	
if __name__ == "__main__": translation_rotation(sys.argv[1])
