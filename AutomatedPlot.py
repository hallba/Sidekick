#!/usr/bin/python

from pylab import *
import numpy
import matplotlib.cm as cm
import matplotlib.colors as colors

def tilt_data_plot(filename,image_output,graph_title=''):
	time = []
	tilt = []
	for line in open(filename):
		if line[0] == "#" or line[0] == "@": continue
		line_content = line.split()
		time.append(line_content[0])
		if float(line_content[1]) <= 90: tilt.append(float(line_content[1]))
		else: tilt.append(180 - float(line_content[1]))
	fig = figure()
	plot(time,tilt)
	xlabel('Time (ps)')
	ylabel('Tilt angle (degrees)')
	title(graph_title)
	savefig(image_output)
	#fig.close()
	
def hist_tilt_data_plot(filename,image_output,graph_title=''):
	tilt = []
	for line in open(filename):
		if line[0] == "#" or line[0] == "@": continue
		line_content = line.split()
		if float(line_content[1]) <= 90: tilt.append(float(line_content[1]))
		else: tilt.append(180 - float(line_content[1]))
	fig = figure()
	N,bins,patches = hist(tilt,range(0,90,5))
	fracs = N.astype(float)/N.max()
	norm = colors.normalize(fracs.min(), fracs.max())

	for thisfrac, thispatch in zip(fracs, patches):
		color = cm.jet(norm(thisfrac))
		thispatch.set_facecolor(color)
	xlabel('Tilt angle (degrees)')
	ylabel('Frequency')
	title(graph_title)
	savefig(image_output)
	#fig.close()

def hist_com_difference_plot(filename,image_output,graph_title=''):
	distance = []
	for line in open(filename):
                if line[0] == "#" or line[0] == "@": continue
                line_content = line.split()
		distance.append(float(line_content[1]))
	fig = figure()
	N,bins,patches = hist(distance,range(-30,30))
        fracs = N.astype(float)/N.max()
        norm = colors.normalize(fracs.min(), fracs.max())

	for thisfrac, thispatch in zip(fracs, patches):
        	color = cm.jet(norm(thisfrac))
        	thispatch.set_facecolor(color)

	xlabel('Z Distance between Protein and Bilayer centers')
	ylabel('Frequency')
	title(graph_title)
	savefig(image_output)

def bilayer_stats(filename):
	boxes = []
	line_counter = 0
	for line in open(filename):
		line_content =  line.split()
		if len(boxes) == 0: 
			boxes = [float(item) for item in line_content]
		else:
			for target in range(len(line_content)):
				boxes[target] += float(line_content[target])
		line_counter += 1
	#average
	averages = [item/line_counter for item in boxes]
	#st dev
	boxes = [0.0 for item in boxes]
	for line in open(filename):
			line_content = line.split()
			data = [float(item) for item in line_content]
			for i in range(len(data)):
				boxes[i] += (data[i] - averages[i])**2
	sd = [(item/(line_counter-1))**0.5 for item in boxes]
	error = [2*value for value in sd]
	return [averages, sd, error]

def bilayer_deformation(upper_filename,lower_filename,image_output,graph_title=''):
	[upper_ave, upper_sd, upper_error] = bilayer_stats(upper_filename)
	[lower_ave, lower_sd, lower_error] = bilayer_stats(lower_filename)

	x = [left_edge+2.5 for left_edge in range(0,len(upper_sd)*5,5)]
	fig = figure()
	#plot(x,averages)
	#print "Averages\n", averages, "\nSD\n", sd
	#error is 95.4% confidence limits based on 2*sd
	#ylabel('Distance of phophate groups from center of bilayer')
	subplot(2,1,1)

	errorbar(x,upper_ave,yerr=upper_error)
	v = axis()
	axis([0,len(upper_sd)*5,v[2],v[3]])
	#xlabel('Distance from protein center of mass')
	#ylabel('Distance of phophate groups from center of bilayer')
	ylabel('Leaflet thickness')
	title(graph_title+"- Upper Leaflet")
	
	subplot(2,1,2)

	errorbar(x,lower_ave,yerr=lower_error)
	v = axis()
	axis([0,len(lower_sd)*5,v[2],v[3]])
	xlabel('Distance from protein center of mass')
	#ylabel('Distance of phophate groups from center of bilayer')
	ylabel('Leaflet thickness')
	title(graph_title+"- Lower Leaflet")
	savefig(image_output)

def bilayer_deformation_bully(upper_filename,lower_filename,image_output,graph_title=''):
	[upper_ave, upper_sd, upper_error] = bilayer_stats(upper_filename)
	[lower_ave, lower_sd, lower_error] = bilayer_stats(lower_filename)

	#x = [left_edge+2.5 for left_edge in range(0,len(upper_sd)*5,5)]
	
	theta = linspace(0.,2.*pi,101)
	R = linspace(0.,2*pi,len(upper_ave))
	Y,X = meshgrid(R, theta)
	
	fig = figure()
	#plot(x,averages)
	#print "Averages\n", averages, "\nSD\n", sd
	#error is 95.4% confidence limits based on 2*sd
	#ylabel('Distance of phophate groups from center of bilayer')
	
	top = subplot(2,1,1,polar=True,frameon=True)
	
	upperz = numpy.array([upper_ave]*101)
	BullyT = top.pcolormesh(X,Y,upperz,cmap=cm.hot,norm=normalize(15,25))	# make sure aspect ratio preserved	top.set_aspect('equal')
	
	CBt = colorbar(BullyT, shrink=0.8, extend='both',spacing='proportional')

	title(graph_title+"- Upper Leaflet")
	
	bottom = subplot(2,1,2,polar=True,frameon=True)
	
	lowerz = numpy.array([lower_ave]*101)
	BullyB = bottom.pcolormesh(X,Y,lowerz,cmap=cm.hot,norm=normalize(15,25))	# make sure aspect ratio preserved	bottom.set_aspect('equal')
	
	CBb = colorbar(BullyB, shrink=0.8, extend='both',spacing='proportional')

	title(graph_title+"- Lower Leaflet")
	savefig(image_output)
	
def bilayer_deformation_combo(upper_filename,lower_filename,image_output,graph_title=''):
	[upper_ave, upper_sd, upper_error] = bilayer_stats(upper_filename)
	[lower_ave, lower_sd, lower_error] = bilayer_stats(lower_filename)

	x = [left_edge+2.5 for left_edge in range(0,len(upper_sd)*5,5)]
		
	fig = figure()
	
	fig.text(.5, .95, graph_title, horizontalalignment='center') 
	#plot(x,averages)
	#print "Averages\n", averages, "\nSD\n", sd
	#error is 95.4% confidence limits based on 2*sd
	#ylabel('Distance of phophate groups from center of bilayer')
	
	subplot(2,2,1)

	errorbar(x,upper_ave,yerr=upper_error)
	v = axis()
	axis([0,len(upper_sd)*5,v[2],v[3]])
	#xlabel('Distance from protein center of mass')
	#ylabel('Distance of phophate groups from center of bilayer')
	ylabel('Leaflet thickness')
	title("Upper Leaflet")
	
	subplot(2,2,3)

	errorbar(x,lower_ave,yerr=lower_error)
	v = axis()
	axis([0,len(lower_sd)*5,v[2],v[3]])
	xlabel('Distance from protein center of mass')
	#ylabel('Distance of phophate groups from center of bilayer')
	ylabel('Leaflet thickness')
	title("Lower Leaflet")
	
	#Polar plots
	
	#This is an ugly hack to prevent the last value from being ignored
	upper_ave.append(0)
	lower_ave.append(0)
	
	theta = linspace(0.,2.*pi,101)
	R = linspace(0.,2*pi,len(upper_ave))
	Y,X = meshgrid(R, theta)

	
	top = subplot(2,2,2,polar=True,frameon=False)
	
	upperz = numpy.array([upper_ave]*101)
	BullyT = top.pcolormesh(X,Y,upperz,cmap=cm.hot,norm=normalize(15,25))
	thetagrids([])
	#print upperz[0]	# make sure aspect ratio preserved	top.set_aspect('equal')
	
	CBt = colorbar(BullyT, shrink=0.8, extend='both',spacing='proportional')

	title("Upper Leaflet")
	
	bottom = subplot(2,2,4,polar=True,frameon=False)
	
	lowerz = numpy.array([lower_ave]*101)
	BullyB = bottom.pcolormesh(X,Y,lowerz,cmap=cm.hot,norm=normalize(15,25))
	thetagrids([])	# make sure aspect ratio preserved	bottom.set_aspect('equal')
	
	CBb = colorbar(BullyB, shrink=0.8, extend='both',spacing='proportional')

	title("Lower Leaflet")
	savefig(image_output)
	
def rotation_plot(filename,image_output,graph_title=''):
	angles = [float(line.split()[1]) for line in open(filename)]
	#ok, a bit of manual trickery for the radial plot histogram
	number_of_bins = 24
	bins = [0 for i in range(number_of_bins)]
	for frame in angles:
		#frame += 180
		target_bin = int(frame*number_of_bins/360)
		if frame < 0: target_bin -= 1
		bins[target_bin] += 1
	#bins = bins[18:] + bins[:18]
	#bins = range(36)
	print bins
	fig = figure()
	
	theta = linspace(0.,2.*pi,(number_of_bins+1))
	R = linspace(0.,2*pi,2)
	Y,X = meshgrid(R, theta)
	Z_data = [[item] * 2 for item in bins]
	Z = numpy.array(Z_data)
	top = subplot(1,1,1,polar=True,frameon=True)
	Radiant = top.pcolormesh(X,Y,Z,cmap=cm.jet)
	thetagrids(range(0,360,10))
	top.set_aspect('equal')
	CBb = colorbar(Radiant, shrink=0.8, extend='both',spacing='proportional')
	title(graph_title)
	savefig(image_output)