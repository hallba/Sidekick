#!/usr/bin/python

import sys,os,pdbio,numpy,math

from HAConf import configuration,programs,debug_status
from LipidBox import add_solvent_dppc
from MDPGenerator import replace_seed, replace_steps
import AutomatedPlot
import xtcio,ndxio,CartesianToolkit
from GenerateCGSystem import allowed_lipids

gromacs = configuration['gromacs_location']
#pymol = programs['pymol']
#vmd = programs['vmd']
scripts = configuration['misc_scripts']

lipid_atoms = 		{	"POPC"	:	13,
						"POPE"	:	13,
						"POPG"	:	13,
						"DOPC"	:	14,
						"DOPE"	:	14,
						"DOPG"	:	14,
						"DLPC"	:	10,
						"DLPE"	:	10,
						"DLPG"	:	10,
						"DPPC"	:	12,
						"DPPE"	:	12,
						"DPPG"	:	12	}
phosphate_atoms = 	{	"POPC"	:	1,
						"POPE"	:	1,
						"POPG"	:	1,
						"DOPC"	:	1,
						"DOPE"	:	1,
						"DOPG"	:	1,
						"DLPC"	:	1,
						"DLPE"	:	1,
						"DLPG"	:	1,
						"DPPC"	:	1,
						"DPPE"	:	1,
						"DPPG"	:	1	}

def index_cartesians(index_list,frame):
	return [frame[item] for item in index_list]

def lipid_list_breakdown(lipid_positions,number_of_atoms):
	return [lipid_positions[i*number_of_atoms:(i+1)*number_of_atoms] for i in range(len(lipid_positions)/number_of_atoms)]

def lipid_list_breakdown_mix(lipid_map,lipid_positions):
	lipid_list = []
	viewpoint = 0
	for item in lipid_map:
		atoms_per_lipid = item[1]
		lipid_subtotal = item[2]
		tmp_positions = lipid_positions[viewpoint:(viewpoint+atoms_per_lipid*lipid_subtotal)]
		viewpoint += atoms_per_lipid*lipid_subtotal
		for i in range(len(tmp_positions)/atoms_per_lipid):
			lipid_list.append(tmp_positions[i*atoms_per_lipid:(i+1)*atoms_per_lipid])
	return lipid_list

def bilayer_center(individual_lipid_positions,lipid_com):
	individual_lipid_comz = [CartesianToolkit.center(lipid)[2] for lipid in individual_lipid_positions]
	upper_list = []
	lower_list = []
	upper_lipids = []
	lower_lipids = []
	for i in range(len(individual_lipid_comz)):
		if individual_lipid_comz[i] > lipid_com[2]:
			upper_list += individual_lipid_positions[i]
			upper_lipids.append(individual_lipid_positions[i])
		else:
			lower_list += individual_lipid_positions[i]
			lower_lipids.append(individual_lipid_positions[i])
	
	upper_list_comz = CartesianToolkit.center(upper_list)[2]
	lower_list_comz = CartesianToolkit.center(lower_list)[2]
	return (((upper_list_comz + lower_list_comz)/2.0),upper_lipids,lower_lipids)

def lipid_list_from_topology(filename):
	lipid_map = []
	input_fh = open(filename)
	for line in input_fh:
		contents = line.split()
		if len(contents) == 0: continue
		#print contents
		#print contents[0]
		if contents[0] in lipid_atoms.keys():
			#print "success", contents
			lipid_map.append([contents[0],lipid_atoms[contents[0]],int(contents[1])])
	input_fh.close()
	return lipid_map

def simulation_statistics(bilayer_normal=2,trajfile="xy_fit.xtc",ndxfile="system.ndx",topfile="ioned_topol.top"):
	trajectory = xtcio.read_xtc(trajfile)
	groups = ndxio.read_ndx(ndxfile)
	bilayer_position_file = open("bilayer_position.dat","w")
	upper_leaflet_file = open("up.dat","w")
	lower_leaflet_file = open("down.dat","w")
	rotation_file = open("rot.dat","w")
	tng_rotation_file = open("comparison.dat","w")
	tilt_file = open("tilt.dat","w")
	
	frame = 0
	helix_insertion = 0.
	lipid_map = lipid_list_from_topology(topfile)
	print "Lipid map" 
	print lipid_map
	
	while(trajectory.next_frame()):
		trajectory.convert_coordinates()
		
		#Do I need to rotate the cartesians to account for non-xy bilayers.
		if bilayer_normal==2:
			pass
		elif bilayer_normal == 1:
			#rotate by 90' in x- rotate y to x, then z to x. There must be a better way
			bilayer_rotmat_first = CartesianToolkit.rotation_matrix_to_x([0,1,0])
			bilayer_rotmat_second = CartesianToolkit.rotation_matrix_to_x([0,0,1])
			#for cart in trajectory.cartesian:
			#	cart = CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat_first)
			#	cart = CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat_second)
			trajectory.cartesian = [CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat_first) for cart in trajectory.cartesian]
			trajectory.cartesian = [CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat_second) for cart in trajectory.cartesian]
		elif bilayer_normal == 0:
			#rotate by 90' in y: x becomes z, z becomes -x
			bilayer_rotmat = CartesianToolkit.rotation_matrix_to_x([0,0,1])
			#for cart in trajectory.cartesian:
			#	cart = CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat)
			trajectory.cartesian = [CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat) for cart in trajectory.cartesian]
		else:
			print "Bilayer normal not feasible. Aborting"
			sys.exit()
		#Calculate atom positions of interest
		try:
			backbone_positions = index_cartesians(groups["B*"],trajectory.cartesian)
		except:
			backbone_positions = index_cartesians(groups["B*_CA"],trajectory.cartesian)
		lipid_positions = []
		for lipid in allowed_lipids:
			try:
				lipid_positions += index_cartesians(groups[lipid],trajectory.cartesian)
				#Need an additional list of lipid groups to calculate the lipid list
			except:
				pass
		'''Meta Code
		lipid_numbers = [ ["DPPC",200], ["DPPG",100] ]
		lipid_numbers = report_lipids_from_topology("ioned_topol.top"<-needs a default?)
		lipid_numbers = report_lipids_from_ndx()
		lipid_numbers = report_lipids_from_gro()
		'''
		
		water_positions = index_cartesians(groups["W"],trajectory.cartesian)
		top_comz = backbone_positions[0][2]
		
		#Calculate centers of interest
		backbone_com = CartesianToolkit.center(backbone_positions)
		lipid_com = CartesianToolkit.center(lipid_positions)
		#(bilayer_midpoint,upper_leaflet,lower_leaflet) = bilayer_center(lipid_list_breakdown(lipid_positions,lipid_atoms["DPPC"]),lipid_com)
		(bilayer_midpoint,upper_leaflet,lower_leaflet) = bilayer_center(lipid_list_breakdown_mix(lipid_map,lipid_positions),lipid_com)
		if top_comz < bilayer_midpoint:
			[upper_leaflet,lower_leaflet] = [lower_leaflet,upper_leaflet]
		
		
		#Calculate bilayer distance
		bkbone_bilayer_distance = backbone_com[2] - bilayer_midpoint
		if top_comz < backbone_com[2]:
			bkbone_bilayer_distance = 0 - bkbone_bilayer_distance
		print >> bilayer_position_file, frame, 10*bkbone_bilayer_distance
		
		
		#Report phosphate positions relative to the bilayer center
		if frame == 0:
			#This is the radius divided by 5 Angstroms. Or, for nm (from xtc) 2*0.5
			if bilayer_normal == 2:
				upper_leaflet_bins = [ [] for i in range(int(trajectory.box[0])) ]
				lower_leaflet_bins = [ [] for i in range(int(trajectory.box[0])) ]
			else:
				upper_leaflet_bins = [ [] for i in range(int(trajectory.box[8])) ]
				lower_leaflet_bins = [ [] for i in range(int(trajectory.box[8])) ]
		for lipid in upper_leaflet:
			phosphate = lipid[phosphate_atoms["DPPC"]]
			distance = CartesianToolkit.vecdist([phosphate[0],phosphate[1],0],[backbone_com[0],backbone_com[1],0])
			target_bin = int(distance/0.5)
			thickness = lipid[phosphate_atoms["DPPC"]][2] - bilayer_midpoint
			#ignore phosphates outside of the boxes (ie at values greater than the radius of the box)
			if thickness < 0: thickness = 0 - thickness
			if target_bin < len(upper_leaflet_bins):
				upper_leaflet_bins[target_bin].append(thickness)
			else:
				continue

		for lipid in lower_leaflet:
			phosphate = lipid[phosphate_atoms["DPPC"]]
			distance = CartesianToolkit.vecdist([phosphate[0],phosphate[1],0],[backbone_com[0],backbone_com[1],0])
			target_bin = int(distance/0.5)
			thickness = lipid[phosphate_atoms["DPPC"]][2] - bilayer_midpoint
			#ignore phosphates outside of the boxes (ie at values greater than the radius of the box)
			if thickness < 0: thickness = 0 - thickness
			if target_bin < len(lower_leaflet_bins):
				lower_leaflet_bins[target_bin].append(thickness)
			else:
				continue
		
		#Calculate rotation
		
		##Define start and end turns- 3 residues of helical backbone. Use to calculate the helix axis vector and the helix angle wrt z
		try:
			start_positions = index_cartesians(groups["B*"][1:4],trajectory.cartesian)
			end_positions = index_cartesians(groups["B*"][-4:-1],trajectory.cartesian)
		except:
			start_positions = index_cartesians(groups["B*_CA"][1:4],trajectory.cartesian)
			end_positions = index_cartesians(groups["B*_CA"][-4:-1],trajectory.cartesian)
		
		start_com = CartesianToolkit.center(start_positions)
		end_com = CartesianToolkit.center(end_positions)
		
		helix_axis = CartesianToolkit.vecsub(end_com,start_com)
		helix_angle = CartesianToolkit.vecangle(helix_axis,[0,0,1])
		
		#Alternative technique for helix rotation calculation
		#Here we calculate a reference vector in the direction of the tilt of the helix
		
		ref_axis = CartesianToolkit.veccross(CartesianToolkit.veccross(CartesianToolkit.vecnorm(helix_axis),[0,0,1]),CartesianToolkit.vecnorm(helix_axis))
		
		#print helix_angle to the new tilt file
		print >> tilt_file, frame, 180*helix_angle/math.pi
		##The reference is the 3rd residue backbone atom
		reference_residue = start_positions[1]
		
		##Calculate the unit vector of the residue with the helix.
		##S/E are the start/end, R is the reference, N is the reference distance on SE
		##       ---------R
		##     /          |
		##   S------------N-----------E
		
		SR = CartesianToolkit.vecsub(reference_residue,start_com)
		SE = CartesianToolkit.vecsub(end_com,start_com)
		SN = CartesianToolkit.vecscale(SE,(CartesianToolkit.veclength(SR)/CartesianToolkit.veclength(SE) * CartesianToolkit.vecscaler(CartesianToolkit.vecnorm(SE),CartesianToolkit.vecnorm(SR))))
		point_on_axis = CartesianToolkit.vecadd(SN,start_com)
		vector_from_axis = CartesianToolkit.vecsub(reference_residue,point_on_axis)
		unit_vector_from_axis = CartesianToolkit.vecnorm(vector_from_axis)
		
		#Calculate new angle
		tng_screw_angle = CartesianToolkit.vecangle(ref_axis,unit_vector_from_axis)
		#if vector_from_xaxis[1] < 0: tng_screw_angle = 0 - tng_screw_angle
		
		updown = CartesianToolkit.veccross(ref_axis,unit_vector_from_axis)
		helix_dot_rehelix = CartesianToolkit.vecangle(updown,SE)
		
		if ( helix_dot_rehelix < math.pi/2 and helix_dot_rehelix >= 0 )or helix_dot_rehelix <-math.pi/2:
			tng_screw_angle = 0 - tng_screw_angle
			#print "Same     ", helix_dot_rehelix*180/math.pi
		else:
			#print "Different", helix_dot_rehelix*180/math.pi
			pass

		##Correct rotation- 0' should always point to upper leaflet (the leaflet with the 1st protein residue in)
		if top_comz < bilayer_midpoint:
			tng_screw_angle += math.pi
			tng_screw_angle = CartesianToolkit.wrapangle(tng_screw_angle)

		
		print >> tng_rotation_file, frame, 180*tng_screw_angle/math.pi
		
		##Homebrew rotation matrix calculations
		##To replicate VMD results, first rotate the matrix around y, then around z
		
		#Now safe against div by 1 errors
		if helix_axis[1] != 0.:
			rotation_matrix_around_z = CartesianToolkit.rotation_matrix_to_x(CartesianToolkit.vecnorm([helix_axis[0],helix_axis[1],0]))
			partially_rotated_helix_vector = CartesianToolkit.rotate_by_matrix(helix_axis,rotation_matrix_around_z)
		else:
			partially_rotated_helix_vector = helix_axis
			rotation_matrix_around_z = None
		
		
		if partially_rotated_helix_vector[2] != 0:
			rotation_matrix_around_y = CartesianToolkit.rotation_matrix_to_x(CartesianToolkit.vecnorm([partially_rotated_helix_vector[0], 0, partially_rotated_helix_vector[2] ]))
		else:
			rotation_matrix_around_y = None
		
		if rotation_matrix_around_z:
			partially_rotated_screw_vector = CartesianToolkit.rotate_by_matrix(unit_vector_from_axis,rotation_matrix_around_z)
		else:
			partially_rotated_screw_vector = unit_vector_from_axis
		if rotation_matrix_around_y:
			vector_from_xaxis = CartesianToolkit.rotate_by_matrix(partially_rotated_screw_vector,rotation_matrix_around_y)
		else:
			vector_from_xaxis = partially_rotated_screw_vector

		#print helix_axis, partially_rotated_helix_vector, vector_from_xaxis

		##Calculate rotation
		
		screw_angle = CartesianToolkit.vecangle(vector_from_xaxis,[0,0,1])
		
		if vector_from_xaxis[1] < 0: screw_angle = 0 - screw_angle
		
		##Correct rotation- 0' should always point to upper leaflet (the leaflet with the 1st protein residue in)
		if top_comz < bilayer_midpoint:
			screw_angle += math.pi
			screw_angle = CartesianToolkit.wrapangle(screw_angle)
			
		print >> rotation_file, frame, 180*screw_angle/math.pi
		
		#Test for helical insertion
		
		#print "Testing", "Comz =", bkbone_bilayer_distance, "(1/-1) Angle =", helix_angle, "(<",5*math.pi/12,")"
		if ( bkbone_bilayer_distance < 1 and bkbone_bilayer_distance > -1 ) and ( helix_angle < (5*math.pi/12) or helix_angle > (7*math.pi/12)):
			helix_insertion += 1
			#print True
		
		#Update frame
		frame += 1

	#For the final frame, test to see what the bilayer looks like
			
	##Most important bilayer tests are the separation of water compartments and the degree of water lipid mixing
	##At present I test the bilayer properties by dividing the box into 5 Angstrom cells
	
	#Add extra tests for bilayers in xy, xz, zy
	
	bilayers = 0
	unusual_mixing = 0
		
	#box_x = trajectory.box[0]
	#box_y = trajectory.box[4]
	box_z = trajectory.box[bilayer_normal*4]
	
	#box = [box_x,box_y,box_z]
	
	waters = [0 for i in range(int(box_z/0.5)+1)]
	lipids = [0 for i in range(int(box_z/0.5)+1)]
	
	#waters = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
	#lipids = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
		
	for water_particle in water_positions:
		bin = int(water_particle[2]/0.5)
		#print water_particle[2], bin, len(waters), box_z, trajectory.box
		if bin >= len(waters):
				#print "Wrapping", water_particle[2], bin, bin - len(waters)
				bin -= len(waters)
		elif bin <= -len(waters):
				#print "Wrapping", water_particle[2], bin, bin + len(waters)
				bin += len(waters)
		waters[bin] += 1
		
	for lipid_atom in lipid_positions:
		bin = int(lipid_atom[2]/0.5)
		if bin >= len(lipids):
				#print "Wrapping", water_particle[2], bin, bin - len(lipids)
				bin -= len(lipids)
		elif bin <= -len(lipids):
				#print "Wrapping", water_particle[2], bin, bin + len(lipids)
				bin += len(lipids)
		lipids[bin] += 1
	water_lipid_cells = 0
	water_lipid_mix = zip(waters,lipids)
	for cell in water_lipid_mix:
		if cell[0] > 0 and cell[1] > 0: water_lipid_cells += 1
	if water_lipid_cells > 11:
		unusual_mixing = True

	if waters[0]: previous_state = 1
	else: previous_state = 0
	for cell in waters[1:]:
		if cell: current_state = 1
		else: current_state = 0
		if current_state < previous_state:
			bilayers += 1
		previous_state = current_state
		
	efficiency_file = open("efficiency.dat","w")
	print >> efficiency_file, "Number of bilayers =", bilayers
	if unusual_mixing: print >> efficiency_file, "Unusual water/lipid mixing. Possible deformation."
	print >> efficiency_file, "Helix insertion =", 100*helix_insertion/frame, "%"

		
	#Finalise- close files and print lipid positions
	for bin in upper_leaflet_bins:
		for item in bin:
			print >> upper_leaflet_file, 10*item,
		print >> upper_leaflet_file, ""
	for bin in lower_leaflet_bins:
		for item in bin:
			print >> lower_leaflet_file, 10*item,
		print >> lower_leaflet_file, ""
	lower_leaflet_file.close()
	upper_leaflet_file.close()
	bilayer_position_file.close()
	rotation_file.close()
	efficiency_file.close()
	tilt_file.close()

def test_bilayer_normal(structure="t_0_center.pdb",index="final_system.ndx",trajectory="t_0.xtc"):
	endpoint = pdbio.read_pdb(structure)
	groups = ndxio.read_ndx(index)
	endpoint.convert_coordinates()
	
	box= [endpoint.box[0], endpoint.box[4], endpoint.box[8]]
	print box
	
	waters = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
	lipids = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
	
	#lipid_positions = index_cartesians(groups["DPPC"],endpoint.cartesian)
	lipid_positions = []
	for lipid in allowed_lipids:
			try:
				lipid_positions += index_cartesians(groups[lipid],trajectory.cartesian)
				#Need an additional list of lipid groups to calculate the lipid list
			except:
				pass
	water_positions = index_cartesians(groups["W"],endpoint.cartesian)
	
	bilayers = [0, 0, 0]
	#waters = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
	#lipids = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
		
	for  i in range(3):
		for water_particle in water_positions:
			bin = int(water_particle[i]/0.5)
			if bin >= len(waters[i]):
				print "Wrapping", water_particle[i], bin, bin - len(waters[i])
				bin -= len(waters[i])
			elif bin <= -len(waters[i]):
				print "Wrapping", water_particle[i], bin, bin + len(waters[i])
				bin += len(waters[i])
			
			#print water_particle[i]
			waters[i][bin] += 1
		print waters[i]
		'''
		for lipid_atom in lipid_positions:
			bin = int(lipid_atom[i]/0.5)
			lipids[i][bin] += 1
		'''
		if waters[i][0]: previous_state = 1
		else: previous_state = 0
		for cell in waters[i][1:]:
			if cell: current_state = 1#; print "Dim", i, "Wet!", cell
			else: current_state = 0#; print "Dim", i, "Dry!", cell
			if current_state < previous_state:
				#print "We have a winner", i
				bilayers[i] += 1
			previous_state = current_state
	top_dim = bilayers.index(max(bilayers))
	print bilayers, top_dim
	'''
	if top_dim == 2:
		print "XY bilayer"
		#xy bilayer- do nothing
		pass
	elif top_dim == 1:
		print "XZ bilayer"
		#xz bilayer- rotate 90' in x
		os.rename("t_0.xtc","initial_t_0.xtc")
		trjconv_command = gromacs + "trjconv -f initial_t_0.xtc -o t_0.xtc -rotate 90 0 0"
		os.system(trjconv_command)
	elif top_dim == 0:
		print "YZ bilayer"
		#yz bilayer- rotate 90' in y
		os.rename("t_0.xtc","initial_t_0.xtc")
		trjconv_command = gromacs + "trjconv -f initial_t_0.xtc -o t_0.xtc -rotate 0 90 0"
		os.system(trjconv_command)
	'''
	return top_dim
	
def visualise(initial_sequence,charge,bias=True,number_lipid_types=1):
	'''gromacs analysis'''
	
	length = len(initial_sequence)
	
	if length %2 == 0:
	        middle_residues = str(length/2 - 1) + "-" + str(length/2 +2)
	else:
	        middle_residues = str(length/2 - 1) + "-" + str(length/2 +1)
	
	group_mod = number_lipid_types - 1
	if charge != 0: group_mod += 1
	
	'''
	if charge != 0:
	        make_ndx_command = 'echo "a b* CA\nr 1-7\n r ' + str(length-6) + '-' + str(length) +'\n r ' + middle_residues + '\n 16 & 17\n 16 & 18\n 16 & 19\n q\n" | '  + gromacs + 'make_ndx -f t_0 -o system.ndx'
	        g_bundle_command = "echo '20\n21\n22\n'| " + gromacs + "g_bundle -f t_0 -s t_0 -na 1 -z -ok -n system.ndx "
	else:
	        make_ndx_command = 'echo "a b* CA\nr 1-7\n r ' + str(length-6) + '-' + str(length) +'\n r ' + middle_residues + '\n 15 & 16\n 15 & 17\n 15 & 18\n q\n" | '  + gromacs + 'make_ndx -f t_0 -o system.ndx'
	        g_bundle_command = "echo '19\n20\n21\n'| " + gromacs + "g_bundle -f t_0 -s t_0 -na 1 -z -ok -n system.ndx "
	'''
	base = 15 + group_mod
	make_ndx_grps = str(base) + " & " + str(base + 1) + "\n " + str(base) + " & " + str(base + 2) + "\n " + str(base) + " & " + str(base + 3) + "\n "
	make_ndx_command = 'echo "a b* CA\nr 1-7\n r ' + str(length-6) + '-' + str(length) +'\n r ' + middle_residues + '\n ' +make_ndx_grps+ 'q\n" | '  + gromacs + 'make_ndx -f t_0 -o system.ndx'
	g_bundle_grps = str(base + 4) + "\n" + str(base + 5) + "\n" + str(base + 6) + "\n"
	g_bundle_command = "echo '" + g_bundle_grps + "'| " + gromacs + "g_bundle -f t_0 -s t_0 -na 1 -z -ok -n system.ndx "

	
	os.system(make_ndx_command)
	#print make_ndx_command
	os.system(g_bundle_command)
	
	AutomatedPlot.tilt_data_plot("bun_tilt.xvg","tilt_V_time.png",initial_sequence)
	AutomatedPlot.hist_tilt_data_plot("bun_tilt.xvg","hist_tilt.png",initial_sequence)
	
	
	AutomatedPlot.tilt_data_plot("bun_kink.xvg","kink_V_time.png",initial_sequence)
	AutomatedPlot.hist_tilt_data_plot("bun_kink.xvg","hist_kink.png",initial_sequence)

	'''VMD analysis'''
	
	#If unbiased, test to find the normal of the bilayer, and rotate the trajectory to account for it
	if not bias:
		editconf_command = "echo '1\n0\n' | " + gromacs + "editconf -f t_0.gro -c -o t_0_center.pdb "
		make_ndx_command = "echo 'q\n' | " +gromacs + "make_ndx -f t_0.gro -o final_system.ndx"
		os.system(editconf_command)
		os.system(make_ndx_command)
		bilayer_normal = test_bilayer_normal(structure="t_0_center.pdb",index="final_system.ndx",trajectory="t_0.xtc")
		print "Bilayer normal", bilayer_normal
	else:
		bilayer_normal = 2
	
	if bilayer_normal == 2:
		AutomatedPlot.tilt_data_plot("bun_tilt.xvg","tilt_V_time.png",initial_sequence)
		AutomatedPlot.hist_tilt_data_plot("bun_tilt.xvg","hist_tilt.png",initial_sequence)
	else:
		AutomatedPlot.tilt_data_plot("bun_tilt.xvg","tilt_V_time.png",initial_sequence,modifier=90.)
		AutomatedPlot.hist_tilt_data_plot("bun_tilt.xvg","hist_tilt.png",initial_sequence,modifier=90.)

	#center the trajectory to get all residues in the box, and the helix in the middle of the box
	trjconv_center_command = "echo '1\n0\n' | " + gromacs + "trjconv -f t_0.xtc -s t_0.tpr -pbc mol -center -o center.xtc "
	os.system(trjconv_center_command)
	if bias:
		#Nasty hack due to tempremental problem for cvs trjconv
		os.system("touch electroneg.dat; touch elements.dat")
		trjconv_xy_command = "echo '1\n0\n' | " + gromacs + "trjconv_xy -f center.xtc -s em.gro -o xy_fit.xtc -fit rotxy+transxy"
		os.system(trjconv_xy_command)
	else:
		#Don't fit to xy if unbiased- there is no way of doing this at present
		os.rename("center.xtc","xy_fit.xtc")
	
	'''
	vmd_command = vmd + " em.gro xy_fit.xtc t_0.gro -dispdev text -eofexit <" + configuration['vmd_state_files'] + "automated_analysis.tcl > vmd.log"
	os.system(vmd_command)
	'''
	simulation_statistics(bilayer_normal)
	
	AutomatedPlot.hist_com_difference_plot("bilayer_position.dat","bilayer_position.png",initial_sequence)
	
	AutomatedPlot.bilayer_deformation_refined("up.dat","down.dat","both-bully.png",graph_title=initial_sequence)
	
	AutomatedPlot.rotation_plot("rot.dat","rotation.png",initial_sequence)
