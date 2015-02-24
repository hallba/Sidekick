#!/usr/bin/python

from __future__ import with_statement

import sys,os,pdbio,numpy,math

from HAConf import configuration,programs,debug_status
from LipidBox import add_solvent_dppc
from MDPGenerator import replace_seed, replace_steps
import AutomatedPlot
import xtcio,ndxio,CartesianToolkit
from GenerateCGSystem import allowed_lipids
import MatrixIO

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

class AnalyseSystem:
	def __init__(self):
		pass
	def AnalysisPipeline(self):
		self.process_trajectory()
		self.pre_analyse()
		self.analytical_loop()
		self.final_analysis()
		self.plot_graphs()
	def process_trajectory(self):
		'''Do some fitting'''
		pass
	def pre_analyse(self):
		'''Any important details about the simulation we need to know before starting the analysis?'''
		pass
	def analytical_loop(self):
		self.trajectory = xtcio.read_xtc(self.analytic_trajfile)
		self.groups = ndxio.read_ndx(self.analytic_ndxfile)
		if "PW" in self.groups.keys():
			self.groups["W"] = self.groups["PW"]
		self.frame = 0
		#self.helix_insertion = 0.
		self.lipid_map = self.lipid_list_from_topology(self.system.topology_name)
		while(self.trajectory.next_frame()):
			self.analysis_core()
	def analysis_core(self):
		pass
	def final_analysis(self):
		pass


	def index_cartesians(self,index_list,frame):
		return [frame[item] for item in index_list]

	def lipid_list_breakdown_mix(self,lipid_map,lipid_positions):
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

	def test_bilayer_normal(self,structure="t_0_center.pdb",index="final_system.ndx",trajectory="t_0.xtc"):
		endpoint = pdbio.read_pdb(structure)
		groups = ndxio.read_ndx(index)
		if "PW" in groups.keys():
			groups["W"] = groups["PW"]
		endpoint.convert_coordinates()
		
		box= [endpoint.box[0], endpoint.box[4], endpoint.box[8]]
		print box
		
		waters = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
		lipids = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
		
		#lipid_positions = index_cartesians(groups["DPPC"],endpoint.cartesian)
		self.lipid_positions = []
		for lipid in allowed_lipids:
				try:
					self.lipid_positions += self.index_cartesians(groups[lipid],trajectory.cartesian)
					#Need an additional list of lipid groups to calculate the lipid list
				except:
					pass
		self.water_positions = self.index_cartesians(groups["W"],endpoint.cartesian)
		
		bilayers = [0, 0, 0]
		#waters = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
		#lipids = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
			
		for  i in range(3):
			for water_particle in self.water_positions:
				bin = int(water_particle[i]/0.5)
				if bin >= len(waters[i]):
					print "Wrapping", water_particle[i], bin, bin - len(waters[i])
					bin -= len(waters[i])
				elif bin <= -len(waters[i]):
					print "Wrapping", water_particle[i], bin, bin + len(waters[i])
					bin += len(waters[i])
				
				#print water_particle[i]
				waters[i][bin] += 1
			#print waters[i]
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
		#print bilayers, top_dim
		return top_dim

	def bilayer_position(self,backbone_com,bilayer_midpoint,top_comz):
		#Calculate bilayer distance
		bkbone_bilayer_distance = backbone_com[2] - bilayer_midpoint
		if top_comz < backbone_com[2]:
			bkbone_bilayer_distance = 0 - bkbone_bilayer_distance
		return bkbone_bilayer_distance

	def bilayer_center(self,individual_lipid_positions,lipid_com):
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

	def lipid_deformation(self,bilayer_midpoint,upper_leaflet,lower_leaflet,backbone_com):
		if self.frame == 0:
			#Bin size is 5 Angstroms, therefore there are the same number of bins as nm in the total box length
			if self.bilayer_normal == 2:
				self.upper_leaflet_bins = [ [] for i in range(int(self.trajectory.box[0])) ]
				self.lower_leaflet_bins = [ [] for i in range(int(self.trajectory.box[0])) ]
			else:
				self.upper_leaflet_bins = [ [] for i in range(int(self.trajectory.box[8])) ]
				self.lower_leaflet_bins = [ [] for i in range(int(self.trajectory.box[8])) ]
		for lipid in upper_leaflet:
			phosphate = lipid[phosphate_atoms["DPPC"]]
			distance = CartesianToolkit.vecdist([phosphate[0],phosphate[1],0],[backbone_com[0],backbone_com[1],0])
			target_bin = int(distance/0.5)
			thickness = lipid[phosphate_atoms["DPPC"]][2] - bilayer_midpoint
			#ignore phosphates outside of the boxes (ie at values greater than the radius of the box)
			if thickness < 0: thickness = 0 - thickness
			if target_bin < len(self.upper_leaflet_bins):
				self.upper_leaflet_bins[target_bin].append(thickness)
			else:
				continue

		for lipid in lower_leaflet:
			phosphate = lipid[phosphate_atoms["DPPC"]]
			distance = CartesianToolkit.vecdist([phosphate[0],phosphate[1],0],[backbone_com[0],backbone_com[1],0])
			target_bin = int(distance/0.5)
			thickness = lipid[phosphate_atoms["DPPC"]][2] - bilayer_midpoint
			#ignore phosphates outside of the boxes (ie at values greater than the radius of the box)
			if thickness < 0: thickness = 0 - thickness
			if target_bin < len(self.lower_leaflet_bins):
				self.lower_leaflet_bins[target_bin].append(thickness)
			else:
				continue

	def helix_rotation(self,backbone_positions,top_comz,bilayer_midpoint):
		start_positions = backbone_positions[1:4]
		end_positions = backbone_positions[-4:-1]
		
		start_com = CartesianToolkit.center(start_positions)
		end_com = CartesianToolkit.center(end_positions)
		
		helix_axis = CartesianToolkit.vecsub(end_com,start_com)
		helix_angle = CartesianToolkit.vecangle(helix_axis,[0,0,1])
		
		#Here we calculate a reference vector in the direction of the tilt of the helix
		
		ref_axis = CartesianToolkit.veccross(CartesianToolkit.veccross(CartesianToolkit.vecnorm(helix_axis),[0,0,1]),CartesianToolkit.vecnorm(helix_axis))

		#Internal helix tilt
		#helix_angle = 180*helix_angle/math.pi
		
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
		
		tng_screw_angle = CartesianToolkit.vecangle(ref_axis,unit_vector_from_axis)

		updown = CartesianToolkit.veccross(ref_axis,unit_vector_from_axis)
		helix_dot_rehelix = CartesianToolkit.vecangle(updown,SE)

		if ( helix_dot_rehelix < math.pi/2 and helix_dot_rehelix >= 0 )or helix_dot_rehelix <-math.pi/2:
			tng_screw_angle = 0 - tng_screw_angle
			#print "Same     ", helix_dot_rehelix*180/math.pi
		else:
			#print "Different", helix_dot_rehelix*180/math.pi
			pass
		
		if top_comz < bilayer_midpoint:
			tng_screw_angle += math.pi
			tng_screw_angle = CartesianToolkit.wrapangle(tng_screw_angle)
		
		return (helix_angle, tng_screw_angle)
		
	def lipid_list_from_topology(self,filename):
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
class AnalyseHelix(AnalyseSystem):
	def __init__(self, system,
			input_name="t_0",
			original_system = "em.gro",
			bilayer_position_file="bilayer_position.dat",
			upper_leaflet_file = "up.dat",
			lower_leaflet_file = "down.dat",
			rotation_file = "rot.dat",
			tilt_file = "tilt.dat",
			efficiency_file = "efficiency.dat"):
			
		self.system = system
		
		self.input_trajectory=input_name + ".xtc"
		self.input_name = input_name
		self.original_system = original_system
		self.bilayer_position_file=bilayer_position_file
		self.upper_leaflet_file = upper_leaflet_file
		self.lower_leaflet_file = lower_leaflet_file
		self.rotation_file = rotation_file
		self.tilt_file = tilt_file
		self.efficiency_file = efficiency_file
		
		backup_file(self.bilayer_position_file)
		backup_file(self.upper_leaflet_file)
		backup_file(self.lower_leaflet_file)
		backup_file(self.rotation_file)
		backup_file(self.tilt_file)
		backup_file(self.efficiency_file)
		
		self.analytic_trajfile="xy_fit.xtc"
		self.analytic_ndxfile="system.ndx"
		
		self.helix_insertion = 0.
		
		self.AnalysisPipeline()
		
	
	def process_trajectory(self):
			
			
			length = len(self.system.initial_sequence)
	
			if length %2 == 0:
			        middle_residues = str(length/2 - 1) + "-" + str(length/2 +2)
			else:
			        middle_residues = str(length/2 - 1) + "-" + str(length/2 +1)
			
			group_mod = self.system.number_lipid_types - 1
			if self.system.charge != 0: group_mod += 1
			
			base = 15 + group_mod
			make_ndx_grps = str(base) + " & " + str(base + 1) + "\n " + str(base) + " & " + str(base + 2) + "\n " + str(base) + " & " + str(base + 3) + "\n "
			make_ndx_command = 'echo "a b* CA\nri 1-7\n ri ' + str(length-6) + '-' + str(length) +'\n ri ' + middle_residues + '\n ' +make_ndx_grps+ 'q\n" | '  + gromacs + 'make_ndx -f t_0 -o system.ndx'
			g_bundle_grps = str(base + 4) + "\n" + str(base + 5) + "\n" + str(base + 6) + "\n"
			g_bundle_command = "echo '" + g_bundle_grps + "'| " + gromacs + "g_bundle -f " + self.input_name + " -s " + self.input_name + " -na 1 -z -ok -n system.ndx "
		
			
			os.system(make_ndx_command)
			#print make_ndx_command
			os.system(g_bundle_command)
		
			#center the trajectory to get all residues in the box, and the helix in the middle of the box
			trjconv_center_command = "echo '1\n0\n' | " + gromacs + "trjconv -f " + self.input_name + ".xtc -s " + self.input_name + ".tpr -pbc mol -center -o center.xtc "
			os.system(trjconv_center_command)
			if self.system.bias:
				#Nasty hack due to tempremental problem for cvs trjconv
				os.system("touch electroneg.dat; touch elements.dat")
				trjconv_xy_command = "echo '1\n0\n' | " + gromacs + "trjconv -f center.xtc -s " + self.original_system + " -o xy_fit.xtc -fit rotxy+transxy"
				os.system(trjconv_xy_command)
			else:
				#Don't fit to xy if unbiased- there is no way of doing this at present
				os.rename("center.xtc","xy_fit.xtc")
				
			os.rename("xy_fit.xtc",self.analytic_trajfile)
			os.rename("system.ndx",self.analytic_ndxfile)
			#self.analytic_trajfile = analytic_trajfile
			#self.analytic_ndxfile = analytic_ndxfile

	def pre_analyse(self):
			#If unbiased, test to find the normal of the bilayer, and rotate the trajectory to account for it
			if not self.system.bias:
				editconf_command = "echo '1\n0\n' | %seditconf -f %s.gro -c -o %s_center.pdb " %(gromacs, self.input_name, self.input_name)
				make_ndx_command = "echo 'q\n' | %smake_ndx -f %s.gro -o final_system.ndx"  %(gromacs, self.input_name)
				os.system(editconf_command)
				os.system(make_ndx_command)
				self.bilayer_normal = self.test_bilayer_normal(structure=self.input_name+"_center.pdb",index="final_system.ndx",trajectory=self.input_name+".xtc")
				print "Bilayer normal", self.bilayer_normal
			else:
				self.bilayer_normal = 2

	def analysis_core(self):
		self.trajectory.convert_coordinates()

		#Do I need to rotate the cartesians to account for non-xy bilayers.
		if self.bilayer_normal==2:
			pass
		elif self.bilayer_normal == 1:
			#rotate by 90' in x- rotate y to x, then z to x. There must be a better way
			bilayer_rotmat_first = CartesianToolkit.rotation_matrix_to_x([0,1,0])
			bilayer_rotmat_second = CartesianToolkit.rotation_matrix_to_x([0,0,1])
			#for cart in trajectory.cartesian:
			#	cart = CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat_first)
			#	cart = CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat_second)
			self.trajectory.cartesian = [CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat_first) for cart in self.trajectory.cartesian]
			self.trajectory.cartesian = [CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat_second) for cart in self.trajectory.cartesian]
		elif self.bilayer_normal == 0:
			#rotate by 90' in y: x becomes z, z becomes -x
			bilayer_rotmat = CartesianToolkit.rotation_matrix_to_x([0,0,1])
			#for cart in trajectory.cartesian:
			#	cart = CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat)
			self.trajectory.cartesian = [CartesianToolkit.rotate_by_matrix(cart,bilayer_rotmat) for cart in self.trajectory.cartesian]
		else:
			print "Bilayer normal not feasible. Aborting"
			sys.exit()

		#Calculate atom positions of interest
		try:
			backbone_positions = self.index_cartesians(self.groups["B*"],self.trajectory.cartesian)
		except:
			backbone_positions = self.index_cartesians(self.groups["B*_CA"],self.trajectory.cartesian)
		self.lipid_positions = []
		for lipid in allowed_lipids:
			try:
				self.lipid_positions += self.index_cartesians(self.groups[lipid],self.trajectory.cartesian)
				#Need an additional list of lipid groups to calculate the lipid list
			except:
				pass

		self.water_positions = self.index_cartesians(self.groups["W"],self.trajectory.cartesian)
		top_comz = backbone_positions[0][2]
		
		#Calculate centers of interest
		backbone_com = CartesianToolkit.center(backbone_positions)
		lipid_com = CartesianToolkit.center(self.lipid_positions)
		#(bilayer_midpoint,upper_leaflet,lower_leaflet) = bilayer_center(lipid_list_breakdown(lipid_positions,lipid_atoms["DPPC"]),lipid_com)
		(bilayer_midpoint,upper_leaflet,lower_leaflet) = self.bilayer_center(self.lipid_list_breakdown_mix(self.lipid_map,self.lipid_positions),lipid_com)
		
		if top_comz < bilayer_midpoint:
			[upper_leaflet,lower_leaflet] = [lower_leaflet,upper_leaflet]
			
		#Now onto the core analysis: position, tilt, screw rotation, bilayer deformation
		
		#Bilayer position
		bkbone_bilayer_distance = self.bilayer_position(backbone_com,bilayer_midpoint,top_comz)
		
		with open(self.bilayer_position_file,"a") as output:
			print >> output, self.frame, 10*bkbone_bilayer_distance
		
		#Lipid deformation
		self.lipid_deformation(bilayer_midpoint,upper_leaflet,lower_leaflet,backbone_com)
		
		#Rotation
		##Uses the new algorithm avoiding manual rotations
		
		(internal_tilt, screw_angle) = self.helix_rotation(backbone_positions,top_comz,bilayer_midpoint)
		
		with open(self.rotation_file,"a") as output:
			print >> output, self.frame, 180*screw_angle/math.pi
		with open(self.tilt_file,"a") as output:
			print >> output, self.frame, 180*internal_tilt/math.pi
			
		if ( bkbone_bilayer_distance < 1 and bkbone_bilayer_distance > -1 ) and ( internal_tilt < (5*math.pi/12) or internal_tilt > (7*math.pi/12)):
			self.helix_insertion += 1
			
		self.frame += 1

	def final_analysis(self):
		#For the final frame, test to see what the bilayer looks like
				
		##Most important bilayer tests are the separation of water compartments and the degree of water lipid mixing
		##At present I test the bilayer properties by dividing the box into 5 Angstrom cells
		
		#Add extra tests for bilayers in xy, xz, zy
		
		bilayers = 0
		unusual_mixing = 0
			
		#box_x = trajectory.box[0]
		#box_y = trajectory.box[4]
		box_z = self.trajectory.box[self.bilayer_normal*4]
		
		#box = [box_x,box_y,box_z]
		
		waters = [0 for i in range(int(box_z/0.5)+1)]
		lipids = [0 for i in range(int(box_z/0.5)+1)]
		
		#waters = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
		#lipids = [[0 for i in range(int(box_cart/0.5)+1)] for box_cart in box]
			
		for water_particle in self.water_positions:
			bin = int(water_particle[2]/0.5)
			#print water_particle[2], bin, len(waters), box_z, trajectory.box
			if bin >= len(waters):
					#print "Wrapping", water_particle[2], bin, bin - len(waters)
					bin -= len(waters)
			elif bin <= -len(waters):
					#print "Wrapping", water_particle[2], bin, bin + len(waters)
					bin += len(waters)
			waters[bin] += 1
			
		for lipid_atom in self.lipid_positions:
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
			
		with open(self.efficiency_file,"w") as efficiency_file:
			print >> efficiency_file, "Number of bilayers =", bilayers
			if unusual_mixing: print >> efficiency_file, "Unusual water/lipid mixing. Possible deformation."
			print >> efficiency_file, "Helix insertion =", 100*self.helix_insertion/self.frame, "%"
	
		#Finalise- close files and print lipid positions
		with open(self.upper_leaflet_file,"w") as output:
			for bin in self.upper_leaflet_bins:
				for item in bin:
					print >> output, 10*item,
				print >> output, ""
		with open(self.lower_leaflet_file,"w") as output:
			for bin in self.lower_leaflet_bins:
				for item in bin:
					print >> output, 10*item,
				print >> output, ""
	
	def plot_graphs(self):
			#First plot the gromacs graphs
			if self.bilayer_normal == 2:
				AutomatedPlot.tilt_data_plot("bun_tilt.xvg","tilt_V_time.png",self.system.initial_sequence)
				AutomatedPlot.hist_tilt_data_plot("bun_tilt.xvg","hist_tilt.png",self.system.initial_sequence)
			else:
				AutomatedPlot.tilt_data_plot("bun_tilt.xvg","tilt_V_time.png",self.system.initial_sequence,modifier=90.)
				AutomatedPlot.hist_tilt_data_plot("bun_tilt.xvg","hist_tilt.png",self.system.initial_sequence,modifier=90.)

			AutomatedPlot.hist_com_difference_plot(self.bilayer_position_file,"bilayer_position.png",self.system.initial_sequence)
			AutomatedPlot.bilayer_deformation_refined(self.upper_leaflet_file,self.lower_leaflet_file,"both-bully.png",graph_title=self.system.initial_sequence)
			AutomatedPlot.rotation_plot(self.rotation_file,"rotation.png",self.system.initial_sequence)

class AnalyseHelixDimer(AnalyseHelix):
	def __init__(self, system,
			input_trajectory="t_0.xtc",
			bilayer_position_file=".posi.dat",
			rotation_file = ".rota.dat",
			tilt_file = ".tilt.dat",
			distance_file = "distance.dat",
			efficiency_file = "efficiency.dat",
			cc_eigenvalue_file = "cc_eigenvalue.dat",
			ca_eigenvalue_file = "ca_eigenvalue.dat",
			handedness_file = "handedness.dat",
			vector_relation_file = "vector.dat",
			contact_file = "contact.dat",
			tip_distance_file = "hairpin.dat"
			):
			
		self.system = system
		
		self.input_trajectory=input_trajectory
		self.bilayer_position_file=bilayer_position_file
		self.rotation_file = rotation_file
		self.tilt_file = tilt_file
		self.efficiency_file = efficiency_file
		self.distance_file = distance_file
		self.cc_eigenvalue_file = cc_eigenvalue_file
		self.ca_eigenvalue_file = ca_eigenvalue_file
		self.handedness_file = handedness_file
		self.vector_relation_file = vector_relation_file
		self.contact_file = contact_file
		self.tip_distance_file = tip_distance_file
		
		backup_file(self.bilayer_position_file)
		backup_file(self.rotation_file)
		backup_file(self.tilt_file)
		backup_file(self.efficiency_file)
		backup_file(self.distance_file)
		backup_file(self.cc_eigenvalue_file)
		backup_file(self.ca_eigenvalue_file)
		backup_file(self.handedness_file)
		backup_file(self.vector_relation_file)
		backup_file(self.contact_file)
		backup_file(self.tip_distance_file)
		
		self.analytic_trajfile="xy_fit.xtc"
		self.analytic_ndxfile="system.ndx"
		
		self.helix_insertion_A = 0.
		self.helix_insertion_B = 0.
		
		self.AnalysisPipeline()

	def test_handedness(self,helix_A_backbone,helix_B_backbone,backbone_A_com,backbone_B_com):
		def define_helix_axis(backbone_positions):
			start_positions = backbone_positions[1:4]
			end_positions = backbone_positions[-4:-1]
			
			start_com = CartesianToolkit.center(start_positions)
			end_com = CartesianToolkit.center(end_positions)
			
			helix_axis = CartesianToolkit.vecsub(end_com,start_com)
			return helix_axis
			
		A_B_comvector = CartesianToolkit.vecsub(backbone_A_com,backbone_B_com)
		axis_A = CartesianToolkit.vecnorm(define_helix_axis(helix_A_backbone))
		axis_B = CartesianToolkit.vecnorm(define_helix_axis(helix_B_backbone))
				
		if self.system.parallel != True:
			axis_B = [-item for item in axis_B]
		
		AB_Cross = CartesianToolkit.veccross(axis_A,axis_B)
		
		angle = CartesianToolkit.vecangle(axis_A,axis_B)
		
		shortest_distance = CartesianToolkit.vecdot(AB_Cross,CartesianToolkit.vecnorm(A_B_comvector))
		
		if shortest_distance > 0:
			angle = -angle
		
		return angle

	def gnm(self,cutoff=0.7):
		backbone_particles = self.groups["B*_CA"]
		nresidues = len(backbone_particles)
		natoms = len(self.protein_positions)
		matrix = numpy.zeros((nresidues,nresidues),"float")
		
		for i in range(nresidues):
			for j in range(i+1,nresidues):
				distance =CartesianToolkit.veclength(CartesianToolkit.vecsub(self.backbone_positions[i],self.backbone_positions[j]))
				#print i,j, distance, cutoff
				if distance < cutoff:
					#print "contact!"
					matrix[i][j] = matrix[i][j] - 1
					matrix[j][i] = matrix[j][i] - 1
					matrix[j][j] = matrix[j][j] + 1
					matrix[i][i] = matrix[i][i] + 1
		
		#print matrix
		u,w,vt = numpy.linalg.svd(numpy.matrix(matrix))
		return w,vt
		
	def close_contact_gnm(self,cutoff=0.45):
		backbone_particles = self.groups["B*_CA"]
		nresidues = len(backbone_particles)
		natoms = len(self.protein_positions)
		matrix = numpy.zeros((nresidues,nresidues),"float")
		'''
		for i in range(nresidues):
			for j in range(i+1,nresidues):
				distance =CartesianToolkit.veclength(CartesianToolkit.vecsub(self.backbone_positions[i],self.backbone_positions[j]))
				#print i,j, distance, cutoff
				if distance < cutoff:
					#print "contact!"
					matrix[i][j] = matrix[i][j] - 1
					matrix[j][i] = matrix[j][i] - 1
					matrix[j][j] = matrix[j][j] + 1
					matrix[i][i] = matrix[i][i] + 1
		'''
		#generate a residue hash
		residue_lookup = {}
		previous = 0
		for particle_id in range(1,len(backbone_particles)):
			for i in range(previous,backbone_particles[particle_id]):
				residue_lookup[i] = particle_id -1
			previous = backbone_particles[particle_id]
		for i in range(backbone_particles[-1],natoms):
			residue_lookup[i] = nresidues -1
		
		for i in range(natoms):
			
			#What is the current residue?
			a_residue = residue_lookup[i]
			
			for j in range(i+1,natoms):
				if CartesianToolkit.veclength(CartesianToolkit.vecsub(self.protein_positions[i],self.protein_positions[j])) < cutoff:
					b_residue = residue_lookup[j]

					matrix[a_residue][b_residue] -= 1.
					matrix[b_residue][a_residue] -= 1.
					matrix[a_residue][a_residue] += 1.
					matrix[b_residue][b_residue] += 1.
		
		#print matrix
		u,w,vt = numpy.linalg.svd(numpy.matrix(matrix))
		return w,vt
		
	def contact_analysis(self):
		#useful stuff
		backbone_particles = self.groups["B*_CA"]
		nresidues = len(backbone_particles)
		natoms = len(self.protein_positions)
		#generate a residue list: for each residue create a list of the atoms (counting from zero) in it
		residue_lookup = [ range(backbone_particles[resid],backbone_particles[resid+1]) for resid in range(nresidues-1)  ]
		residue_lookup.append(range(backbone_particles[-1],natoms))
		
		#what am I trying to generate here? what is the best way of representing this data?
		#Alternative ways to represent data:
		#a matrix for each time point
		#a pair of contact distances varying across time (ie a vs b vs time). this could get complicated as more species appear (ie distances between A, and B and C and D...)
		
		#For now, lets use the matrices annotated in a pseudo-xfarbe format: first line is the size/shape of the matrix
		
		distances = [ [min([CartesianToolkit.vecdist(self.protein_positions[Aatom],self.protein_positions[Batom]) for Aatom in residue for Batom in comparison_residue])  for residue in residue_lookup ] for comparison_residue  in residue_lookup ]		
		#distances is a complete matrix. It is not the most readable form. 
		
		return distances
		
	def final_analysis(self):
		pass
	def plot_graphs(self):
		pass
	def process_trajectory(self,analytic_trajfile="xy_fit.xtc",analytic_ndxfile="system.ndx"):
			length = len(self.system.initial_sequence)
	
			if length %2 == 0:
			        middle_residues = str(length/2 - 1) + "-" + str(length/2 +2)
			else:
			        middle_residues = str(length/2 - 1) + "-" + str(length/2 +1)
			
			group_mod = self.system.number_lipid_types - 1
			if self.system.charge != 0: group_mod += 1
			
			base = 15 + group_mod
			make_ndx_grps = str(base) + " & " + str(base + 1) + "\n " + str(base) + " & " + str(base + 2) + "\n " + str(base) + " & " + str(base + 3) + "\n " +str(base) + " & " + str(base + 4) + "\n "
			make_ndx_command = 'echo "a b* CA\nri 1-7\n ri ' + str(length-6) + '-' + str(length) +'\n ri ' + middle_residues +'\n ri 1-' + str(length) + '\n ' +make_ndx_grps+ 'q\n" | '  + gromacs + 'make_ndx -f t_0 -o system.ndx'
			g_bundle_grps = str(base + 5) + "\n" + str(base + 6) + "\n" + str(base + 7) + "\n"
			g_bundle_command = "echo '" + g_bundle_grps + "'| " + gromacs + "g_bundle -f t_0 -s t_0 -na 1 -z -ok -n system.ndx "
		
			
			os.system(make_ndx_command)
			#print make_ndx_command
			os.system(g_bundle_command)
		
			#center the trajectory to get all residues in the box, and the helix in the middle of the box
			trjconv_center_command = "echo '" + str(23+group_mod) + "\n0\n' | " + gromacs + "trjconv -f t_0.xtc -s t_0.tpr -pbc mol -center -o center.xtc -n system.ndx"
			os.system(trjconv_center_command)
			if self.system.bias:
				#Nasty hack due to tempremental problem for cvs trjconv
				os.system("touch electroneg.dat; touch elements.dat")
				trjconv_xy_command = "echo '" + str(23+group_mod) + "\n0\n' | " + gromacs + "trjconv -f center.xtc -s em.gro -o xy_fit.xtc -fit rotxy+transxy -n system.ndx"
				os.system(trjconv_xy_command)
			else:
				#Don't fit to xy if unbiased- there is no way of doing this at present
				os.rename("center.xtc","xy_fit.xtc")
				
			os.rename("xy_fit.xtc",analytic_trajfile)
			os.rename("system.ndx",analytic_ndxfile)
			self.analytic_trajfile = analytic_trajfile
			self.analytic_ndxfile = analytic_ndxfile

			
	def analysis_core(self):
		self.trajectory.convert_coordinates()

		#Calculate atom positions of interest
		try:
			backbone_positions = self.index_cartesians(self.groups["B*"],self.trajectory.cartesian)
		except:
			backbone_positions = self.index_cartesians(self.groups["B*_CA"],self.trajectory.cartesian)
		self.lipid_positions = []
		for lipid in allowed_lipids:
			try:
				self.lipid_positions += self.index_cartesians(self.groups[lipid],self.trajectory.cartesian)
				#Need an additional list of lipid groups to calculate the lipid list
			except:
				pass
		
		self.backbone_positions = backbone_positions
		
		helix_A_backbone = backbone_positions[:len(self.system.initial_sequence_A)]
		helix_B_backbone = backbone_positions[len(self.system.initial_sequence_A):]
		
		self.water_positions = self.index_cartesians(self.groups["W"],self.trajectory.cartesian)
		top_comz_A = helix_A_backbone[0][2]
		top_comz_B = helix_B_backbone[0][2]
		
		self.protein_positions = self.index_cartesians(self.groups["Protein"],self.trajectory.cartesian)
		
		#Calculate centers of interest
		backbone_A_com = CartesianToolkit.center(helix_A_backbone)
		backbone_B_com = CartesianToolkit.center(helix_B_backbone)
		lipid_com = CartesianToolkit.center(self.lipid_positions)
		#(bilayer_midpoint,upper_leaflet,lower_leaflet) = bilayer_center(lipid_list_breakdown(lipid_positions,lipid_atoms["DPPC"]),lipid_com)
		(bilayer_midpoint,upper_leaflet,lower_leaflet) = self.bilayer_center(self.lipid_list_breakdown_mix(self.lipid_map,self.lipid_positions),lipid_com)
		
		if top_comz_A < bilayer_midpoint:
			[upper_leaflet,lower_leaflet] = [lower_leaflet,upper_leaflet]
			
		#Now onto the core analysis: first position, tilt, screw rotation for each helix
		
		#Bilayer position
		bkbone_bilayer_distanceA = self.bilayer_position(backbone_A_com,bilayer_midpoint,top_comz_A)
		bkbone_bilayer_distanceB = self.bilayer_position(backbone_B_com,bilayer_midpoint,top_comz_B)
		
		with open(self.system.initial_sequence_A + ".0." + self.bilayer_position_file,"a") as output:
			print >> output, self.frame, 10*bkbone_bilayer_distanceA
		with open(self.system.initial_sequence_B + ".1." + self.bilayer_position_file,"a") as output:
			print >> output, self.frame, 10*bkbone_bilayer_distanceB
		
		#Rotation
		##Uses the new algorithm avoiding manual rotations
		
		(internal_tiltA, screw_angleA) = self.helix_rotation(helix_A_backbone,top_comz_A,bilayer_midpoint)
		(internal_tiltB, screw_angleB) = self.helix_rotation(helix_B_backbone,top_comz_B,bilayer_midpoint)
		
		with open(self.system.initial_sequence_A + ".0." + self.rotation_file,"a") as output:
			print >> output, self.frame, 180*screw_angleA/math.pi
		with open(self.system.initial_sequence_A + ".0." + self.tilt_file,"a") as output:
			print >> output, self.frame, 180*internal_tiltA/math.pi

		with open(self.system.initial_sequence_B + ".1." + self.rotation_file,"a") as output:
			print >> output, self.frame, 180*screw_angleB/math.pi
		with open(self.system.initial_sequence_B + ".1." + self.tilt_file,"a") as output:
			print >> output, self.frame, 180*internal_tiltB/math.pi
			
		if ( bkbone_bilayer_distanceA < 1 and bkbone_bilayer_distanceA > -1 ) and ( internal_tiltA < (5*math.pi/12) or internal_tiltA > (7*math.pi/12)):
			self.helix_insertion_A += 1
		if ( bkbone_bilayer_distanceB < 1 and bkbone_bilayer_distanceB > -1 ) and ( internal_tiltB < (5*math.pi/12) or internal_tiltB > (7*math.pi/12)):
			self.helix_insertion_B += 1
		
		#Now dimer specific analysis
		
		#Distance
		
		com_distance = CartesianToolkit.veclength(CartesianToolkit.vecsub(backbone_A_com,backbone_B_com))
		with open(self.distance_file,"a") as output:
			print >> output, self.frame, com_distance*10
		
		#Handedness; limited to when distance is below a cutoff
		
		if com_distance < 1.5:
			handed_angle = self.test_handedness(helix_A_backbone,helix_B_backbone,backbone_A_com,backbone_B_com)
			with open(self.handedness_file,"a") as output:
				print >> output, self.frame, handed_angle*180/math.pi
		
		#GNM eigenvalue analysis
		try:
			(eigenvalues,eigenvectors) = self.close_contact_gnm(cutoff=0.7)
			#for now, record the second lowest eigenvalue. Also, as a sanity check the lowest, which should always be zero
			with open(self.cc_eigenvalue_file,"a") as output:
				print >> output, self.frame, round(eigenvalues[-2],8), round(eigenvalues[-1],8)
		except:
			pass
		
		try:	
			(eigenvalues,eigenvectors) = self.gnm(cutoff=1.05)
			#for now, record the second lowest eigenvalue
			with open(self.ca_eigenvalue_file,"a") as output:
				print >> output, self.frame, round(eigenvalues[-2],8), round(eigenvalues[-1],8)
		except:
			pass
		#Orientational variation of B around A
		
		A_2_B_xy = CartesianToolkit.vecsub(backbone_B_com,backbone_A_com)[:2]
		
		with open(self.vector_relation_file,"a") as output:
			print >> output, self.frame, A_2_B_xy[0], A_2_B_xy[1]
		
		#Nearest neighbour contact
		
		distances = self.contact_analysis()
		with MatrixIO.append(self.contact_file) as output:
			output.store(distances,self.frame)
		
		#C-to-N distance (for determining feasible loop distances)
		
		with open(self.tip_distance_file,"a") as output:
			print >> output, self.frame, CartesianToolkit.vecdist(helix_A_backbone[-1],helix_B_backbone[0]), CartesianToolkit.vecdist(helix_A_backbone[0],helix_B_backbone[-1])
		
		self.frame += 1

class AnalseHelixOligomer(AnalyseHelixDimer):
	def analysis_core(self):
		self.trajectory.convert_coordinates()

		#Calculate atom positions of interest
		try:
			backbone_positions = self.index_cartesians(self.groups["B*"],self.trajectory.cartesian)
		except:
			backbone_positions = self.index_cartesians(self.groups["B*_CA"],self.trajectory.cartesian)
		self.lipid_positions = []
		for lipid in allowed_lipids:
			try:
				self.lipid_positions += self.index_cartesians(self.groups[lipid],self.trajectory.cartesian)
				#Need an additional list of lipid groups to calculate the lipid list
			except:
				pass
		
		self.backbone_positions = backbone_positions
		
		helix_backbones = [] #[backbone_positions for item in self.system.sequences]
		curr_res = 0
		for item in self.system.sequences:
			helix_backbones.append(backbone_positions[curr_res:curr_res+len(item)])
			curr_res += len(item)
		#helix_A_backbone = backbone_positions[:len(self.system.initial_sequence_A)]
		#helix_B_backbone = backbone_positions[len(self.system.initial_sequence_A):]
		
		self.water_positions = self.index_cartesians(self.groups["W"],self.trajectory.cartesian)
		top_comz = [item[2] for item in helix_backbones]
		#top_comz_A = helix_A_backbone[0][2]
		#top_comz_B = helix_B_backbone[0][2]
		
		self.protein_positions = self.index_cartesians(self.groups["Protein"],self.trajectory.cartesian)
		
		#Calculate centers of interest
		backbone_com = [CartesianToolkit.center(item) for item in helix_backbones]
		#backbone_A_com = CartesianToolkit.center(helix_A_backbone)
		#backbone_B_com = CartesianToolkit.center(helix_B_backbone)
		lipid_com = CartesianToolkit.center(self.lipid_positions)
		#(bilayer_midpoint,upper_leaflet,lower_leaflet) = bilayer_center(lipid_list_breakdown(lipid_positions,lipid_atoms["DPPC"]),lipid_com)
		(bilayer_midpoint,upper_leaflet,lower_leaflet) = self.bilayer_center(self.lipid_list_breakdown_mix(self.lipid_map,self.lipid_positions),lipid_com)
		
		if top_comz[0] < bilayer_midpoint:
			[upper_leaflet,lower_leaflet] = [lower_leaflet,upper_leaflet]
			
		#Now onto the core analysis: first position, tilt, screw rotation for each helix
		
		#Bilayer position
		bkbone_bilayer_distance = [self.bilayer_position(backbone_a,bilayer_midpoint,top_comz_a) for (backbone_a,top_comz_a) in zip(backbone_com,top_comz)]
		#bkbone_bilayer_distanceA = self.bilayer_position(backbone_A_com,bilayer_midpoint,top_comz_A)
		#bkbone_bilayer_distanceB = self.bilayer_position(backbone_B_com,bilayer_midpoint,top_comz_B)
		
		for (item,dist,id) in zip(self.system.sequences,bkbone_bilayer_distance,range(len(self.system.sequences))):
			with open(item + "." +str(id) +"." + self.bilayer_position_file,"a") as output:
				print >> output, self.frame, 10*dist

		
		#Rotation
		##Uses the new algorithm avoiding manual rotations
		
		rotation = [ self.helix_rotation(backbone_a,top_comz_a,bilayer_midpoint) for (backbone_a,top_comz_a) in zip(backbone_com,top_comz) ]
		
		#(internal_tiltA, screw_angleA) = self.helix_rotation(helix_A_backbone,top_comz_A,bilayer_midpoint)
		#(internal_tiltB, screw_angleB) = self.helix_rotation(helix_B_backbone,top_comz_B,bilayer_midpoint)
		
		for (item,tiltangles,id) in zip(self.system.sequences,rotation,range(len(self.system.sequences))):
			with open(item + "." + str(id) + "." + self.rotation_file,"a") as output:
				print >> output, self.frame, 180*tiltangles[1]/math.pi
			with open(item + "." + str(id) + "." + self.tilt_file,"a") as output:
				print >> output, self.frame, 180*tiltangles[0]/math.pi
		
		#Need to think of something better than a per helix efficiency....
			
		#if ( bkbone_bilayer_distanceA < 1 and bkbone_bilayer_distanceA > -1 ) and ( internal_tiltA < (5*math.pi/12) or internal_tiltA > (7*math.pi/12)):
		#	self.helix_insertion_A += 1
		#if ( bkbone_bilayer_distanceB < 1 and bkbone_bilayer_distanceB > -1 ) and ( internal_tiltB < (5*math.pi/12) or internal_tiltB > (7*math.pi/12)):
		#	self.helix_insertion_B += 1
		
		#Now oligomer specific analysis
		
		#Distance- calculate a pairwise matrix of interhelix backbone com distances
		
		distance_matrix  = [ [ CartesianToolkit.veclength(CartesianToolkit.vecsub(item,jitem))*10 for jitem in backbone_com] for item in backbone_com]
		
		with MatrixIO.append(self.distance_file) as output:
			output.store(distance_matrix,self.frame)
		
		#Handedness matrix; unlike the dimer analysis, this doesn't involve a distance cutoff... Yet.
		
		handedness_matrix = [ [ self.test_handedness(item,jitem,CartesianToolkit.center(item),CartesianToolkit.center(jitem))*180/math.pi for jitem in backbone_com] for item in backbone_com]
		
		with MatrixIO.append(self.handedness_file) as output:
			output.store(handedness_matrix,self.frame)
		
		#GNM eigenvalue analysis
		try:
			(eigenvalues,eigenvectors) = self.close_contact_gnm(cutoff=0.7)
			#for now, record the second lowest eigenvalue. Also, as a sanity check the lowest, which should always be zero
			with open(self.cc_eigenvalue_file,"a") as output:
				print >> output, self.frame, round(eigenvalues[-2],8), round(eigenvalues[-1],8)
		except:
			pass
		
		try:	
			(eigenvalues,eigenvectors) = self.gnm(cutoff=1.05)
			#for now, record the second lowest eigenvalue
			with open(self.ca_eigenvalue_file,"a") as output:
				print >> output, self.frame, round(eigenvalues[-2],8), round(eigenvalues[-1],8)
		except:
			pass
		#Orientational variation of B around A
		
		vector_matrix = [ [CartesianToolkit.vecsub(jitem,item)[:2] for jitem in backbone_com] for item in backbone_com ]
		
		with MatrixIO.append_3D(self.vector_relation_file) as output:
			output.store(vector_matrix,self.frame)
		
		#Nearest neighbour contact
		
		distances = self.contact_analysis()
		with MatrixIO.append(self.contact_file) as output:
			output.store(distances,self.frame)
		
		#C-to-N distance (for determining feasible loop distances), only used/sensible in alternating topologies
		
		tip_matrix = [ [ CartesianToolkit.vecdist(helix_backbones[item][-1],helix_backbones[jitem][0]) for jitem in range(len(self.system.sequences))] for item in range(len(self.system.sequences))]
		
		with MatrixIO.append(self.tip_distance_file):
			output.store(tip_matrix,self.frame)
		#with open(self.tip_distance_file,"a") as output:
		#	print >> output, self.frame, CartesianToolkit.vecdist(helix_A_backbone[-1],helix_B_backbone[0]), CartesianToolkit.vecdist(helix_A_backbone[0],helix_B_backbone[-1])
		
		self.frame += 1
		
