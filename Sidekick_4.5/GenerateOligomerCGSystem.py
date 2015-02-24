#!/usr/bin/python

from __future__ import with_statement

import sys,os,pdbio

from HAConf import configuration,programs,debug_status
from LipidBox import add_solvent_dppc, add_solvent_dlpc, add_solvent_dopc, add_solvent_popc
from MDPGenerator import replace_seed, replace_steps
import AutomatedPlot, PeptideBuild, GenerateCGSystem, GenerateDimerCGSystem

gromacs = configuration['gromacs_location']
scripts = configuration['misc_scripts']
#Lookup table for sizes
defined_system_sizes = {
			"S"		:	{
				"xy":		10,
				"lipid":	4,
				"xyz":		11.5
			}
		}

tail_lookup = 	{	"DP" : add_solvent_dppc,
					"DL" : add_solvent_dlpc,
					"PO" : add_solvent_popc,
					"DO" : add_solvent_dopc	
				}

allowed_lipids = [ "DPPC", "DPPG", "DPPE", "DLPC", "DLPG", "DLPE", "DOPC", "DOPG", "DOPE", "POPC", "POPG", "POPE" ]


class generic_oligomer_system(GenerateDimerCGSystem.generic_dimer_system):
	def __init__(self,sequences,lipid_type="DPPC",topology_name="ioned_topol",pdb_name="ioned.pdb",bias=True,seed=None,angle=0,position=0.,system_size="XS",preformed_insertion=True,special=None,rotatebyangle=0,alternating=False):
		self.sequences = sequences #expecting a list
		self.lipid_type="DPPC" #No options at present
		self.pdb_name = pdb_name
		self.topology_name = topology_name
		self.charge = 0
		self.epsilon = 20
		self.bias = bias
		self.seed = seed
		self.angle = angle % 360
		self.rotatebyangle = rotatebyangle % 360
		self.position = position
		#self.system_size = defined_system_sizes[system_size]
		self.system_xy = 9 # No options yet
		self.preformed_insertion = True #No options at present
		self.special = special
		self.apolar = True
		self.pme = False
		
		self.number_lipid_types = 1
		
		self.alternating = alternating
		
		#First generate n helices and their itps
		
		#Where to put everything? Best known packings of equal circles in a square from http://hydra.nat.uni-magdeburg.de/packing/csq/csq.html
		#All are based on a 1AU square, with an origin in the center
		#If you need more, you need to add them explicitly
		
		self.best_packings ={	3	:	{	"cart"	:	[[-0.246,  -0.246],[0.246,  -0.114],[-0.114,   0.246]],
											"radi"	:	0.254
										},
								4	:	{	"cart"	:	[[-0.25,-0.25],[0.25,-0.25],[-0.25,0.25],[0.25,0.25]],
											"radi"	:	0.25
										},
								5	:	{	"cart"	:	[[-0.293,  -0.293],[0.293,  -0.293],[0,0],[-0.293,  0.293],[0.293,  0.293]],
											"radi"	:	0.207
										},
								6	:	{	"cart"	:	[[-0.312, -0.312],[0.312, -0.312],[0, -0.104],[-0.312, 0.104],[0.312, 0.104],[0, 0.312]],
											"radi"	:	0.188
										},
								7	:	{	"cart"	:	[[-0.326, -0.326],[0.023, -0.326],[0.326, -0.151],[-0.326, 0.023],[0.023, 0.023],[0.3,0.3],[-0.151, 0.326]],
											"radi"	:	0.174
										},
								8	:	{	"cart"	:	[[-0.330,-0.330],[0.330,-0.330],[0,-0.241],[-0.241,0],[0.241,0],[0,0.241],[-0.330,0.330],[0.330,0.330]],
											"radi"	:	0.170
										},
								9	:	{	"cart"	:	[[-0.333,-0.333],[0,-0.333],[0.333,-0.333],[-0.333,0],[0,0],[0.333,0],[-0.333,0.333],[0,0.333],[0.333,0.333]],
											"radi"	:	0.167
									}
							}
		
		if len(self.sequences) not in self.best_packings.keys():
			print "No best packing has been coded for", len(self.sequences), "helices. Exiting"
			sys.exit()
		
		best_packing = self.best_packings[len(self.sequences)]
		self.topology_header_generator()
		self.radius = best_packing["radi"] * self.system_xy
		if self.radius < 1.2:
			print "WARNING: shortest distance between helices is shorter than the cut off"
		else:
			print "Shortest distance between helices is", self.radius

		for item in range(len(self.sequences)):
			self.initial_sequence = sequences[item]
			self.generate_atomistic_helix()
			self.CoarseGrain()
			rotate = 180 * (item%2) if self.alternating else 0
			self.PeptideEM(xrotate=rotate,zrotate=self.rotatebyangle,x_displacement=best_packing["cart"][item][0]*self.system_xy,y_displacement=best_packing["cart"][item][1]*self.system_xy)
			self.RenameCGFiles("Helix"+str(item))
		
		self.initial_sequence = self.sequences[0] #Initial sequence is the first sequence (A)
		
		#Combine the peptides into a single PDB and prepare for a dimer
		
		os.system("cat Helix*.pdb | grep ATOM >> prot_em.pdb")
		self.oligomer_topology_header_generator()
		with open("preion_topol.top","w") as preion_topology:
			print >> preion_topology, self.topology_header
		
		#Now reposition them, insert into a preformed bilayer/water system and add counter ions
		
		self.BuildMe()

	def PeptideEM(self,xrotate=0,zrotate=0,x_displacement=0,y_displacement=0):
		topology = open("preion_topol.top","w")
		topology.write(self.topology_header)
		topology.flush()
		
		os.system(gromacs+"grompp -f "+ configuration['mdp_files'] + "cg-em_extended.mdp -c protein_cg_box.pdb -o prot_em -p preion_topol >& grompp_protem.log")
		os.system(gromacs+"mdrun -deffnm prot_em")
		
		#os.system(gromacs+"genbox -cp prot_em.gro -cs "+ configuration['martini_systems'] + "solvent.pdb -o deionised_system.pdb -vdwd 0.24 >& genbox.log")
		os.system(gromacs+"editconf -f prot_em.gro -o prot_em.gro -center 0 0 0")
		os.system(gromacs+"editconf -f prot_em.gro -o prot_em_post_rot.pdb -rotate " + str(xrotate) + " 0 " + str(zrotate))		
		os.system(gromacs+"editconf -f prot_em_post_rot.pdb -o prot_em.pdb -translate " + str(x_displacement) + " " + str(y_displacement) + " " + str(self.position))
		topology.close()

class generic_BondOligomer_System(generic_oligomer_system):
	def oligomer_topology_header_generator(self):
		self.topology_header = """#include "ff_v1.4_polar-backbone.itp"
#include "bond-dppg.itp"
"""
		for item in range(len(self.sequences)):
			self.topology_header += '#include "Helix'+ str(item) +'.itp"\n'
		self.topology_header += """#include "special_ion.itp"
[ system ]
; Name
TM Helix Dimer
		
[ molecules ]
; Compound        #mols
"""
		for item in range(len(self.sequences)):
			self.topology_header += 'Helix'+ str(item) +'	1\n'
	
class generic_MARTINIOligomer_System(generic_oligomer_system):
	def oligomer_topology_header_generator(self):
		self.topology_header = self.martini_base() + """#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_ions.itp"
#include "martini-dppg-tieleman.itp"
"""
		for item in range(len(self.sequences)):
			self.topology_header += '#include "Helix'+ str(item) +'.itp"\n'
		self.topology_header += """
[ system ]
; Name
TM Helix Dimer
		
[ molecules ]
; Compound        #mols
"""
		for item in range(len(self.sequences)):
			self.topology_header += 'Helix'+ str(item) +'	1\n'

class MARTINIOligomer_System(generic_MARTINIOligomer_System,GenerateCGSystem.MARTINI_system):
	pass

class BondOligomer_System(generic_BondOligomer_System,GenerateCGSystem.BOND_system):
	pass

class BondOligomer_system_0_9_5(generic_BondOligomer_System,GenerateCGSystem.BOND_system_0_9_5):
	pass

class BondOligomer_system_development(generic_BondOligomer_System,GenerateCGSystem.BOND_system_development):
	pass
	
class BondOligomer_uncharged_DE(generic_BondOligomer_System,GenerateCGSystem.BOND_uncharged_DE):
	pass

class MARTINIOligomer_system_1_1_1(generic_MARTINIOligomer_System,GenerateCGSystem.MARTINI_system_1_1_1):
	pass
		
class MARTINIOligomer_system_1_1_2(generic_MARTINIOligomer_System,GenerateCGSystem.MARTINI_system_1_1_2):
	pass

class MARTINIOligomer_system_1_1_2_b(generic_MARTINIOligomer_System,GenerateCGSystem.MARTINI_system_1_1_2_b):
	pass

class MARTINIOligomer_system_112b_ENM(generic_MARTINIOligomer_System,GenerateCGSystem.MARTINI_system_112b_ENM):
	pass


installed_models = {	"MARTINI_1.1.1"			:	MARTINIOligomer_system_1_1_1,
						"MARTINI_1.1.2"			:	MARTINIOligomer_system_1_1_2,
						"MARTINI_1.1.2.b"		:	MARTINIOligomer_system_1_1_2_b,
						"Bond"					:	BondOligomer_System,
						"Bond0.9.5"				:	BondOligomer_system_0_9_5,
						"MARTINI"				:	MARTINIOligomer_System,
						"BondDev"				:	BondOligomer_system_development,
						"BondUnchargedDE"		:	BondOligomer_uncharged_DE,
						"MARTINI_1.1.2.b.ENM"	:	MARTINIOligomer_system_112b_ENM,
						"Latest-Stable-MARTINI"	:	MARTINIOligomer_system_1_1_2_b,
						"Latest-Stable-Bond"	:	BondOligomer_System
					}