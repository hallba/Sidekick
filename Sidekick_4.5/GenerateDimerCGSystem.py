#!/usr/bin/python

from __future__ import with_statement

import sys,os,pdbio

from HAConf import configuration,programs,debug_status
from LipidBox import add_solvent_dppc, add_solvent_dlpc, add_solvent_dopc, add_solvent_popc
from MDPGenerator import replace_seed, replace_steps
import AutomatedPlot, PeptideBuild, GenerateCGSystem
import GromacsInterface

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


class generic_dimer_system(GenerateCGSystem.generic_system):
	def __init__(self,sequenceA,sequenceB,lipid_type="DPPC",topology_name="ioned_topol",pdb_name="ioned.pdb",bias=True,seed=None,angle=0,position=0.,system_size="XS",preformed_insertion=True,special=None,parallel=True,rotatebyangle=0):
		self.initial_sequence_A=sequenceA
		self.initial_sequence_B=sequenceB
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
		self.preformed_insertion = True #No options at present
		self.special = special
		self.parallel = parallel
		self.apolar = True
		self.pme = False
		
		self.number_lipid_types = 1
		
		self.topology_header_generator()
		
		#First generate 2 helices and their itps
		self.initial_sequence = self.initial_sequence_A
		self.generate_atomistic_helix()
		self.CoarseGrain()
		self.PeptideEM(zrotate=self.rotatebyangle,xy_displacement=2.25)
		#if self.rotatebyangle != 0:
		#	self.RotatePeptideZ(self.rotatebyangle)
		#self.RepositionPeptide(2.25,2.25)
		self.RenameCGFiles("HelixA")
		
		self.initial_sequence = self.initial_sequence_B
		self.generate_atomistic_helix()
		self.CoarseGrain()
		#self.PeptideEM()
		#self.RepositionPeptide(-2.25,-2.25)
		if self.parallel != True:
			print "Antiparallel helices"
			#self.RotatePeptide(180)
			self.PeptideEM(xrotate=180,zrotate=0,xy_displacement=-2.25)
		else:
			self.PeptideEM(xrotate=0,zrotate=0,xy_displacement=-2.25)
		self.RenameCGFiles("HelixB")
		
		self.initial_sequence = self.initial_sequence_A #Initial sequence is the first sequence (A)
		
		#Combine the peptides into a single PDB and prepare for a dimer
		
		os.system("cat HelixA.pdb HelixB.pdb | grep ATOM >> prot_em.pdb")
		self.dimer_topology_header_generator()
		with open("preion_topol.top","w") as preion_topology:
			print >> preion_topology, self.topology_header
		
		#Now reposition them, insert into a preformed bilayer/water system and add counter ions
		
		self.BuildMe()
		
	def RenameCGFiles(self,name):
		os.rename("prot_em.pdb",name+".pdb")
		with open("protein-cg.itp","r") as input_itp:
			with open(name+".itp","w") as output_itp:
				for line in input_itp:
					if line[:7] == "PROTEIN" or line[:7] == "Protein":
						print >> output_itp, name, 1
					else:
						print >> output_itp, line,
		os.rename("protein-cg.itp",name+"_orig.itp")
	
	def RepositionPeptide(self,x,y):
		command = gromacs+"editconf -f prot_em.pdb -o prot_em.pdb -translate %(x)8.3f %(y)8.3f 0" % {"x":x,"y":y }
		os.system(command)
	def RotatePeptide(self,degrees,y=-4.5,z=-7.5):
		x = 0
		move_to_origin_command = gromacs+"editconf -f prot_em.pdb -o prot_em.pdb -translate %(x)8.3f %(y)8.3f %(z)8.3f" % {"x":x,"y":y, "z":z}
		move_back_command = gromacs+"editconf -f prot_em.pdb -o prot_em.pdb -translate %(x)8.3f %(y)8.3f %(z)8.3f" % {"x":-x,"y":-y, "z":-z }
		command = gromacs+"editconf -f prot_em.pdb -o prot_em.pdb -rotate %(x)8.3f 0 0" % {"x":degrees}
		os.system(move_to_origin_command)
		os.system(command)
		os.system(move_back_command)
	def RotatePeptideZ(self,degrees,x=-4.5,y=-4.5):
		z = 0
		move_to_origin_command = gromacs+"editconf -f prot_em.pdb -o prot_em.pdb -translate %(x)8.3f %(y)8.3f %(z)8.3f" % {"x":x,"y":y, "z":z}
		move_back_command = gromacs+"editconf -f prot_em.pdb -o prot_em.pdb -translate %(x)8.3f %(y)8.3f %(z)8.3f" % {"x":-x,"y":-y, "z":-z }
		command = gromacs+"editconf -f prot_em.pdb -o prot_em.pdb -rotate 0 0 %(z)8.3f" % {"z":degrees}
		os.system(move_to_origin_command)
		os.system(command)
		os.system(move_back_command)
	def PeptideEM(self,xrotate=0,zrotate=0,xy_displacement=0):
		topology = open("preion_topol.top","w")
		topology.write(self.topology_header)
		topology.flush()
		
		#os.system(gromacs+"grompp -f "+ configuration['mdp_files'] + "cg-em_extended.mdp -c protein_cg_box.pdb -o prot_em -p preion_topol >& grompp_protem.log")
		#os.system(gromacs+"mdrun -deffnm prot_em")
		
                GromacsInterface.grompp(options={"f":configuration['mdp_files'] + "cg-em_extended.mdp", "c":"protein_cg_box.pdb", "o":"prot_em", "p":"preion_topol", "maxwarn":1},logfile="grompp_protem.log" )
                GromacsInterface.mdrun(options={"deffnm":"prot_em"})

		#os.system(gromacs+"genbox -cp prot_em.gro -cs "+ configuration['martini_systems'] + "solvent.pdb -o deionised_system.pdb -vdwd 0.24 >& genbox.log")
		os.system(gromacs+"editconf -f prot_em.gro -o prot_em.gro -center 0 0 0")
		os.system(gromacs+"editconf -f prot_em.gro -o prot_em_post_rot.pdb -rotate " + str(xrotate) + " 0 " + str(zrotate))		
		os.system(gromacs+"editconf -f prot_em_post_rot.pdb -o prot_em.pdb -translate " + str(xy_displacement) + " " + str(xy_displacement) + " " + str(self.position))
		topology.close()

class generic_BondDimer_System(generic_dimer_system):
	def dimer_topology_header_generator(self):
		self.topology_header = """#include "ff_v1.4_polar-backbone.itp"
#include "bond-dppg.itp"
#include "HelixA.itp"
#include "HelixB.itp"
#include "special_ion.itp"
[ system ]
; Name
TM Helix Dimer
		
[ molecules ]
; Compound        #mols
HelixA           1
HelixB           1\n"""
	
class generic_MARTINIDimer_System(generic_dimer_system):
	def dimer_topology_header_generator(self):
		self.topology_header = self.martini_base() + """#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_ions.itp"
#include "martini-dppg-tieleman.itp"
#include "HelixA.itp"
#include "HelixB.itp"
		
[ system ]
; Name
TM Helix Dimer
		
[ molecules ]
; Compound        #mols
HelixA           1
HelixB           1\n"""

class MARTINIDimer_System(generic_MARTINIDimer_System,GenerateCGSystem.MARTINI_system):
	pass

class BondDimer_System(generic_BondDimer_System,GenerateCGSystem.BOND_system):
	pass

class BONDDimer_system_0_9_5(generic_BondDimer_System,GenerateCGSystem.BOND_system_0_9_5):
	pass

class BONDDimer_system_development(generic_BondDimer_System,GenerateCGSystem.BOND_system_development):
	pass
	
class BONDDimer_uncharged_DE(generic_BondDimer_System,GenerateCGSystem.BOND_uncharged_DE):
	pass

class MARTINIDimer_system_1_1_1(generic_MARTINIDimer_System,GenerateCGSystem.MARTINI_system_1_1_1):
	pass
		
class MARTINIDimer_system_1_1_2(generic_MARTINIDimer_System,GenerateCGSystem.MARTINI_system_1_1_2):
	pass

class MARTINIDimer_system_1_1_2_b(generic_MARTINIDimer_System,GenerateCGSystem.MARTINI_system_1_1_2_b):
	pass

class MARTINIDimer_system_112b_ENM(generic_MARTINIDimer_System,GenerateCGSystem.MARTINI_system_112b_ENM):
	pass


installed_models = {	"MARTINI_1.1.1"			:	MARTINIDimer_system_1_1_1,
						"MARTINI_1.1.2"			:	MARTINIDimer_system_1_1_2,
						"MARTINI_1.1.2.b"		:	MARTINIDimer_system_1_1_2_b,
						"Bond"					:	BondDimer_System,
						"Bond0.9.5"				:	BONDDimer_system_0_9_5,
						"MARTINI"				:	MARTINIDimer_System,
						"BondDev"				:	BONDDimer_system_development,
						"BondUnchargedDE"		:	BONDDimer_uncharged_DE,
						"MARTINI_1.1.2.b.ENM"	:	MARTINIDimer_system_112b_ENM,
						"Latest-Stable-MARTINI"	:	MARTINIDimer_system_1_1_2_b,
						"Latest-Stable-Bond"	:	BondDimer_System
					}
