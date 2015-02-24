#!/usr/bin/python

import sys,os,pdbio

from HAConf import configuration,programs,debug_status
from LipidBox import add_solvent_dppc, add_solvent_dlpc, add_solvent_dopc, add_solvent_popc
from MDPGenerator import replace_seed, replace_steps
import AutomatedPlot, PeptideBuild, GromacsInterface

gromacs = configuration['gromacs_location']
#pymol = programs['pymol']
#vmd = programs['vmd']
scripts = configuration['misc_scripts']
#Lookup table for sizes
defined_system_sizes = {
			"XS"	:	{
				"xy":		7,
				"lipid":	2,
				"xyz":		9
			},
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

class generic_system:
	'''DPPC lipids, neutralised. No method specified'''
	def __init__(self,sequence,lipid_type="DPPC",topology_name="ioned_topol",pdb_name="ioned.pdb",bias=True,seed=None,angle=0,position=0.,system_size="XS",preformed_insertion=False,special=None):
		self.initial_sequence=sequence
		self.lipid_type=lipid_type
		self.pdb_name = pdb_name
		self.topology_name = topology_name
		self.charge = 0
		self.epsilon = 20
		self.bias = bias
		self.seed = seed
		self.angle = angle
		self.position = position
		self.system_size = defined_system_sizes[system_size]
		self.preformed_insertion = preformed_insertion
		self.special = special
		self.apolar = True
		self.pme = False
		
		self.topology_header_generator()
		#
		#First generate a helix and an itp
		self.generate_atomistic_helix()
		self.CoarseGrain()
		self.PeptideEM()
		self.BuildMe()
		
	def generate_atomistic_helix(self):
		'''
		run_command = pymol + """ -cq """ + scripts + "BuildHelix.py -d """
		run_command = run_command + """'BuildMe("""
		run_command = run_command + '"' + self.initial_sequence + '"' + """)'"""
		os.system(run_command)
		'''
		PeptideBuild.BuildHelix(self.initial_sequence)
	def fill_solvent_box(self,topology):
		'''
		process_lipid_mix
		add_dp_lipids
		add_dl_lipids
		write_topology
		'''

		if len(self.lipid_type.split("/")) == 1:
			self.number_lipid_types = 1
		else:
			self.number_lipid_types = len(self.lipid_type.split("/")[0].split(":"))
			for item in self.lipid_type.split("/")[1].split(":"):
				if float(item) == 0:
					self.number_lipid_types -= 1

		#What is our lipid tail? Only allow one type of lipid tail in the mix
		lipid_tail_mix = [item[:2] for item in self.lipid_type.split("/")[0].split(":")]
		if lipid_tail_mix[0] in tail_lookup.keys():
			add_lipids = tail_lookup[lipid_tail_mix[0]]
		#Check if all lipid tails are the same. If not, you've got to correct it to run
		for item in lipid_tail_mix:
			if item != lipid_tail_mix[0]:
				print "Only a single lipid tail type allowed (ie bilayers must be pure DP, DL, PO etc.)"
				sys.exit()
		
		if self.seed:
			if self.bias:
				add_lipids(self.system_size["xy"],15,9,self.system_size["lipid"], "prot_em.pdb", "deionised_system.pdb",seed=self.seed,position=self.position,apolar=self.apolar)
			else:
				add_lipids(self.system_size["xyz"],self.system_size["xyz"],self.system_size["xyz"], self.system_size["lipid"],"prot_em.pdb","deionised_system.pdb",seed=self.seed,position=self.position,apolar=self.apolar)
		else:
			if self.bias:
				add_lipids(self.system_size["xy"],15,9,self.system_size["lipid"],"prot_em.pdb","deionised_system.pdb",position=self.position,apolar=self.apolar)
			else:
				add_lipids(self.system_size["xyz"],self.system_size["xyz"],self.system_size["xyz"],self.system_size["lipid"],"prot_em.pdb","deionised_system.pdb",position=self.position,apolar=self.apolar)
			
		lipids = 0
		water = 0
		
		for line in open("deionised_system.pdb","r"):
			if line[17:20] == "DPP" and line[13:16] == "NC3": lipids += 1
			if line[17:20] == "DLP" and line[13:16] == "NC3": lipids += 1
			if line[17:20] == "DOP" and line[13:16] == "NC3": lipids += 1
			if line[17:20] == "POP" and line[13:16] == "NC3": lipids += 1
			if line[17:20] == "  W" : water += 1
			elif line[17:20] == " PW" and line[13:15] == "W " : water += 1
		
		if self.apolar:
			water_name = "W"
		else:
			water_name = "PW"
		
		if self.lipid_type == "DPPC":
			#old style behaviour
			topology.write("DPPC        " + str(lipids) + "\n" + water_name + "         " + str(water))
			lipid_topology = "DPPC        " + str(lipids) + "\n"
		else:
			#process line
			lipid_mix = self.lipid_type.split("/")
			if len(lipid_mix) == 1:
				#Pure bilayer- check component and write as appropriate
				if lipid_mix[0] in allowed_lipids:
					lipid_string = lipid_mix[0] + "        "
					topology.write(lipid_string + str(lipids) + "\n" + water_name + "         " + str(water))
					#lipid_topology = lipid_string + str(lipids) + "\nW         " + str(water)
					lipid_topology = lipid_string + str(lipids) + "\n"
				else:
					print "Lipid not allowed! Exiting"
					sys.exit()
			else:
				#Mixed bilayer- first check components
				lipid_contents = lipid_mix[0].split(":")
				for item in lipid_contents:
					if item not in allowed_lipids:
						print "Lipid", item, "not allowed! Exiting"
						sys.exit()
				lipid_ratio = lipid_mix[1].split(":")
				lipid_ratio = [float(item) for item in lipid_ratio]
				lipid_total = sum(lipid_ratio)
				#Reminder- lipids is the actual number of lipids
				remaining_lipids = lipids
				lipid_topology = ""
				for quantity,lipid_choice in zip(lipid_ratio, lipid_contents):
					current_number = int(round(quantity * lipids / lipid_total))
					if current_number > remaining_lipids: current_number = remaining_lipids
					lipid_string = lipid_choice + "        " + str(current_number) + "\n"
					topology.write(lipid_string )
					lipid_topology += lipid_string
					remaining_lipids -= current_number
				if remaining_lipids > 0:
					topology.write(lipid_string + str(remaining_lipids) + "\n")
				topology.write(water_name + "         " + str(water))
				
		topology.close()
		return (lipids,water,lipid_topology)		
	def preformed_box(self,topology):
		if not self.apolar:
			print "No preformed apolar systems available"
			sys.exit()
		os.system(gromacs+"editconf -f prot_em.pdb -o positioned.pdb -box 9 9 15 -center 4.5 4.5 " + str( (7.5 + self.position))+ " >& editconf_position.log")
		os.system(gromacs+"genbox -cp positioned.pdb -cs " + configuration['pmf'] + "PMF_membrane_system.pdb -o deionised_system.pdb >& genbox.log")
			
		lipids = 0
		water = 0
		
		for line in open("deionised_system.pdb","r"):
			if line[17:20] == "DPP" and line[13:16] == "NC3": lipids += 1
			if line[17:20] == "  W" : water += 1
		
		topology.write("DPPC        " + str(lipids) + "\nW         " + str(water))
		lipid_topology = "DPPC        " + str(lipids) + "\n"
		topology.close()
		return (lipids,water,lipid_topology)
	def PeptideEM(self):
		topology = open("preion_topol.top","w")
		topology.write(self.topology_header)
		topology.flush()
		
		#os.system(gromacs+"grompp -f "+ configuration['mdp_files'] + "cg-em_extended.mdp -c protein_cg_box.pdb -o prot_em -p preion_topol >& grompp_protem.log")
		#os.system(gromacs+"mdrun -deffnm prot_em")
		GromacsInterface.grompp(options={"f":configuration['mdp_files'] + "cg-em_extended.mdp", "c":"protein_cg_box.pdb", "o":"prot_em", "p":"preion_topol", "maxwarn":1},logfile="grompp_protem.log" )
		GromacsInterface.mdrun(options={"deffnm":"prot_em"})		


		#os.system(gromacs+"genbox -cp prot_em.gro -cs "+ configuration['martini_systems'] + "solvent.pdb -o deionised_system.pdb -vdwd 0.24 >& genbox.log")
		os.system("echo '1\n' |"+gromacs+"editconf -f prot_em.gro -o prot_em_pre_rot.pdb -princ")
		os.system(gromacs+"editconf -f prot_em_pre_rot.pdb -o prot_em_post_rot.pdb -c -rotate " + str(self.angle) + " 0 0")
		os.system(gromacs+"editconf -f prot_em_post_rot.pdb -o prot_em.pdb -translate 0 0 " + str(self.position))
		
		topology.close()
		
	def BuildMe(self):
		#Start the simulation process
		
		#self.PeptideEM()
		
		topology = open("preion_topol.top","a")
		
		if self.preformed_insertion:
			(lipids,water,lipid_topology) = self.preformed_box(topology)
		else:
			(lipids,water,lipid_topology) = self.fill_solvent_box(topology)
		
		self.AddCounterIons(lipids,water,lipid_topology)
		
class MARTINI_type_system(generic_system):
	def topology_header_generator(self):
		self.topology_header = self.martini_base() + """#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_ions.itp"
#include "martini-dppg-tieleman.itp"
#include "protein-cg.itp"
		
[ system ]
; Name
TM Helix
		
[ molecules ]
; Compound        #mols
Protein           1\n"""
	def patch(self):
		return None
	def GenerateSSFile(self):
		secondary_structure = open("ssdump.ssd","w")
		if len(self.initial_sequence) >= 2:
			secondary_structure.write(str(len(self.initial_sequence))+"\n~"+"H"*(len(self.initial_sequence)-2) +"~")
		else:
			secondary_structure.write(str(len(self.initial_sequence))+"\n~")
		secondary_structure.close()
	def CoarseGrain(self):
		#Generate an itp
		sequence = open("sequence.seq","w")
		sequence.write(""">chainA\n"""+self.initial_sequence)
		sequence.close()
		
		self.GenerateSSFile()
		
		martini_script = self.script_name()
		
		os.system("ln -s " + configuration['martini_setup'] + "* .")
		cg_command = "perl " + martini_script + " -itp protein-cg.itp"
		if self.special:
			cg_command += " --" + self.special 
		os.system("echo "+cg_command)
		os.system(cg_command)
		GromacsInterface.pdb2gmx(options={"f":"tm.pdb", "p":"tm", "o":"atomistic.pdb", "ignh":True}, selections=["13","3"], logfile="pdb2gmx.log")
		#os.system("echo '5\n' |" + gromacs + "pdb2gmx -f tm.pdb -o atomistic.pdb -ignh >& pdb2gmx.log")
		os.system(gromacs+"editconf -f atomistic.pdb -o protein_chain.pdb -label A >& editconf_label.log")
		os.system(configuration['martini_setup']+"atom2cg_v2.1tryout.awk protein_chain.pdb > protein_cg.pdb")
		os.system(gromacs+"editconf -f protein_cg.pdb -o protein_cg_box.pdb -box 7 7 15 -center 3.5 3.5 7.5 >& editconf_box.log")

		#Run any script independent "patches"
		patchme = self.patch()
		if patchme != None:
			#Apply patch
			print "Applying patch '",patchme,"'"
			os.system(patchme)
		
	def AddCounterIons(self,lipids,water,lipid_topology):
		#os.system(gromacs+"grompp -f " + configuration['mdp_files'] + "cg-em_extended.mdp -c deionised_system.pdb -p preion_topol -o solvated >& ionisation_state.log")
		GromacsInterface.grompp(options={"f":configuration['mdp_files'] + "cg-em_extended.mdp", "c":"deionised_system.pdb", "o":"solvated", "p":"preion_topol", "maxwarn":1},logfile="ionisation_state.log" )
	
		charge = 0
		for line in open("ionisation_state.log"):
			if line[:36] == "  System has non-zero total charge: ": charge = float(line[36:])
		
		if charge == 0:
			os.system("cp preion_topol.top ioned_topol.top")
			os.system("cp deionised_system.pdb ioned.pdb")
		else:
			ion_topology = open("ioned_topol.top","w")
			ion_topology.write(self.topology_header)
			if charge > 0:
				charge_flag = "-nn"
				charge = int(charge)
				ion_name = "CL-"
			else:
				charge_flag = "-np"
				charge = -int(charge)
				ion_name = "NA+"
			#How many lipids? + 1 to the group for each additional lipid
			if len(self.lipid_type.split("/")) == 1:
				#gromacs 4+ has decided to move water's default index number
				os.system("echo '14\n' | " + gromacs+"genion -s solvated.tpr -o " + self.pdb_name + " " + charge_flag + " " + str(charge))
			else:
				groupid = 14 + len(self.lipid_type.split("/")[0].split(":")) - 1
				for item in self.lipid_type.split("/")[1].split(":"):
					if float(item) == 0:
						groupid -= 1
				os.system("echo '" + str(groupid) + "\n' | " + gromacs+"genion -s solvated.tpr -o " + self.pdb_name + " " + charge_flag + " " + str(charge))
			
			if self.apolar:
				water_name = "W"
			else:
				water_name = "PW"
			ion_topology.write(lipid_topology)
			ion_topology.write(water_name + "         " + str(water-charge) + "\n")
			ion_topology.write(ion_name + " "*7 + str(charge))
			ion_topology.flush()
			ion_topology.close()
		os.system("cp ioned_topol.top " + self.topology_name)
		self.charge = charge
		#replace epsilon for prefered value
		self.epsilon = 15

class MARTINI_system(MARTINI_type_system):
	def martini_base(self):
		return """#include "martini_v2.1tryout.itp"\n"""
	def script_name(self):
		#return "breakme"
		return "./seq2cgtop_martini_v2.1tryout.pl"

class MARTINI_system_1_1_1(MARTINI_type_system):
	def martini_base(self):
		return """#include "martini_v2.1.itp"\n"""
	def script_name(self):
		return "./seq2itp_martini2.1.1.pl sequence.seq ssdump.ssd"
		
class MARTINI_system_1_1_2(MARTINI_type_system):
	def martini_base(self):
		return """#include "martini_v2.1.itp"\n"""
	def script_name(self):
		return "./seq2itp_martini2.1.2.pl sequence.seq ssdump.ssd"

class MARTINI_system_1_1_2_b(MARTINI_type_system):
	def martini_base(self):
		return """#include "martini_v2.1.itp"\n"""
	def script_name(self):
		return "./seq2itp_martini2.1.2.b.pl sequence.seq ssdump.ssd"

class MARTINI_system_112b_unchargedLYS(MARTINI_type_system):
	def patch(self):
		return configuration["python"] + " ./uncharged_lys_patch.py"
	def martini_base(self):
		return """#include "martini_v2.1.itp"\n"""
	def script_name(self):
		return "./seq2itp_martini2.1.2.b.pl sequence.seq ssdump.ssd"

class MARTINI_system_112b_methylatedLYS(MARTINI_type_system):
	def patch(self):
		return configuration["python"] + " ./methylated_lys_patch.py"
	def martini_base(self):
		return """#include "martini_v2.1.itp"\n"""
	def script_name(self):
		return "./seq2itp_martini2.1.2.b.pl sequence.seq ssdump.ssd"

class MARTINI_system_112b_ENM(MARTINI_type_system):
	def patch(self):
		return configuration["python"] + " ./gnm_patch.py"
		#return None
	def martini_base(self):
		return """#include "martini_v2.1.itp"\n"""
	def script_name(self):
		return "./seq2itp_martini2.1.2.b.pl sequence.seq ssdump.ssd"

class MARTINI_system_1_1_2_b_super_helix(MARTINI_system_1_1_2_b):
	def GenerateSSFile(self):
		secondary_structure = open("ssdump.ssd","w")
		if len(self.initial_sequence) >= 2:
			secondary_structure.write(str(len(self.initial_sequence))+"\n"+"H"*(len(self.initial_sequence)) +"")
		else:
			secondary_structure.write(str(len(self.initial_sequence))+"\nH")
		secondary_structure.close()	

class PMARTINI_system(MARTINI_type_system):
	def martini_base(self):
		self.apolar = False
		self.epsilon = 2.5
		return """#include "martini_v2.P.itp"\n"""
	def script_name(self):
		return "./seq2itp_martini2.P.pl sequence.seq ssdump.ssd"

class BOND_system(generic_system):
	def topology_header_generator(self):
		self.topology_header = """#include "ff_v1.4_polar-backbone.itp"
#include "bond-dppg.itp"
#include "protein-cg.itp"
#include "special_ion.itp"
[ system ]
; Name
TM Helix
		
[ molecules ]
; Compound        #mols
Protein           1\n"""
	
	def script_name(self):
		return "./cg_partitioning_hbonds_v1.pl"
	def CoarseGrain(self):
		secondary_structure = open("structure.txt","w")
		secondary_structure.write('"~",\n'+'"H",\n'*(len(self.initial_sequence)-2) +'"~"')
		secondary_structure.close()
		
		#Generate a sequence.seq to make redrawing the graphs easier
		sequence = open("sequence.seq","w")
		sequence.write(""">chainA\n"""+self.initial_sequence)
		sequence.close()
		
		bond_script = self.script_name()
		
		#os.system("source " + gromacs + "GMXRC")
		
		#os.system("echo '5\n' |" + gromacs + "pdb2gmx -f tm.pdb -o atomistic.pdb -ignh >& pdb2gmx.log")
		GromacsInterface.pdb2gmx(options={"f":"tm.pdb", "p":"tm", "o":"atomistic.pdb", "ignh":True}, selections=["13","3"], logfile="pdb2gmx.log")
		os.system("ln -s " + configuration['bond_setup'] + "* .")
		cg_command = "perl " + bond_script + " " + gromacs 
		if self.special:
			cg_command += " --" + self.special
		os.system(cg_command)
		
		os.system(gromacs+"editconf -f protein-cg.pdb -o protein_cg_box.pdb -box 7 7 15 -center 3.5 3.5 7.5 >& editconf_box.log")
		temporary_empty_file = open("special_ion.itp","w")
		temporary_empty_file.close()
		
	def AddCounterIons(self,lipids,water,lipid_topology):
				
		#os.system(gromacs+"grompp -f " + configuration['mdp_files'] + "cg-em_extended.mdp -c deionised_system.pdb -p preion_topol -o solvated >& ionisation_state.log")
		GromacsInterface.grompp(options={"f":configuration['mdp_files'] + "cg-em_extended.mdp", "c":"deionised_system.pdb", "o":"solvated", "p":"preion_topol", "maxwarn":1},logfile="ionisation_state.log" )
		
		charge = 0
		for line in open("ionisation_state.log"):
			if line[:36] == "  System has non-zero total charge: ": charge = float(line[36:])
		
		#Charge of "Typical" ions is 0.7 to reflect shielding in Bond models
		
		if charge == 0:
			os.system("cp preion_topol.top ioned_topol.top")
			os.system("cp deionised_system.pdb ioned.pdb")
		else:
			ion_topology = open("ioned_topol.top","w")
			ion_topology.write(self.topology_header)
			if charge > 0:
				charge_flag = "-nn"
				number_of_ions = int(charge/0.7) + 1
				ion_name = "CL"
				remainder = -(charge - round((number_of_ions - 1)*0.7,1))
			else:
				charge_flag = "-np"
				number_of_ions = int(charge/-0.7) + 1
				ion_name = "NA"
				remainder = charge - round((number_of_ions - 1)*-0.7,1)
			
			special_ion = ion_name + "B"
			
			special_ion_itp = open("special_ion.itp","w")
			special_ion_itp.write(""";;;;;; ANOTHER SPECIAL HYDRATED SODIUM ION

[ moleculetype ]
; molname       nrexcl
""" + special_ion + """                 1

[ atoms ]
;id type resnr residu atom cgnr   charge
 1   Qda   1     ION      """ + ion_name + """    1    """ + str(remainder) )
 			special_ion_itp.close()
			
			if len(self.lipid_type.split("/")) == 1:
				os.system("echo '14\n' | " + gromacs+"genion -s solvated.tpr -pq 0.7 -nq -0.7 -o " + self.pdb_name + " " + charge_flag + " " + str(number_of_ions))
			else:
				groupid = 14 + len(self.lipid_type.split("/")[0].split(":")) - 1
				os.system("echo '" + str(groupid) + "\n' | " + gromacs+"genion -s solvated.tpr -pq 0.7 -nq -0.7 -o " + self.pdb_name + " " + charge_flag + " " + str(number_of_ions))
			
			ion_topology.write(lipid_topology)
			ion_topology.write("W         " + str(water-number_of_ions) + "\n")
			ion_topology.write(ion_name + " "*7 + str(number_of_ions-1) + "\n")
			ion_topology.write(special_ion + " "*7 + str(1))
			ion_topology.flush()
			ion_topology.close()
		os.system("cp ioned_topol.top " + self.topology_name)
		self.charge = charge

class BOND_system_0_9_5(BOND_system):
	def script_name(self):
		return "./make-top-v50.pl"

class BOND_system_development(BOND_system):
	def script_name(self):
		return "./cg_partitioning_hbonds_v1.9_H+.pl"

class BOND_uncharged_DE(BOND_system):
	def script_name(self):
		return "./cg_partitioning_hbonds_v1_unchargedDE.pl"

class PME_BOND(BOND_system):
	def script_name(self):
		self.pme = True
		return "./cg_partitioning_hbonds_v1.pl"

installed_models = {	"MARTINI_1.1.1"			:	MARTINI_system_1_1_1,
						"MARTINI_1.1.2"			:	MARTINI_system_1_1_2,
						"MARTINI_1.1.2.b"		:	MARTINI_system_1_1_2_b,
						"PMARTINI"				:	PMARTINI_system,
						"Bond"					:	BOND_system,
						"Bond0.9.5"				:	BOND_system_0_9_5,
						"BondPME"				:	PME_BOND,
						"MARTINI"				:	MARTINI_system,
						"BondDev"				:	BOND_system_development,
						"BondUnchargedDE"		:	BOND_uncharged_DE,
						"MARTINI_1.1.2.b.ENM"	:	MARTINI_system_112b_ENM,
						"MARTINI_1.1.2.b.SHelix":	MARTINI_system_1_1_2_b_super_helix,
						"MARTINI_1.1.2.b.unchargedLYS":	MARTINI_system_112b_unchargedLYS,
						"MARTINI_1.1.2.b.methylatedLYS": MARTINI_system_112b_methylatedLYS,
						"Latest-Stable-MARTINI"	:	MARTINI_system_1_1_2_b,
						"Latest-Stable-Bond"	:	BOND_system
					}
