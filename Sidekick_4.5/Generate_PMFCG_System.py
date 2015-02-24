#!/usr/bin/python

import os,sys
import HAConf
from HAConf import configuration
from GenerateCGSystem import generic_system
import pdbio,ndxio,CartesianToolkit

gromacs = configuration['gromacs_location']

def lipid_list_breakdown(lipid_positions,number_of_atoms):
	return [lipid_positions[i*number_of_atoms:(i+1)*number_of_atoms] for i in range(len(lipid_positions)/number_of_atoms)]
def index_cartesians(index_list,frame):
	return [frame[item] for item in index_list]
	
class MARTINI_system(generic_system):
	def __init__(self,sequence,lipid_type="DPPC",topology_name="ioned_topol",pdb_name="ioned.pdb",bias=True,window=0.,force=1000.):
		self.initial_sequence=sequence
		self.lipid_type=lipid_type
		self.pdb_name = pdb_name
		self.topology_name = topology_name
		self.charge = 0
		self.bias = bias
		self.window = float(window)
		self.force = float(force)
		if self.window > 70 or self.window < -70:
			print "Window is out of the default range"
			exit()
		#First generate a helix and an itp
		self.generate_atomistic_helix()
		self.CoarseGrain()
	def script_name(self):
		return "./seq2cgtop_martini_v2.1tryout.pl"
	def CoarseGrain(self):
		#Generate an itp
		sequence = open("sequence.seq","w")
		sequence.write(""">chainA\n"""+self.initial_sequence)
		sequence.close()
		
		secondary_structure = open("ssdump.ssd","w")
		if len(self.initial_sequence) >= 2:
			secondary_structure.write(str(len(self.initial_sequence))+"\n~"+"H"*(len(self.initial_sequence)-2) +"~")
		else:
			secondary_structure.write(str(len(self.initial_sequence))+"\n~")
		secondary_structure.close()
		
		martini_script = self.script_name()
		
		os.system("ln -s " + configuration['martini_setup'] + "* .")
		os.system("perl " + martini_script + " -itp protein-cg.itp")
		os.system("echo '5\n' |" + gromacs + "pdb2gmx -f tm.pdb -o atomistic.pdb -ignh >& pdb2gmx.log")
		os.system(gromacs+"editconf -f atomistic.pdb -o protein_chain.pdb -label A >& editconf_label.log")
		os.system(configuration['martini_setup']+"atom2cg_v2.1tryout.awk protein_chain.pdb > protein_cg.pdb")
		os.system(gromacs+"editconf -f protein_cg.pdb -o protein_cg_box.pdb -box 9 9 15 -center 4.5 4.5 7.5 >& editconf_box.log")

		#Start the simulation process
		
		topology_header = """#include "martini_v2.1tryout.itp"
#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_ions.itp"
#include "protein-cg.itp"

; Include Position restraint file
#ifdef POSRES_HELIX
#include "posre_helix.itp"
#endif
		
[ system ]
; Name
TM Helix
		
[ molecules ]
; Compound        #mols
Protein           1\n"""
		
		topology = open("preion_topol.top","w")
		topology.write(topology_header)
		topology.flush()
		
		os.system(gromacs+"grompp -f "+ configuration['mdp_files'] + "cg-em_extended.mdp -c protein_cg_box.pdb -o prot_em -p preion_topol >& grompp_protem.log")
		os.system(gromacs+"mdrun -deffnm prot_em")
		
		#os.system(gromacs+"genbox -cp prot_em.gro -cs "+ configuration['martini_systems'] + "solvent.pdb -o deionised_system.pdb -vdwd 0.24 >& genbox.log")
		os.system("echo '1\n' |"+gromacs+"editconf -f prot_em.gro -o prot_em.pdb -princ")
		
		#The CG helix is now in minimised and ready for introduction to the system
		##The system is the default, in the PMF directory PMF_membrane_system with 131 lipids in the upper, 125 in the lower
		##Center in Z is at 75A, so each window is an addition or subtraction to this
		##The (predetermined) range is from 70 to -70; in principle this should allow for a 50 residue helix
		##Use shorter to be on the safe side
		
		os.system(gromacs+"editconf -f prot_em.pdb -o positioned.pdb -box 9 9 15 -center 4.5 4.5 " + str( (75. + self.window) /10)+ " >& editconf_position.log")
		os.system(gromacs+"genbox -cp positioned.pdb -cs " + HAConf.configuration['pmf'] + "PMF_membrane_system.pdb -o deionised_system.pdb >& genbox.log")
			
		lipids = 0
		water = 0
		
		for line in open("deionised_system.pdb","r"):
			if line[17:20] == "DPP" and line[13:16] == "NC3": lipids += 1
			if line[17:20] == "  W" : water += 1
		
		topology.write("DPPC        " + str(lipids) + "\nW         " + str(water))
		topology.close()
		
		os.system(gromacs+"grompp -f " + configuration['mdp_files'] + "cg-em_extended.mdp -c deionised_system.pdb -p preion_topol -o solvated >& ionisation_state.log")
		
		charge = 0
		for line in open("ionisation_state.log"):
			if line[:36] == "  System has non-zero total charge: ": charge = float(line[36:])
		
		if charge == 0:
			os.system("cp preion_topol.top ioned_topol.top")
			os.system("cp deionised_system.pdb ioned.pdb")
		else:
			ion_topology = open("ioned_topol.top","w")
			ion_topology.write(topology_header)
			if charge > 0:
				charge_flag = "-nn"
				charge = int(charge)
				ion_name = "CL-"
			else:
				charge_flag = "-np"
				charge = -int(charge)
				ion_name = "NA+"
			os.system("echo '13\n' | " + gromacs+"genion -s solvated.tpr -o " + self.pdb_name + " " + charge_flag + " " + str(charge))
			ion_topology.write("DPPC        " + str(lipids) + "\nW         " + str(water-charge) + "\n")
			ion_topology.write(ion_name + " "*7 + str(charge))
			ion_topology.flush()
			ion_topology.close()
		os.system("cp ioned_topol.top " + self.topology_name)
		self.charge = charge
		
		#Generate some backbone position restraints
		
		if charge != 0:
			make_ndx_command = 'echo "a b* CA\nkeep 16\nq\n" | '  + gromacs + 'make_ndx -f ioned.pdb -o backbone.ndx'
		else:
			make_ndx_command = 'echo "a b* CA\nkeep 15\nq\n" | '  + gromacs + 'make_ndx -f ioned.pdb -o backbone.ndx'
		
		os.system(make_ndx_command)
		genpr_command = gromacs + "genpr -f ioned.pdb -n backbone.ndx -o posre_helix.itp"
		os.system(genpr_command)
		
		#Now generate an index/ pull file for the run itself
		
		make_ndx_command = 'echo "a b* CA\nq\n" | '  + gromacs + 'make_ndx -f ioned.pdb -o ionised_system.ndx'
		os.system(make_ndx_command)
		
		##Now choose a subset of lipids to calculate the center of mass
		##
		
		structure = pdbio.read_pdb("ioned.pdb")
		groups = ndxio.read_ndx("ionised_system.ndx")
		
		structure.convert_coordinates()
		lipid_positions = index_cartesians(groups["DPP"],structure.cartesian)
		lipid_com = CartesianToolkit.center(lipid_positions)
		individual_lipid_positions = lipid_list_breakdown(lipid_positions,12)
		
		leaflet_id = {	"upper"	:	[],
						"lower"	:	[] }
		
		individual_lipid_comz = [CartesianToolkit.center(lipid)[2] - lipid_com[2] for lipid in individual_lipid_positions]
		#This needs to be set to the first non protein residue
		residue = 1 + len(self.initial_sequence)
		for lipid_z in individual_lipid_comz:
			if lipid_z > 0:
				leaflet_id["upper"].append(residue)
			else:
				leaflet_id["lower"].append(residue)
			residue += 1
		
		if len(leaflet_id["upper"]) == len(leaflet_id["lower"]):
			reference_group = "DPPC"
			make_ndx_command = 'echo "a b* CA\nq\n" | '  + gromacs + 'make_ndx -f ioned.pdb -o pull_system.ndx'
		elif len(leaflet_id["upper"]) > len(leaflet_id["lower"]):
			difference = len(leaflet_id["upper"]) - len(leaflet_id["lower"])
			reference_group = "symm"
			make_ndx_group = "r"
			symmetric_bilayer = leaflet_id["upper"][difference:] + leaflet_id["lower"]
			for residue in symmetric_bilayer:
				#reference_group += "_" + str(residue)
				make_ndx_group += " " + str(residue)
			if charge != 0:
				make_ndx_command = 'echo "a b* CA\n' + make_ndx_group + '\nname 17 symm\nq\n" | '  + gromacs + 'make_ndx -f ioned.pdb -o pull_system.ndx'
			else:
				make_ndx_command = 'echo "a b* CA\n' + make_ndx_group + '\nname 16 symm\nq\n" | '  + gromacs + 'make_ndx -f ioned.pdb -o pull_system.ndx'
		else:
			difference = len(leaflet_id["lower"]) - len(leaflet_id["upper"])
			reference_group = "symm"
			make_ndx_group = "r"
			symmetric_bilayer = leaflet_id["upper"] + leaflet_id["lower"][difference:]
			for residue in symmetric_bilayer:
				#reference_group += "_" + str(residue)
				make_ndx_group += " " + str(residue)
			if charge != 0:
				make_ndx_command = 'echo "a b* CA\n' + make_ndx_group + '\nname 17 symm\nq\n" | '  + gromacs + 'make_ndx -f ioned.pdb -o pull_system.ndx'
			else:
				make_ndx_command = 'echo "a b* CA\n' + make_ndx_group + '\nname 16 symm\nq\n" | '  + gromacs + 'make_ndx -f ioned.pdb -o pull_system.ndx'
				
		os.system(make_ndx_command)
		
		umbrella = open("umbrella.ppa","w")
		print >> umbrella, """verbose = no
runtype = umbrella
group_1 = B*_CA
reference_group = """ + reference_group + """
reftype = com_t0 ; center of mass, corrected for molecules that cross the simulation box
pulldim = N N Y ; only along z-axis
k1 = """ + str(self.force) + """
pos1= 0.0 0.0 """  + str(self.window) + """ ; the first 2 values should *not* matter
"""
		umbrella.close()

class MARTINI_system_2_1_1(MARTINI_system):
	def script_name(self):
		return "./seq2itp_martini2.1.1.pl sequence.seq ssdump.ssd"
		
		
class MARTINI_system_2_1_2(MARTINI_system):
	def script_name(self):
		return "./seq2itp_martini2.1.2.pl sequence.seq ssdump.ssd"