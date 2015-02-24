#!/usr/bin/python

import sys,os,pdbio

from HAConf import configuration,programs,debug_status
from LipidBox import add_solvent_dppc
from MDPGenerator import replace_seed, replace_steps
import AutomatedPlot

import AnalyseTrajectory
import GromacsInterface

initial_sequence=sys.argv[1]
seed = sys.argv[2]

gromacs = configuration['gromacs_location']
pymol = programs['pymol']
vmd = programs['vmd']
scripts = configuration['misc_scripts']

#First generate a helix and an itp

run_command = pymol + """ -cq """ + scripts + "BuildHelix.py -d """
run_command = run_command + """'BuildMe("""
run_command = run_command + '"' + initial_sequence + '"' + """)'"""
os.system(run_command)

sequence = open("sequence.seq","w")
sequence.write(""">chainA\n"""+initial_sequence)
sequence.close()

secondary_structure = open("ssdump.ssd","w")
secondary_structure.write(str(len(initial_sequence))+"\n~"+"H"*(len(initial_sequence)-2) +"~")
secondary_structure.close()

os.system("ln -s " + configuration['martini_setup'] + "* .")
os.system("perl ./seq2cgtop_martini_v2.1tryout.pl -itp protein-cg.itp")
os.system("echo '5\n' |" + gromacs + "pdb2gmx -f tm.pdb -o atomistic.pdb -ignh >& pdb2gmx.log")
os.system(gromacs+"editconf -f atomistic.pdb -o protein_chain.pdb -label A >& editconf_label.log")
os.system(configuration['martini_setup']+"atom2cg_v2.1tryout.awk protein_chain.pdb > protein_cg.pdb")
os.system(gromacs+"editconf -f protein_cg.pdb -o protein_cg_box.pdb -box 7 7 15 -center 3.5 3.5 7.5 >& editconf_box.log")

#Start the simulation process

topology_header = """#include "martini_v2.1tryout.itp"
#include "martini_v2.0_lipids.itp"
#include "martini_v2.0_ions.itp"
#include "protein-cg.itp"

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
add_solvent_dppc(7,15,9,2,"prot_em.pdb","deionised_system.pdb")

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
	os.system("echo '13\n' | " + gromacs+"genion -s solvated.tpr -o ioned.pdb " + charge_flag + " " + str(charge))
	ion_topology.write("DPPC        " + str(lipids) + "\nW         " + str(water-charge) + "\n")
	ion_topology.write(ion_name + " "*7 + str(charge))
	ion_topology.flush()
	ion_topology.close()


GromacsInterface.SDMinimise(topology_file="ioned_topol",input_structure="ioned.pdb",output_name="em")

GromacsInterface.RobustMolecularDynamics(topology_file="ioned_topol",input_structure="em",output_name="t_0",seed=seed,steps=2500000)

#analyse

AnalyseTrajectory.visualise(initial_sequence,charge)

#Cleanup code
##Bzips the tachyon input files
##Removes all trr files (if everything went well, we shouldn't need them...)

os.system("bzip2 *0.dat")
os.system("rm *0.dat")
os.system("rm *.trr")

