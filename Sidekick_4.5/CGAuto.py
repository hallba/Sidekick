#!/usr/bin/python

import sys,os,pdbio

from HAConf import configuration,programs,debug_status
from LipidBox import add_solvent_dppc
from MDPGenerator import replace_seed, replace_steps
import AutomatedPlot

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

os.system(gromacs+"grompp -f " + configuration['mdp_files'] + "cg-em_extended -c ioned.pdb -p ioned_topol -o em >& grompp_em.log")
os.system(gromacs+"mdrun -deffnm em")

if debug_status:
	replace_steps(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","runfile.mdp")
else:
	replace_seed(configuration['mdp_files'] + "cg-mdrun-50ns_sts.mdp","seeded.mdp",seed)
	replace_steps("seeded.mdp","runfile.mdp",10000000)

os.system(gromacs+"grompp -f runfile -c em -p ioned_topol -o t_0 >& grompp_md.log")
os.system(gromacs+"mdrun -deffnm t_0")

#analyse

length = len(initial_sequence)

if length %2 == 0:
	middle_residues = str(length/2 - 1) + "-" + str(length/2 +2)
else:
	middle_residues = str(length/2 - 1) + "-" + str(length/2 +1)

if charge != 0:
	make_ndx_command = 'echo "a b*\nr 1-7\n r ' + str(length-6) + '-' + str(length) +'\n r ' + middle_residues + '\n 16 & 17\n 16 & 18\n 16 & 19\n q\n" | '  + gromacs + 'make_ndx -f t_0 -o system.ndx'
	g_bundle_command = "echo '20\n21\n22\n'| " + gromacs + "g_bundle -f t_0 -s t_0 -na 1 -z -ok -n system.ndx "
else:
	make_ndx_command = 'echo "a b*\nr 1-7\n r ' + str(length-6) + '-' + str(length) +'\n r ' + middle_residues + '\n 15 & 16\n 15 & 17\n 15 & 18\n q\n" | '  + gromacs + 'make_ndx -f t_0 -o system.ndx'
        g_bundle_command = "echo '19\n20\n21\n'| " + gromacs + "g_bundle -f t_0 -s t_0 -na 1 -z -ok -n system.ndx "

os.system(make_ndx_command)
#print make_ndx_command
os.system(g_bundle_command)

AutomatedPlot.tilt_data_plot("bun_tilt.xvg","tilt_V_time.png",initial_sequence)
AutomatedPlot.hist_tilt_data_plot("bun_tilt.xvg","hist_tilt.png",initial_sequence)


AutomatedPlot.tilt_data_plot("bun_kink.xvg","kink_V_time.png",initial_sequence)
AutomatedPlot.hist_tilt_data_plot("bun_kink.xvg","hist_kink.png",initial_sequence)

#Nasty hack due to tempremental problem for cvs trjconv
os.system("touch electroneg.dat; touch elements.dat")
trjconv_xy_command = "echo '1\n1\n0\n' | " + gromacs + "trjconv_xy -f t_0.xtc -s em.gro -center -o xy_fit.xtc -fit rotxy+transxy"
os.system(trjconv_xy_command)

vmd_command = vmd + " em.gro xy_fit.xtc t_0.gro -dispdev text -eofexit <" + configuration['vmd_state_files'] + "automated_analysis.tcl > vmd.log"
os.system(vmd_command)

AutomatedPlot.hist_com_difference_plot("bilayer_position.dat","bilayer_position.png",initial_sequence)

AutomatedPlot.bilayer_deformation_refined("up.dat","down.dat","both-bully.png",graph_title=initial_sequence)

AutomatedPlot.rotation_plot("rot.dat","rotation.png",initial_sequence)

#Cleanup code
##Bzips the tachyon input files
##Removes all trr files (if everything went well, we shouldn't need them...)

os.system("bzip2 *0.dat")
os.system("rm *.trr")
