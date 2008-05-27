#/usr/bin/python

from os import system
from HAConf import configuration

def add_solvent_dppc(xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file):
	vdw_lipid = 0.25
	vdw_water = 0.25

	gromacs = configuration['gromacs_location']
	lipids = configuration['lipids']
	#system('rm " + lipids + "temp.pdb")
	editconf_command = gromacs+"editconf -f "+ pdb_file_name + " -o temp.pdb -box " + str(xy) + " " + str(xy) + " " + str(shortz) + " -c"
	system(editconf_command)
	
	
	for i in range(lipid_multiples):
		for j in range(64):
			genbox_command = gromacs + "genbox -cp temp.pdb -nmol 1 -nice 0 -vdwd " + str(vdw_lipid) + " -try 1000 -ci " + str(lipids) + "mol_" + str(j) + ".pdb -o temp.pdb >& genbox.log"
			#print genbox_command
			system(genbox_command)

	final_editconf = gromacs + "editconf -f temp.pdb -o prot+lipid.pdb -box " + str(xy) + " " + str(xy) + " " + str(z) + " -c"
	system(final_editconf)
	final_genbox = gromacs + "genbox -cp prot+lipid.pdb -cs " + lipids + "water-minimized.pdb -vdwd " + str(vdw_water) + " -o " + output_pdb_file
	system(final_genbox)
	system("rm -r \#*")
	
def write_dppg_itp():
	dppg = open("dppg.itp","w")
	dppg.write(""";;;;;; DIPALMITOYL PG  
	dppg.close()