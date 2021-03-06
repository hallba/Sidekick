#/usr/bin/python

from os import system,remove
from HAConf import configuration

class add_solvent_dppc:
	def __init__(self,xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file,seed=None,position=0,apolar=True):	
		self.lipids = configuration['DPlipids']
		self.apolar = apolar
		self.fill_box(xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file,seed,position)

	def add_water(self,lipids,vdw_water,output_pdb_file,apolar=True):
		gromacs = configuration['gromacs_location']
		if apolar:
			final_genbox = gromacs + "genbox -cp prot+lipid.pdb -cs " + lipids + "water-minimized.pdb -vdwd " + str(vdw_water) + " -o " + output_pdb_file
		else:
			final_genbox = gromacs + "genbox -cp prot+lipid.pdb -cs " + lipids + "polarize-water.gro -vdwd " + str(0.16) + " -o " + output_pdb_file
		system(final_genbox)
	def fill_box(self,xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file,seed=None,position=0):
		vdw_lipid = 0.25
		vdw_water = 0.25
	
		gromacs = configuration['gromacs_location']
		lipids = self.lipids
		#system('rm " + lipids + "temp.pdb")
		#editconf_command = gromacs+"editconf -f "+ pdb_file_name + " -o temp.pdb -box " + str(xy) + " " + str(xy) + " " + str(shortz) + " -c"
		editconf_command = gromacs+"editconf -f "+ pdb_file_name + " -o temp_ini.pdb -box " + str(xy) + " " + str(xy) + " " + str(shortz) + " -c"
		system(editconf_command)
		editconf_command = gromacs+"editconf -f temp_ini.pdb -o temp.pdb -translate 0 0 " + str(position)
		system(editconf_command)
		
		
		for i in range(lipid_multiples):
			for j in range(64):
				if seed:
					genbox_command = gromacs + "genbox -cp temp.pdb -nmol 1 -nice 0 -seed " + str(seed) + " -vdwd " + str(vdw_lipid) + " -try 10000 -ci " + str(lipids) + "mol_" + str(j) + ".pdb -o temp.pdb >& genbox.log"
				else:
					genbox_command = gromacs + "genbox -cp temp.pdb -nmol 1 -nice 0 -vdwd " + str(vdw_lipid) + " -try 10000 -ci " + str(lipids) + "mol_" + str(j) + ".pdb -o temp.pdb >& genbox.log"
				#print genbox_command
				system(genbox_command)
			#Remove any unnecessary backups
			system("rm -r \#*")
	
		final_editconf = gromacs + "editconf -f temp.pdb -o prot+lipid.pdb -box " + str(xy) + " " + str(xy) + " " + str(z) + " -c"
		system(final_editconf)
		#final_genbox = gromacs + "genbox -cp prot+lipid.pdb -cs " + lipids + "water-minimized.pdb -vdwd " + str(vdw_water) + " -o " + output_pdb_file
		#system(final_genbox)
		self.add_water(lipids,vdw_water,output_pdb_file,self.apolar)
		system("rm -r \#*")
	
class add_solvent_dlpc(add_solvent_dppc):
	def __init__(self,xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file,seed=None,position=0,apolar=True):	
		self.lipids = configuration['DLlipids']
		self.apolar = apolar
		self.fill_box(xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file,seed,position)

class add_solvent_dopc(add_solvent_dppc):
	def __init__(self,xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file,seed=None,position=0,apolar=True):	
		self.lipids = configuration['DOlipids']
		self.apolar = apolar
		self.fill_box(xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file,seed,position)

class add_solvent_popc(add_solvent_dppc):
	def __init__(self,xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file,seed=None,position=0,apolar=True):	
		self.lipids = configuration['POlipids']
		self.apolar = apolar
		self.fill_box(xy,z,shortz,lipid_multiples,pdb_file_name,output_pdb_file,seed,position)

def write_dppg_itp():
	dppg = open("dppg.itp","w")
	dppg.write(""";;;;;; DIPALMITOYL PG  
;
; in general models PGs with saturated tail lengths C15-18

[moleculetype]
; molname 	nrexcl
  DPPG 		1

[atoms]
; id 	type 	resnr 	residu 	atom 	cgnr 	charge
  1 	Nda	1 	DPPG 	GLH 	1 	0.0 
  2 	Qa 	1 	DPPG 	PO4 	2 	-1.0 
  3 	Na 	1 	DPPG 	GL1 	3 	0 
  4 	Na 	1 	DPPG 	GL2 	4 	0 
  5 	C1 	1 	DPPG 	C1A 	5 	0 
  6 	C1 	1 	DPPG 	C2A 	6 	0 
  7 	C1 	1 	DPPG 	C3A 	7 	0 
  8 	C1 	1 	DPPG 	C4A 	8 	0 
  9 	C1 	1 	DPPG 	C1B 	9 	0 
  10 	C1 	1 	DPPG 	C2B 	10 	0 
  11 	C1 	1 	DPPG 	C3B 	11 	0 
  12 	C1 	1 	DPPG 	C4B 	12 	0 

[bonds]
; i j 	funct 	length 	force.c.
  1 2 	1 	0.470 	1250
  2 3 	1 	0.470 	1250
  3 4 	1 	0.370 	1250
  3 5 	1 	0.470 	1250
  5 6 	1 	0.470 	1250
  6 7 	1 	0.470 	1250
  7 8 	1 	0.470 	1250
  4 9 	1 	0.470 	1250
  9 10 	1 	0.470 	1250
  10 11 1 	0.470 	1250
  11 12 1 	0.470 	1250

[angles]
; i j k 	funct 	angle 	force.c.
  2 3 4 	2 	120.0 	25.0 
  2 3 5 	2 	180.0 	25.0 
  3 5 6 	2 	180.0 	25.0 
  5 6 7 	2 	180.0 	25.0 
  6 7 8 	2 	180.0 	25.0 
  4 9 10 	2 	180.0 	25.0 
  9 10 11 	2 	180.0 	25.0 
 10 11 12 	2 	180.0 	25.0 
""")
	dppg.close()