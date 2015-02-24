#!/usr/bin/python

import sys,os,copy
#sys.path.append('/Users/Shared/Sidekick')
import CartesianToolkit, HAConf
import GromacsInterface
'''
transform...
{-0.157036110759 -0.987567722797 -0.00704677961767 107.50088501}
{0.98730045557 -0.157159239054 0.0232125539333 7.39267349243}
{-0.0240314360708 -0.00331207970157 0.999705731869 2.89353561401}
{0.0 0.0 0.0 1.0}

transform = [	[ -0.157,	0.987,	-0.024,	0.000],
				[ -0.988,	-0.157,	-0.003,	0.000],
				[ -0.007,	0.023, 1.000,	0.000],
				[ 107.501,	7.393,	2.894,	1.000]	]

    3ALA      N   15   5.058   4.831   2.222
    3ALA      H   16   5.078   4.811   2.126
    3ALA     CA   17   5.174   4.815   2.313
    3ALA     CB   18   5.302   4.809   2.228
    3ALA      C   19   5.188   4.924   2.421
    3ALA      O   20   5.190   4.893   2.540
   translation = [-0.370002746582 3.061668396 1.54166793823]
   
ATOM    117  N   ALA    20      -1.510   0.238  -1.553  1.00  0.00
ATOM    118  H   ALA    20      -1.427  -0.026  -2.514  1.00  0.00
ATOM    119  CA  ALA    20      -2.249  -0.602  -0.635  1.00  0.00
ATOM    120  CB  ALA    20      -2.900  -1.809  -1.331  1.00  0.00
ATOM    121  C   ALA    20      -1.320  -1.141   0.409  1.00  0.00
ATOM    122  O   ALA    20      -1.684  -1.233   1.585  1.00  0.00

rotation
{-0.153550609946 -0.988060593605 -0.0125884236768 110.347702026}
{0.988135278225 -0.15357992053 0.00139002955984 8.75713348389}
{-0.00330676254816 -0.0122256251052 0.999919772148 0.801406860352}
{0.0 0.0 0.0 1.0}

100 deg
{-0.173648177663 -0.984807753013 0.0 0.0}
{0.984807753013 -0.173648177663 0.0 0.0}
{0.0 0.0 1.0 0.0}
{0.0 0.0 0.0 1.0}

transform (for a centered peptide)
{-0.153930529952 -0.988046526909 -0.00833086576313 -0.111345767975}
{0.987738966942 -0.154093712568 0.0250375568867 0.0237468481064}
{-0.0260220058262 -0.00437467638403 0.999651789665 1.50333404541}
{0.0 0.0 0.0 1.0}

reverse transform (for a centered peptide)
{-0.153928115964 0.987739622593 -0.026011519134 -0.00148594379425}
{-0.988046824932 -0.154091536999 -0.00438735820353 -0.0997638106346}
{-0.00834172219038 0.025025261566 0.999652028084 -1.50434374809}
{0.0 0.0 0.0 1.0}
atoms	 = { "N":	[50.58,   48.31,   22.22],
			 "H":	[50.78,   48.11,   21.26],
			 "CA":	[51.74,   48.15,   23.13],
			 "CB":  [53.02,   48.09,   22.28],
			 "C":	[51.88,   49.24,   24.21],
			 "O":	[51.90,   48.93,   25.40]
			}

'''	
transform = [	[ -0.157,	-0.987,	-0.07,	0.000],
				[ 0.988,	-0.157,	0.023,	0.000],
				[ -0.024,	-0.003, 1.000,	0.000],
				[ 0.000,	0.000,	0.000,	1.000]	]

#translation = [-0.370, 3.062, 1.542]
translation =  [-0.111,	0.024,	1.503]

atoms	 = { "N":	[-1.510,   0.238,  -1.553,"N"],
			 "H":	[-1.427,  -0.026,  -2.514,"H"],
			 "CA":	[-2.249,  -0.602,  -0.635,"CA"],
			 "CB":  [-2.900,  -1.809,  -1.331,"CB"],
			 "C":	[-1.320,  -1.141,   0.409,"C"],
			 "O":	[-1.684,  -1.233,   1.585,"O"],
			 "CG":	[-4.300,  -1.813,  -1.511,"CG"],
			 "CD2": [-5.098,  -1.119,  -2.413,"CD2"],
			 "ND1":	[-5.069,  -2.545,  -0.611,"ND1"],
			 "CE1": [-6.301,  -2.227,  -1.045,"CE1"],
			 "NE2": [-6.408,  -1.388,  -2.111,"NE2"],
			 "CG1":	[-3.879,  -1.256,  -2.643,"CG1"],
			 "CG2":	[-3.851,  -2.687,  -0.559,"CG2"],
			 "CD1":	[-4.350,  -2.332,  -3.644,"CD1"],
			 "LCG": [-3.918,  -1.368,  -2.596,"CG"],
			 "LCD1":[-4.560,  -2.626,  -3.198,"CD1"],
			 "LCD2":[-5.015,  -0.407,  -2.109,"CD2"],
			 "KCG": [-3.765,  -2.723,  -0.646,"CG"],
			 "KCD": [-4.370,  -3.896,  -1.423,"CD"],
			 "KCE":	[-5.140,  -4.802,  -0.454,"CE"],
			 "KNZ":	[-5.865,  -5.832,  -1.216,"NZ"],
			 "MCG":	[-3.962,  -1.317,  -2.519,"CG"],
			 "MSD": [-4.826,  -2.758,  -3.165,"SD"],
			 "MCE": [-5.646,  -1.974,  -4.559,"CE"],
			 "FCG":	[-3.818,  -1.415,  -2.617,"CG"],
			 "FCD1":[-3.306,  -1.165,  -3.895,"CD1"],
			 "FCD2":[-5.185,  -1.253,  -2.379,"CD2"],
			 "FCE1":[-4.148,  -0.740,  -4.918,"CE1"],
			 "FCE2":[-6.029,  -0.837,  -3.407,"CE2"],
			 "FCZ": [-5.509,  -0.578,  -4.674,"CZ"],
			 "PCB": [-2.781,  -1.674,  -1.634,"CB"], 
			 "PCG":	[-1.729,  -1.706,  -2.752,"CG"],
			 "PCD": [-1.322,  -0.237,  -2.897,"CD"],
			 "OG":	[-2.024,  -2.685,  -2.025,"OG"],
			 "TCG2":[-3.729,  -2.809,  -0.652,"CG2"],
			 "OG1": [-3.856,  -1.250,  -2.401,"OG1"],
			 "WCG": [-3.940,  -1.353,  -2.486,"CG"],
			 "WCD1":[-3.741,  -1.096,  -3.859,"CD1"],
			 "WCD2":[-5.263,  -1.052,  -2.233,"CD2"],
			 "WCE2":[-5.851,  -0.622,  -3.449,"CE2"],
			 "WCE3":[-6.023,  -1.086,  -1.036,"CE3"],
			 "WNE1":[-4.922,  -0.642,  -4.478,"NE1"],
			 "WCZ2":[-7.207,  -0.233,  -3.477,"CZ2"],
			 "WCZ3":[-7.363,  -0.702,  -1.092,"CZ3"],
			 "WCH2":[-7.949,  -0.285,  -2.294,"CH2"],
			 "YCG":	[-3.876,  -1.471,  -2.539,"CG"],
			 "YCD1":[-5.195,  -1.150,  -2.202,"CD1"],
			 "YCD2":[-3.487,  -1.467,  -3.881,"CD2"],
			 "YCE1":[-6.110,  -0.823,  -3.200,"CE1"],
			 "YCE2":[-4.402,  -1.141,  -4.876,"CE2"],
			 "YCZ":	[-5.715,  -0.819,  -4.534,"CZ"],
			 "YOH":	[-6.615,  -0.496,  -5.509,"OH"],
			 "RCG":	[-3.936,  -2.634,  -0.762,"CG"],
			 "RCD":	[-4.497,  -3.749,  -1.662,"CD"],
			 "RNE":	[-5.506,  -4.546,  -0.903,"NE"],
			 "RCZ":	[-6.179,  -5.588,  -1.381,"CZ"],
			 "RNH1":[-6.030,  -6.050,  -2.591,"NH1"],
			 "RNH2":[-7.031,  -6.179,  -0.604,"NH2"],
			 "NCG":	[-3.834,  -1.353,  -2.625,"CG"],
			 "NND2":[-5.053,  -0.986,  -2.330,"ND2"],
			 "NOD1":[-3.478,  -1.334,  -3.795,"OD1"],
			 "DCG":	[-3.789,  -2.730,  -0.649,"CG"],
			 "DOD1":[-4.467,  -2.285,   0.301,"OD1"],
			 "DOD2":[-3.754,  -3.945,  -0.929,"OD2"],
			 "SG":	[-3.985,  -2.757,  -0.489,"SG"],
			 "QCG":	[-3.749,  -2.790,  -0.720,"CG"],
			 "QCD":	[-4.359,  -3.975,  -1.481,"CD"],
			 "QNE2":[-4.159,  -4.098,  -2.768,"NE2"],
			 "QOE1":[-5.012,  -4.827,  -0.900,"OE1"],
			 "ECG":	[-3.759,  -2.785,  -0.671,"CG"],
			 "ECD":	[-4.367,  -3.965,  -1.407,"CD"],
			 "EOE1":[-4.201,  -4.111,  -2.633,"OE1"],
			 "EOE2":[-5.021,  -4.774,  -0.717,"OE2"],
			 
			}
			
amino_converter = {	'ALA': 'A',			'ARG': 'R',			'ASN': 'N',
					'ASP': 'D',			'CYS': 'C',			'GLU': 'E',
					'GLN': 'Q',			'GLY': 'G',			'HIS': 'H',
					'ILE': 'I',			'LEU': 'L',			'LYS': 'K',
					'MET': 'M',			'PHE': 'F',			'PRO': 'P',
					'SER': 'S',			'THR': 'T',			'TRP': 'W',
					'TYR': 'Y',			'VAL': 'V',			'ACE': 'Z'
					}

fasta_converter = {}

for three_letter_code in amino_converter.keys():
	fasta_converter[amino_converter[three_letter_code]] = three_letter_code

residue_atoms = 	{	"A":	["N","CA","C","O","CB","H"],
						"G":	["N","CA","C","O","H"],
						"H":	["N","CA","C","O","CB","CG","CD2","ND1","CE1","NE2","H"],
						"I":	["N","CA","C","O","CB","CG1","CG2","CD1","H"],
						"L":	["N","CA","C","O","CB","LCG","LCD1","LCD2","H"],
						"K":	["N","CA","C","O","CB","KCG","KCD","KCE","KNZ","H"],
						"M":	["N","CA","C","O","CB","MCG","MSD","MCE","H"],
						"F":	["N","CA","C","O","CB","FCG","FCD1","FCD2","FCE1","FCE2","FCZ","H"],
						"P":	["N","CA","C","O","PCB","PCG","PCD"],
						"S":	["N","CA","C","O","CB","OG","H"],
						"T":	["N","CA","C","O","CB","CG2","OG1","H"],
						"W":	["N","CA","C","O","CB","WCG","WCD1","WCD2","WCE2","WCE3","WNE1","WCZ2","WCZ3","WCH2","H"],
						"Y":	["N","CA","C","O","CB","YCG","YCD1","YCD2","YCE1","YCE2","YCZ","YOH","H"],
						"V":	["N","CA","C","O","CB","CG1","CG2","H"],
						"R":	["N","CA","C","O","CB","RCG","RCD","RNE","RCZ","RNH1","RNH2","H"],
						"N":	["N","CA","C","O","CB","NCG","NND2","NOD1","H"],
						"D":	["N","CA","C","O","CB","DCG","DOD2","DOD1","H"],
						"C":	["N","CA","C","O","CB","SG","H"],
						"Q":	["N","CA","C","O","CB","QCG","QCD","QNE2","QOE1","H"],
						"E":	["N","CA","C","O","CB","ECG","ECD","EOE2","EOE1","H"],
						}



def BuildHelix(sequence,filename="tm.pdb"):
	#ln the opls ff here
	os.system("ln -s " + HAConf.sidekick_location + "/gromacs/share/gromacs/top/oplsaa.ff `pwd`" )
	os.system("ln -s " + HAConf.sidekick_location + "/gromacs/share/gromacs/top/residuetypes.dat")
	used_atoms = copy.deepcopy(atoms)
	atom_index = 0
	residue_index = 0
	#First I build a raw, suboptimal helix and write that out
	initial_structure = open("initial_helix.pdb","w")
	print >> initial_structure, "TITLE     GRowing Old MAkes el Chrono Sweat\nREMARK    THIS IS A SIMULATION BOX\nCRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1           1\nMODEL        1\n"

	for residue in sequence:
		residue_index += 1
		for atom in residue_atoms[residue]:
			atom_index += 1
			print >> initial_structure, "ATOM  %5d  %-4s%3s%6d    %8.3f%8.3f%8.3f  1.00  0.00" % (atom_index, used_atoms[atom][3], fasta_converter[residue], residue_index, used_atoms[atom][0], used_atoms[atom][1], used_atoms[atom][2] )
		#update atom positions by multiplication with rotation matrix
		atom_names = used_atoms.keys()
		for name in atom_names:
			used_atoms[name][0:3] = CartesianToolkit.vecadd(used_atoms[name][0:3],translation)
			#used_atoms[name] = CartesianToolkit.rotate_by_matrix(used_atoms[name],rotation)
			used_atoms[name][0:3] = CartesianToolkit.transform_by_rot_trans_matrix(used_atoms[name][0:3],transform)
	initial_structure.close()
	#Now clear up any clashes by 100 steps of steepest descent in gromacs
	mdp_file = open("quick_sd.mdp","w")
	print >> mdp_file, "integrator = steep\nnsteps = 100"
	mdp_file.close()
	#[pdb2gmx_stdin, pdb2gmx_stout_sterr] = os.popen4(HAConf.configuration['gromacs_location'] + "pdb2gmx -f initial_helix.pdb -p initial_helix -o initial_helix")
	#Select gmx53a6
	#print >> pdb2gmx_stdin, "9"
	#print >> pdb2gmx_stdin, "3"
	#pdb2gmx_stdin.flush()
	#for line in pdb2gmx_stout_sterr: pass
	GromacsInterface.pdb2gmx(logfile="pdb2gmx.log")
	GromacsInterface.grompp(logfile="grompp.log")
	#[grompp_stdin, grompp_stout_sterr] = os.popen4(HAConf.configuration['gromacs_location'] + "grompp -f quick_sd.mdp -c initial_helix -p initial_helix -o helix_em")
	#for line in grompp_stout_sterr: pass
	#os.system(HAConf.configuration['gromacs_location'] + "grompp -f quick_sd.mdp -c initial_helix -p initial_helix -o helix_em")
	#for line in grompp_stout_sterr: print
	GromacsInterface.mdrun()
	#[mdrun_stdin, mdrun_stout_sterr] = os.popen4(HAConf.configuration['gromacs_location'] + "mdrun -deffnm helix_em")
	#for line in mdrun_stout_sterr: pass
	if os.path.exists("helix_em.gro"):
		[editconf_stdin, editconf_stout_sterr] = os.popen4(HAConf.configuration['gromacs_location'] + "editconf -f helix_em.gro -o " + filename)
	else:
		print "Problem in atomistic helix generation. Continuing nevertheless."
		[editconf_stdin, editconf_stout_sterr] = os.popen4(HAConf.configuration['gromacs_location'] + "editconf -f helix_em.tpr -o " + filename)
	for line in editconf_stout_sterr: pass
	print "Generated alpha helix with", sequence
	
	
	
if __name__ == "__main__":
	import sys
	try:
		sequence = sys.argv[1]
	except:
		print "Requires a sequence to work on"
		exit()
	BuildHelix(sequence)
