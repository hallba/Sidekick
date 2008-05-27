#!/usr/bin/python
#
#Sticking with nm for now. Remember the conversion

#For reference, Z is a special extra residue for acetyl groups which you might want to include in the analyses
##X isn't recognised

amino_converter = {	'ALA': 'A',			'ARG': 'R',			'ASN': 'N',
			'ASP': 'D',			'CYS': 'C',			'GLU': 'E',
			'GLN': 'Q',			'GLY': 'G',			'HIS': 'H',
			'ILE': 'I',			'LEU': 'L',			'LYS': 'K',
			'MET': 'M',			'PHE': 'F',			'PRO': 'P',
			'SER': 'S',			'THR': 'T',			'TRP': 'W',
			'TYR': 'Y',			'VAL': 'V',			'ACE': 'Z'
			}

class read_pdb:
	def __init__(self,filename):
		self.source = open(filename, "rb")
		self.fasta = ''
		self.natoms = 0
		self.box = 0
		
		self.master_list = []
		self.residue_list = []
		self.chain_list = []
		self.ca_list = []
		self.protein_list = []
		self.frame = []
				
		self.chain = {}
		self.residue = {}
		self.conformation = {}
		
		self.read_structure()
	def read_structure(self):
		residue = []
		#These are explicit counter to follow the indices
		count = 0
		residue_count = 0
		for line in self.source:
			if line[:6] == "ATOM  " or line[:6] == "HETATM":
				#print line
				[count,residue_count] = self.add_atom(line,count,residue_count)

		#Now, update the dictionaries to eg give a mapping of chain(A) to the appropriate list
		for chain_id in self.chain_list:
			self.chain[chain_id[0]["Chain_ID"]]= chain_id
				
	def add_atom(self,line,count,residue_count):
				try:
					residue_name = amino_converter[line[17:20]]
				except:
					residue_name = 'X'
				
				current_atom = {	'Number':	int(line[6:11]),
							'Name':		line[11:16],
							'Conformation':	line[17],
							'Residue_Name':	residue_name,
							'Residue_Number':int(line[22:26]),
							'Chain_ID':	line[21],
							'Index':	count,
							'Residue_Index':residue_count,
							'Beta' :	float(line[60:66])
						}
				self.master_list.append(current_atom)
				#Is this a new residue? If so, add the residue name to the fasta sequence
				if len(self.residue_list)==0:
					self.update_residue(current_atom)
					self.residue_list.append([])
					self.chain_list.append([])
					
				elif self.residue_list[-1][-1]["Residue_Number"] != current_atom["Residue_Number"]:
					self.update_residue(current_atom)
					self.residue_list.append([])
					residue_count += 1
					if self.chain_list[-1][-1]["Chain_ID"] != current_atom["Chain_ID"]:
						self.chain_list.append([])
				
				#Now add everything to the appropriate lists
				self.residue_list[-1].append(current_atom)
				self.chain_list[-1].append(current_atom)
				if current_atom["Name"] == "  CA " and residue_name != 'X': self.ca_list.append(current_atom)
				if residue_name != 'X' : self.protein_list.append(current_atom)
				self.frame.append(float(line[30:38])/10.0)
				self.frame.append(float(line[38:46])/10.0)
				self.frame.append(float(line[46:54])/10.0)
				count += 1				
				return count,residue_count
	
	def update_residue(self,current_atom):
		 self.fasta = self.fasta + current_atom['Residue_Name']
	
	#def get_indices(self,):
	#	[i["Index"] for i in pdb.chain["X"]]
