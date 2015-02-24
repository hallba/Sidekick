class universal_get_options:
	def __init__(self,required_options):
		#global options
		#Sort out the options
		from optparse import OptionParser
	
		class careful_OptionParser (OptionParser):
		   
		  def check_required (self, opt):
		    option = self.get_option(opt)
		   
		    # Assumes the option's 'default' is set to None!
		    if (getattr(self.values, option.dest) == None):
		      self.error("%s option not supplied" % option)
	
		self.parser = careful_OptionParser()
		self.get_generics()
		self.get_specifics()
		(self.options, self.args) = self.parser.parse_args()
		for req_option in required_options:
			self.parser.check_required(req_option)
	def get_generics(self):
		'''All jobs have these options'''
		self.parser.add_option("-t", "--type",
			dest="model_type", default="MARTINI",
			help="Type of CG model to use, or 'search' for all models in a multi-batch job. Options are:\n\tMARTINI = MARTINI 2.1 (default)\n\tMARTINI2.1.1 = MARTINI 2.1.1\n\tBond\n\tBond0.9.5")
		self.parser.add_option("-j","--job_name",
			dest="destination", default=None,
			help="Location of the results folder on the NFS share (defined by HAConf)")
		
	def get_speficics(self):
		'''The universal get_options lacks these'''
		pass

class Analysis_Options(universal_get_options):
	def get_specifics(self):
		self.parser.add_option("-f", "--file", dest="activity_file", default=None,
	         help="File with sequences and activities on seperate lines", metavar="FILE")
		self.parser.add_option("-o", "--output", dest="summary_file", default=None,
	         help="CSV File with sequences and metrics on seperate lines", metavar="FILE")
		
class Helix_Options(universal_get_options):
	def get_specifics(self):
		self.parser.add_option("-s", "--sequence", dest="sequence", default=None,
		     help="FASTA Sequence")
		self.parser.add_option("-b", "--batch-mode",
			 action="store_true", dest="batch", default=False,
			 help="Batch mode creates results in an appropriate directory structure")
		self.parser.add_option("-l", "--length",
			 dest="simlength", default=100,type="int",
			 help="Generate a given number of ns (default = 100ns)")
		self.parser.add_option("-u", "--unbias",
		     dest="bias", action="store_false", default=True,
			 help="Unbias the membrane formation process by spreading lipids through the length of the box")
		self.parser.add_option("-r", "--random-seed",
		     dest="seed", default=5,type="int",
		     help="Random seed (default = 5)")
		self.parser.add_option("--randomize_lipids",
			 dest="random_lipids", default=False, action="store_true",
			 help="Generate different lipid positions for each different seed (by default, only the mdrun seeds are different)")
		self.parser.add_option("-a", "--angle",
		     dest="angle", default=0,type="int",
		     help="Initial angle for helix (default = 0)")
		self.parser.add_option("--change_temperature",
		     dest="temperature", default=None,type="int",
		     help="Temperature in K for simulation to run at (default 323k)")
		self.parser.add_option("-p", "--position",
		     dest="position", default=0,type="float",
		     help="Initial position for helix relative to box center in nm (default = 0)")
		self.parser.add_option("--system_size",
			 dest="system_size", default="XS",
			 help="System size: XS=128 lipids, S=256 lipids")
		self.parser.add_option("--preformed_bilayer",
			 dest="preformed_insertion", default=False, action="store_true",
			 help="Uses a preformed DPPC bilayer and places the helix at 90' to it, 35A away from the center of the bilayer. This option overrides -p, -u and -a")
		self.parser.add_option("--special",dest='special',default=None,help="Special options to pass to the script. Failure to ensure that the scripts accept the options can cause the jobs not to run")
		self.parser.add_option("--lipid_headgroup_mix",
				dest='headgroup_mix',default="DPPC",
				help="Mix of DP lipids in the system. Format 'DPPG:DPPC/10:90'. Defaults to DPPC")
		self.parser.add_option("--wallclock",
			 dest="wallclock", default=48,type="int",
			 help="Maximum length simulation should run for in hours (default=48)")
		     
class Multi_Helix_Options(universal_get_options):
	def get_specifics(self):
		self.parser.add_option("-s", "--sequence", dest="sequence", default=None,
		          help="FASTA Sequence")
		self.parser.add_option("-b", "--batch-mode",
			 action="store_true", dest="batch", default=False,
			 help="Batch mode creates results in an appropriate directory structure")
		self.parser.add_option("-l", "--length",
			 dest="simlength", default=100,type="int",
			 help="Generate a given number of ns (default = 100ns)")
		self.parser.add_option("-u", "--unbias",
		     dest="bias", action="store_false", default=True,
			 help="Unbias the membrane formation process by spreading lipids through the length of the box")
		self.parser.add_option("-m", "--mutate",
			 dest="mutant",default=None,
			 help = "Residue to scan through the given sequence")
		self.parser.add_option("-f", "--file", dest="sequence_file", default=None,
	         help="File with sequences to be simulated on seperate lines", metavar="FILE")
		self.parser.add_option("-d", "--duplicates",
	         dest="duplicates", default=10,type="int",
	         help="Number of duplicates (default = 5)")
		self.parser.add_option("--stripe",
			 action="store_true", dest="stripe", default=False,
			 help="Stripe mode performs all sequences a seed at a time, instead of each sequence in order")
		self.parser.add_option("--randomize_lipids",
			 dest="random_lipids", default=False, action="store_true",
			 help="Generate different lipid positions for each different seed (by default, only the mdrun seeds are different)")
		self.parser.add_option("--vary_angle",
			 dest="vary_angle", default=False, action="store_true",
			 help="Generate variable initial helix angles (in 15' increments)")
		self.parser.add_option("--vary_position",
			 dest="vary_position", default=False, action="store_true",
			 help="Generate variable initial helix positions (in 1.5nm increments)")
		self.parser.add_option("--system_size",
			 dest="system_size", default="XS",
			 help="System size: XS=128 lipids, S=256 lipids")
		self.parser.add_option("--art",
			 dest="art", default=False, action="store_true",
			 help="Use ART tools for grid computations")
		self.parser.add_option("--special",dest='special',default=None,help="Special options to pass to the script")
		self.parser.add_option("--wallclock",
			 dest="wallclock", default=48,type="int",
			 help="Maximum length simulation should run for in hours (default=48)")
		self.parser.add_option("--binary_lipid_scan",
			 dest="binary_lipid_scan", default=None,
			 help="Binary lipid scan. Format= DPPC:DPPG/min_%DPPC:step:max_%DPPC. If set, the lipid headgroup mix is ignored")
		self.parser.add_option("--lipid_headgroup_mix",
				dest='headgroup_mix',default="DPPC",
				help="Mix of DP lipids in the system. Format 'DPPG:DPPC/10:90'. Defaults to DPPC")
		self.parser.add_option("--change_temperature",
		     dest="temperature", default=None,type="int",
		     help="Temperature in K for simulation to run at (default 323k)")

class Multi_UmbrellaSA_Helix_Options(universal_get_options):
	def get_specifics(self):
		self.parser.add_option("-s", "--sequence", dest="sequence", default=None,
		          help="FASTA Sequence")
		self.parser.add_option("-b", "--batch-mode",
			 action="store_true", dest="batch", default=False,
			 help="Batch mode creates results in an appropriate directory structure")
		self.parser.add_option("-l", "--length",
			 dest="simlength", default=100,type="int",
			 help="Generate a given number of ns (default = 100ns)")
		self.parser.add_option("-u", "--unbias",
		     dest="bias", action="store_false", default=True,
			 help="Unbias the membrane formation process by spreading lipids through the length of the box")
		self.parser.add_option("-m", "--mutate",
			 dest="mutant",default=None,
			 help = "Residue to scan through the given sequence")
		self.parser.add_option("-f", "--file", dest="sequence_file", default=None,
	         help="File with sequences to be simulated on seperate lines", metavar="FILE")
		self.parser.add_option("-d", "--duplicates",
	         dest="duplicates", default=100,type="int",
	         help="Number of duplicates per window (default = 100)")
		self.parser.add_option("--stripe",
			 action="store_true", dest="stripe", default=False,
			 help="Stripe mode performs all sequences a seed at a time, instead of each sequence in order")
		self.parser.add_option("--randomize_lipids",
			 dest="random_lipids", default=False, action="store_true",
			 help="Generate different lipid positions for each different seed (by default, only the mdrun seeds are different)")
		self.parser.add_option("--dont_vary_angle",
			 dest="vary_angle", default=True, action="store_false",
			 help="Generate variable initial helix angles (in 15' increments)")
		self.parser.add_option("--window_size",
			 dest="window_size", default=0.5, type="float",
			 help="Size of umbrella windows (nm)")
		self.parser.add_option("--window_max",
			 dest="window_max", default=6.50, type="float",
			 help="Maximum distance from bilayer COM (nm)")
		self.parser.add_option("--system_size",
			 dest="system_size", default="XS",
			 help="System size: XS=128 lipids, S=256 lipids")
		self.parser.add_option("--art",
			 dest="art", default=False, action="store_true",
			 help="Use ART tools for grid computations")
		self.parser.add_option("--special",dest='special',default=None,help="Special options to pass to the script")
		self.parser.add_option("--wallclock",
			 dest="wallclock", default=48,type="int",
			 help="Maximum length simulation should run for in hours (default=48)")
		self.parser.add_option("--binary_lipid_scan",
			 dest="binary_lipid_scan", default=None,
			 help="Binary lipid scan. Format= DPPC:DPPG/min_%DPPC:step:max_%DPPC. If set, the lipid headgroup mix is ignored")
		self.parser.add_option("--lipid_headgroup_mix",
				dest='headgroup_mix',default="DPPC",
				help="Mix of DP lipids in the system. Format 'DPPG:DPPC/10:90'. Defaults to DPPC")
		self.parser.add_option("--change_temperature",
		     dest="temperature", default=None,type="int",
		     help="Temperature in K for simulation to run at (default 323k)")


class PMF_Helix_Options(universal_get_options):
	def get_specifics(self):
		self.parser.add_option("-s", "--sequence", dest="sequence", default=None,
		     help="FASTA Sequence")
		self.parser.add_option("-b", "--batch-mode",
			 action="store_true", dest="batch", default=False,
			 help="Batch mode creates results in an appropriate directory structure")
		self.parser.add_option("-l", "--length",
			 dest="simlength", default=100,type="int",
			 help="Generate a given number of ns (default = 100ns)")
		self.parser.add_option("-r", "--random-seed",
		     dest="seed", default=5,type="int",
		     help="Random seed (default = 5)")
		self.parser.add_option("-w", "--window",
			 dest="window", default=0., type="float",
			 help="Window position (Angstroms)")
		self.parser.add_option("--system_size",
			 dest="system_size", default="XS",
			 help="System size: XS=128 lipids, S=256 lipids")
class Multi_PMF_Helix_Options(universal_get_options):
	def get_specifics(self):
		self.parser.add_option("-s", "--sequence", dest="sequence", default=None,
		          help="FASTA Sequence")
		self.parser.add_option("-b", "--batch-mode",
			 action="store_true", dest="batch", default=False,
			 help="Batch mode creates results in an appropriate directory structure")
		self.parser.add_option("-l", "--length",
			 dest="simlength", default=100,type="int",
			 help="Generate a given number of ns (default = 100ns)")
		self.parser.add_option("-u", "--unbias",
		     dest="bias", action="store_false", default=True,
			 help="Unbias the membrane formation process by spreading lipids through the length of the box")
		self.parser.add_option("-m", "--mutate",
			 dest="mutant",default=None,
			 help = "Residue to scan through the given sequence")
		self.parser.add_option("-f", "--file", dest="sequence_file", default=None,
	         help="File with sequences to be simulated on seperate lines", metavar="FILE")
		self.parser.add_option("-w", "--window",
			 dest="window", default=0., type="float",
			 help="Window position (Angstroms)")

class Helix_Dimer_Options(universal_get_options):
	def get_specifics(self):
		self.parser.add_option("--s1", "--sequence_1", dest="sequence1", default=None,
		     help="FASTA Sequence for first helix (order is irrelevant)")
		self.parser.add_option("--s2", "--sequence_2", dest="sequence2", default=None,
		     help="FASTA Sequence for second helix (order is irrelevant)")
		self.parser.add_option("-b", "--batch-mode",
			 action="store_true", dest="batch", default=False,
			 help="Batch mode creates results in an appropriate directory structure")
		self.parser.add_option("-l", "--length",
			 dest="simlength", default=100,type="int",
			 help="Generate a given number of ns (default = 100ns)")
		self.parser.add_option("-r", "--random-seed",
		     dest="seed", default=5,type="int",
		     help="Random seed (default = 5)")
		self.parser.add_option("--change_temperature",
		     dest="temperature", default=None,type="int",
		     help="Temperature in K for simulation to run at (default 323k)")
		self.parser.add_option("-a","--antiparallel",
			 dest="parallel", default=True, action="store_false",
			 help="Switch helices to antiparallel orientation. Defaults to parallel")
		self.parser.add_option("--special",dest='special',default=None,help="Special options to pass to the script. Failure to ensure that the scripts accept the options can cause the jobs not to run")
		self.parser.add_option("--wallclock",
			 dest="wallclock", default=100,type="int",
			 help="Maximum length simulation should run for in hours (default=100)")
		self.parser.add_option("--rotatebyangle",
			 dest="rotatebyangle", default=0,type="int",
		     help="Additional rotation angle for helix about z (default = 0)")

class Helix_Oligomer_Options(universal_get_options):
	def get_specifics(self):
		self.parser.add_option("-s", "--sequences", dest="sequences", default=None,
		     help="Comma separated list of FASTA Sequences of helices ")
		self.parser.add_option("-b", "--batch-mode",
			 action="store_true", dest="batch", default=False,
			 help="Batch mode creates results in an appropriate directory structure")
		self.parser.add_option("-l", "--length",
			 dest="simlength", default=2500,type="int",
			 help="Generate a given number of ns (default = 2500ns)")
		self.parser.add_option("-r", "--random-seed",
		     dest="seed", default=5,type="int",
		     help="Random seed (default = 5)")
		self.parser.add_option("--change_temperature",
		     dest="temperature", default=None,type="int",
		     help="Temperature in K for simulation to run at (default 323k)")
		self.parser.add_option("-a","--alternating",
			 dest="alternating", default=False, action="store_true",
			 help="Alternate orientations of the peptides (ie antiparallel sequences in order, as opposed to all parallel")
		self.parser.add_option("--special",dest='special',default=None,help="Special options to pass to the script. Failure to ensure that the scripts accept the options can cause the jobs not to run")
		self.parser.add_option("--wallclock",
			 dest="wallclock", default=2500,type="int",
			 help="Maximum length simulation should run for in hours (default=2500)")

class Multi_Helix_Dimer_Options(universal_get_options):
	def get_specifics(self):
		self.parser.add_option("--s1", "--sequence1", dest="sequence1", default=None,
		          help="FASTA Sequence for one helix")
		self.parser.add_option("--s2", "--sequence2", dest="sequence2", default=None,
		          help="FASTA Sequence for second helix")
		self.parser.add_option("-b", "--batch-mode",
			 action="store_true", dest="batch", default=False,
			 help="Batch mode creates results in an appropriate directory structure")
		self.parser.add_option("-l", "--length",
			 dest="simlength", default=500,type="int",
			 help="Generate a given number of ns (default = 500ns)")
		self.parser.add_option("-f", "--file", dest="sequence_file", default=None,
	         help="File with sequences to be simulated on seperate lines. Both sequences must be on the same line, space seperated", metavar="FILE")
		self.parser.add_option("-d", "--duplicates",
	         dest="duplicates", default=100,type="int",
	         help="Number of duplicates (default = 100)")
		self.parser.add_option("--stripe",
			 action="store_true", dest="stripe", default=False,
			 help="Stripe mode performs all sequences a seed at a time, instead of each sequence in order")
		self.parser.add_option("--special",dest='special',default=None,help="Special options to pass to the CG script")
		self.parser.add_option("--wallclock",
			 dest="wallclock", default=320,type="int",
			 help="Maximum length simulation should run for in hours (default=320)")
		self.parser.add_option("--change_temperature",
		     dest="temperature", default=None,type="int",
		     help="Temperature in K for simulation to run at (default 323k)")
                self.parser.add_option("-a","--antiparallel",
                         dest="parallel", default=True, action="store_false",
                         help="Switch helices to antiparallel orientation. Defaults to parallel")
		self.parser.add_option("--vary_angle",
			 dest="vary_angle", default=False, action="store_true",
			 help="Generate variable initial helix rotation angles (in 15' increments)")
		
'''Old code- deprecated in favour of code which better differentiates between the running program'''

def get_options(required_options):
	#global options
	#Sort out the options
	from optparse import OptionParser

	class careful_OptionParser (OptionParser):
	   
	  def check_required (self, opt):
	    option = self.get_option(opt)
	   
	    # Assumes the option's 'default' is set to None!
	    if (getattr(self.values, option.dest) == None):
	      self.error("%s option not supplied" % option)

	parser = careful_OptionParser()


	#Input for CG_Helix.py
	parser.add_option("-s", "--sequence", dest="sequence", default=None,
		          help="FASTA Sequence")
	parser.add_option("-r", "--random-seed",
		          dest="seed", default=5,type="int",
		          help="Random seed (default = 5)")
	parser.add_option("-b", "--batch-mode",
			action="store_true", dest="batch", default=False,
			help="Batch mode creates results in an appropriate directory structure")
	parser.add_option("-l", "--length",
			dest="simlength", default=100,type="int",
			help="Generate a given number of ns (default = 100ns)")
	parser.add_option("-t", "--type",
			dest="model_type", default="MARTINI",
			help="Type of CG model to use, or 'search' for all models in a multi-batch job. Options are:\n\tMARTINI = MARTINI 2.1 (default)\n\tMARTINI2.1.1 = MARTINI 2.1.1\n\tBond\n\tBond0.9.5")
	parser.add_option("-j","--job_name",
			dest="destination", default=None,
			help="Location of the results folder on the NFS share (defined by HAConf)")
	parser.add_option("-u", "--unbias",
		    dest="bias", action="store_false", default=True,
			help="Unbias the membrane formation process by spreading lipids through the length of the box")
	
	#Input for xgrid submission programs
	parser.add_option("-f", "--file", dest="sequence_file", default=None,
		          help="File with sequences to be simulated on seperate lines", metavar="FILE")
	parser.add_option("-d", "--duplicates",
		          dest="duplicates", default=10,type="int",
		          help="Number of duplicates (default = 5)")
	parser.add_option("-m", "--mutate",
			dest="mutant",default=None,
			help = "Residue to scan through the given sequence")

	(options, args) = parser.parse_args()
	
	for req_option in required_options:
	   parser.check_required(req_option)
	return options

