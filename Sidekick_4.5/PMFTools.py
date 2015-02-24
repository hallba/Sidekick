#!/usr/bin/PMFTools

from __future__ import with_statement
import HAConf, GenerateCGSystem
import os
gromacs = HAConf.configuration['gromacs_location']

def UmbrellaGenerate(name="umbrella",window=None):
	#Need to generate 2 files; an index file with a group representing lipids and backbone, and a ppa file
	lipid_mixture = ""
	for lipid in GenerateCGSystem.allowed_lipids:
		lipid_mixture += " %s" % lipid
	make_ndx_command = 'echo "keep 0\na b* CA\nr %s\nname 2 Lipid\nq\n" | %s/make_ndx -f em.gro -o %s.ndx' % (lipid_mixture,gromacs,name)
	os.system(make_ndx_command)
	with open(name+".ppa","w") as ppa:
		print >> ppa, "verbose = no"
		print >> ppa, "runtype = umbrella"
		print >> ppa, "group_1 = B*_CA"
		print >> ppa, "reference_group = Lipid"
		print >> ppa, "pulldim = N N Y"
		print >> ppa, "k1 = 1000"
		print >> ppa, "pos1 = 0.0 0.0 %f" % window
	