#!/usr/bin/python

from __future__ import with_statement
import os

os.rename("protein-cg.itp","original_protein-cg.itp")


with open("protein-cg.itp","w") as oup:
	with open("original_protein-cg.itp") as inp:
		for line in inp:
			if "LYS\tSHQd" in line:
			 	contents = line.split("\t")
				contents[6] = "0.000"
				new_line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t" % tuple(contents)
				print >> oup, new_line
			else:
			 	print >> oup, line,
