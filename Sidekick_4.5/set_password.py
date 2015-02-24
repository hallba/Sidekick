#!/usr/bin/python

from __future__ import with_statement
import sys, base64, pickle

password = base64.b64encode(sys.argv[1])
with open("pw.pickle","w") as oup:
	pickle.dump(password,oup)
	
