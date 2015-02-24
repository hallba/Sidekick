#!/usr/bin/python

#gnm_patch- adds extra bonds to an itp file to allow for "stiffening" or more appropriately restraining a system

from __future__ import with_statement

import os

def vecdist(cart_list_a,cart_list_b):
	relating_vector = vecsub(cart_list_a, cart_list_b)
	return veclength(relating_vector)

def vecsub(cart_list_a, cart_list_b):
	zipped_carts = zip(cart_list_a, cart_list_b)
	return [(cart[0] - cart[1]) for cart in zipped_carts]

def veclength(cart_list):
	squares = [cart ** 2 for cart in cart_list]
	return float(sum(squares))**0.5

def backup_file(filename):
	if os.path.exists(filename):
		target_name = "#" + filename
		failure = True
		if not os.path.exists(target_name):
			os.rename(filename,target_name)
			alt_target_name = target_name
			failure = False
		else:
			for i in range(20):
				alt_target_name = target_name + "." + str(i)
				if os.path.exists(alt_target_name):
					continue
				else:
					os.rename(filename,alt_target_name)
					failure = False
					break
		if failure:
			print "ITP did not exist *or* Too many backups. Clean up and try again"
			exit()
		return alt_target_name

def kirchoff_generation(positions,upper_cutoff,lower_cutoff):
   #contacts = 0
   natoms = len(positions)
   #print "natoms", natoms
   contacts = []
   #matrix = [[0.0]*natoms for i in range(natoms)]
   for i in range(natoms):
   	#i+2 to avoid neighbouring atom bonds
   	for j in range(i+2,natoms):
   		distance_angstroms = vecdist(positions[i],positions[j])
   		if distance_angstroms >= lower_cutoff and distance_angstroms <= upper_cutoff:
   			contacts.append([i,j,distance_angstroms])

   return contacts

def construct_gnm(contacts,ca_indices, forceconstant):
	itp_string = ""
	for bond in contacts:
		itp_string += "%-6d %-6d %-6d %-12.8f %-8d; GNM %6d->%6d\n" % (ca_indices[bond[0]], ca_indices[bond[1]], 1, bond[2]/10, forceconstant, bond[0]+1, bond[1]+1 )
	return itp_string

def get_positions(pdbfilename):
	ca_indices = []
	positions = []
	with open(pdbfilename) as pdbdata:
		index = 0
		for line in pdbdata:
			if line[:4] != "ATOM":
				continue
			index += 1
			atom_name = line[11:16].strip()
			if atom_name =="CA" or ((atom_name[:1] == "B" or atom_name[1:2] == "B") and atom_name[:1] != "C"): #Bondini, atomistic or MARTINI
				positions.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
				ca_indices.append(index)

	return ca_indices, positions
	
def itp_append(olditpfilename,itpfilename,itpstring):
	with open(olditpfilename) as olditp:
		with open(itpfilename,"w") as newitp:
			for line in olditp:
				if line[:15] == "[ constraints ]":
					print >> newitp, itpstring
				print >> newitp, line,
	print itpstring

def basic_workflow(lower_cutoff,upper_cutoff,forceconstant,pdbfilename,itpfilename):
	[ca_indices,positions] = get_positions(pdbfilename)
	contacts = kirchoff_generation(positions,upper_cutoff,lower_cutoff)
	itpstring = construct_gnm(contacts, ca_indices, forceconstant)
	olditpfilename = backup_file(itpfilename)
	itp_append(olditpfilename,itpfilename,itpstring)
	print "finished"
	
if __name__ == "__main__":
	basic_workflow(0,7,1000,"protein_cg.pdb","protein-cg.itp")