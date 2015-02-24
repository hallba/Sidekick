#!/usr/bin/python

class base:
	def __iter__(self):
		return self
	def __enter__(self):
		return self
	def __exit__(self, type, value, traceback):
		self.close()
	def close(self):
		self.source.close()
		
class read(base):
	def __init__(self,filename):
		self.source = open(filename, "r")
		self.frame=0
	def next(self):
		self.frame, xshape, yshape = [int(item) for item in self.source.next().split()]
		matrix = [[float(item) for item in self.source.next().split()] for row in range(yshape)]
		length_list = [len(item) for item in matrix]
		if min(length_list) != xshape or max(length_list) != xshape:
			raise
		self.matrix = matrix
		return matrix
	
class write(base):
	def __init__(self,filename):
		self.source = open(filename, "w")
	def store(self,matrix,identifier=0):
		xshape, yshape = len(matrix[0]), len(matrix)
		print >> self.source, identifier, xshape, yshape
		for row in matrix:
			for column in row:
				print >> self.source, "% 08.3f" % column,
			print >> self.source, ""
		
class append(write):
	def __init__(self,filename):
		self.source = open(filename, "a")

class read_3D(read):
	def next(self):
		self.frame, xshape, yshape, zshape = [int(item) for item in self.source.next().split()]
		matrix = [[[float(value) for value in item.split(',')] for item in self.source.next().split()] for row in range(yshape)]
		length_list = [len(item) for item in matrix]
		if min(length_list) != xshape or max(length_list) != xshape:
			raise
		self.matrix = matrix
		return matrix

class write_3D(write):
	def store(self,matrix,identifier=0):
		xshape, yshape, zshape = len(matrix[0]), len(matrix), len(matrix[0][0])
		print >> self.source, identifier, xshape, yshape, zshape
		for row in matrix:
			for column in row:
				stack = ""
				for item in column:
					stack += "% 08.3f," % item
				print >> self.source, stack[:-1],
			print >> self.source, ""

class append_3D(write_3D):
	def __init__(self,filename):
		self.source = open(filename, "a")
