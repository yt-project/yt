#!/usr/bin/env python
# encoding: utf-8
"""
CeverinoReader.py

Created by Christopher on 2010-08-11.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import os.path
import ipdb
import numpy

class BracketAccessMixin:
    """Allow access to attributes via dictionary syntax."""
    def __getitem__(self, name):
        return getattr(self, name)

    def __setitem__(self, name, value):
        return setattr(self, name, value)

#For format info see:
#http://astronomy.nmsu.edu/danielcv/ForPriya/MW6/outputs/DataSet_README

def read_data(file):
	byte_size = os.path.getsize(file)
	
	#ipdb.set_trace()
		
	
	#Is it a DM file? Or is it a Gas file (same number of bytes per row)
	row_size = 4+4+3*8+3*8+8+4 # 68
	if numpy.mod(byte_size,row_size) == 0:
		n_dm = byte_size / row_size
		format = 'i4 ,i4 ,f8,f8,f8 ,f8,f8,f8 ,f8 ,i4'
		data = numpy.core.records.fromfile(file, formats=format, shape=n_dm, 
			byteorder='>', names='arr_length,id,posx,posy,posz,velx,vely,velz,'
								+'mass,arr_length2')
		if numpy.all(data['arr_length']==60): 
			return 0,data
		
	#Is it a Stars file? (Can't distinguish between S or Si unless by filename)
	row_size = 4+4+3*8+3*8+8+4+4 # 72
	if numpy.mod(byte_size,row_size) == 0:
		n_dm = byte_size / row_size
		format = 'i4 ,i4 ,f8,f8,f8 ,f8,f8,f8 ,f8,f4 ,i4'
		data = numpy.core.records.fromfile(file, formats=format, shape=n_dm, 
			byteorder='>', names='arr_length,id,posx,posy,posz,velx,vely,velz,'
							+'mass,age,arr_length2')
		if numpy.all(data['arr_length']==64):
			return 1,data
		
	#Is it a SZ = Stars with Metals info file?
	row_size = 4+4+3*8+3*8+8+4+4+4+4 # 80
	if numpy.mod(byte_size,row_size) == 0:
		n_dm = byte_size / row_size
		format = 'i4 ,i4 ,f8,f8,f8 ,f8,f8,f8 ,f8,f4,f4,f4 ,i4'
		data = numpy.core.records.fromfile(file, formats=format, shape=n_dm, 
			byteorder='>', names='arr_length,id,posx,posy,posz,velx,vely,velz,'
								+'mass,age,metals_snII,metals_snIa,arr_length2')
		if numpy.all(data['arr_length']==72):
			return 2,data
		
	#Is it a Gas file?
	row_size = 4+4+3*4+3*4+4+4+4 # 44
	if numpy.mod(byte_size,row_size) == 0:
		format = 'i4,f4, f4,f4,f4 , f4,f4,f4 ,f4,f4 ,i4'
		data=numpy.core.records.fromfile(file, formats=format, 
						byteorder='>', names='arr_length,cell_size,posx,posy,'
						+'posz,velx,vely,velz,density,'+
						'temperature,arr_length2')
		if numpy.all(data['arr_length']==36):
			return 3,data
		
	#Is it a GZ = Gas file with metals?
	row_size = 4+4+3*4+3*4+4+4+4+4+4 # 52
	if numpy.mod(byte_size,row_size) == 0:
		format = 'i4, f4, f4,f4,f4, f4,f4,f4, f4, f4, f4,f4, i4'
		data=numpy.core.records.fromfile(file, formats=format, 
						byteorder='>', names='arr_length,cell_size,posx,posy,'
						+'posz,velx,vely,velz,density,'+
						'temperature,metals_snII,metals_snIa,arr_length2')
		if numpy.all(data['arr_length']==44):
			return 4,data

def check_file(file):
	if not os.path.exists(file): return False
	
	byte_size = os.path.getsize(file)
	
	byte_multipliers = [68,72,80,44,52]
	for bm in byte_multipliers:
		if numpy.mod(byte_size,bm)==0: return True
	return False

def append_field(rec, name, arr, dtype=None):
    if dtype is None:
        dtype = arr.dtype
    arr = numpy.asarray(arr)
    newdtype = numpy.dtype(rec.dtype.descr + [(name, dtype)])
    newrec = numpy.empty(rec.shape, dtype=newdtype)
    for field in rec.dtype.fields:
        newrec[field] = rec[field]
    newrec[name] = arr
    return newrec

def remove_field(rec, rfield):
	descr = []
	for format in rec.dtype.descr:
		if format[0]==rfield: continue
		descr += [format,]
	newdtype = numpy.dtype(descr)
	newrec = numpy.empty(rec.shape, dtype=newdtype)
	for field in rec.dtype.fields:
		if field==rfield: continue
		newrec[field] = rec[field]
	return newrec

def read(file):
	type,data = read_data(file)
	data = remove_field(data,'arr_length')
	data = remove_field(data,'arr_length2')
	x,y,z = data['posx'],data['posy'],data['posz']
	data['cell_size'] = data['cell_size']/1000.0 
	# scale so that cell size and position same units
	s = data['cell_size']
	us = numpy.unique(s)
	xl,yl,zl = x-0.5*s,y-0.5*s,z-0.5*s
	xr,yr,zr = x+0.5*s,y+0.5*s,z+0.5*s
	dat=[(s==j)*idx for idx,j in enumerate(us)]
	l=numpy.sum(numpy.array(dat),axis=0)
	#le = numpy.array((xl,yl,zl))
	#re = numpy.array((xr,yr,zr))
	le = numpy.array([xl,yl,zl]).transpose()
	dle = numpy.min(le,axis=0)
	li = grid_indices(le,s,dle)
	lix,liy,liz = li[:,0],li[:,1],li[:,2]
	data = append_field(data,'left_edge_x',xl)
	data = append_field(data,'left_edge_y',yl)
	data = append_field(data,'left_edge_z',zl)
	data = append_field(data,'left_index_x',lix)
	data = append_field(data,'left_index_y',liy)
	data = append_field(data,'left_index_z',liz)
	data = append_field(data,'right_edge_x',xr)
	data = append_field(data,'right_edge_y',yr)
	data = append_field(data,'right_edge_z',zr)
	data = append_field(data,'level',l)
	return data
	
def parent_cells(dat):
	#Odd indices are never on the parent's LE 
	#So find where all the indices 
	
	data = append_field(dat,'parent_id',numpy.zeros(len(dat))-1)
	data = append_field(data,'is_parent',numpy.zeros(len(dat)))
	
	#Go through each cell skipping cells with the parent field != -1
	#For every cell find the parent level and integer index
	#Find parent cell by subtracting -1 off of odd integer indices
	#    and keeping even indices the same. Level is L-1
	#Find parent cell, or create if it isn't there
	#Assign parent index to the child
	
	#iterating a list that's appended to
	xl,yl,zl = data['left_edge_x'],data['left_edge_y'],data['left_edge_z']
	dlex,dley,dlez = numpy.min(xl),numpy.min(yl),numpy.min(zl),
	
	
	l,lix   = data['level'],data['left_index_x']
	liy,liz = data['left_index_y'],data['left_index_z']
	plist  = numpy.array([])
	search = numpy.array([])
	tot= len(data)
	for j,row in enumerate(data):
		if row['parent_id'] > -1: continue
		# l is level or left, i is index, e is edge
		# p is parent, cs is cell size
		# d is domain
		
		lix,liy,liz = row['left_index_x'],row['left_index_y'],row['left_index_z']
		pil = row['level']-1
		#This may give us negative level values... we'll have to normalize later
		pcs = row['cell_size']*2.0
		pix,piy,piz = lix-numpy.mod(lix,2),liy-numpy.mod(liy,2),liz-numpy.mod(liz,2)
		pex,pey,pez = dlex+pcs*pix,dley+pcs*piy,dlez+pcs*piz 
		
		#now find the parent cell given the int index and level
		#ipdb.set_trace()
		if len(plist)>0:
			(index,) = numpy.where(l==pil and lix == pix and liy == piy and liz == piz)
		else:
			index = []
			
		if len(index) == 1: 
			#we only 1 match
			[pidx,] = index
			#print "found parent:   %d %d %d (%d cells of %d)" % (pix,piy,piz,j,tot)
		elif len(index) > 2:
			#...otherwise we have a problem with two matches
			raise Exception()
		else:
			#create the parent row since we didn't find one
			prow = numpy.zeros(1,dtype=data.dtype)
			prow['parent_id'] = len(data)
			prow['level'] = pil
			prow['cell_size'] = pcs
			prow['left_index_x'] = pix
			prow['left_index_y'] = piy
			prow['left_index_z'] = piz
			prow['left_edge_x'] = pex
			prow['left_edge_y'] = pey
			prow['left_edge_z'] = pez
			prow['is_parent']   = 1
			if len(plist)>0:
		 		plist = numpy.append(plist,prow)
		 	else:
		 		plist = prow
		 	
		 	#data = numpy.append(data,prow)
		 	# yes, append to a list inside the iterated list
		 	# shouldn't break as long as we're adding to the end
		 	# but still bad practice
		 	
		 	#recreate the search list
		 	l,lix   = plist['level'],plist['left_index_x']
			liy,liz = plist['left_index_y'],plist['left_index_z']
			#search  = numpy.array([l,lix,liy,liz]).transpose()
			tot= len(data)
			#print "created parent: %d %d %d (%d cells of %d)" % (pix,piy,piz,j,tot)
		if numpy.mod(j,1000)==0:
			print "%d cells of %d, parents created: %d" % (j,tot,len(plist))  
	return data
	
	#normalize the level values
		
def grid_indices(left_vertices,cell_sizes,domain_left):
	cell_sizes_expanded = numpy.vstack((cell_sizes,)*3).transpose()
	indices = (left_vertices-domain_left)/cell_sizes_expanded
	#check that floating point error are less 0.001
	ints = numpy.rint(indices)
	bad = numpy.abs(ints*1.00-indices)>0.001
	if numpy.any(bad): 
		raise Exception("Floating point error tolerance exceeded")
	return ints
	
def grid_indices_CRAP(left_vertices,left_cell_sizes,dle,max,maxlevel=10):
	#Create one list of 1D nodes along x- y- z-axes
	#Not necessary for each dimension; cells are all cubic
	#for each of the levels
	#Starting with the finest level, move up til you don't match
	# any vertices
	#Remember the index-matches in each dim
	#Make sure all the matched occur at the samedimension
	
	#left edge vertex lists for a given level
	cell_sizes = numpy.unique(left_cell_sizes)
	levels = range(maxlevel,-1,-1) #ierate backwards
	vlistsx = [numpy.array([cell_size*j + dle[0] for j in range(round(max/cell_size))]) 
				for cell_size in cell_sizes] 
	vlistsy = [numpy.array([cell_size*j + dle[1] for j in range(round(max/cell_size))]) 
				for cell_size in cell_sizes] 
	vlistsz = [numpy.array([cell_size*j + dle[2] for j in range(round(max/cell_size))]) 
				for cell_size in cell_sizes] 
	
	j,precision = 0,0.0001
	left_edge_index = numpy.zeros(left_edges.shape)
	left_edge_level = numpy.zeros(left_edges.shape[0])
	for x,y,z in left_edges:
		for level, vlistx,vlisty,vlistz in zip(levels,vlistsx,vlistsy,vlistsz):
			(temp_x,) = numpy.where(numpy.abs(vlistx-x)<precision) 
			(temp_y,) = numpy.where(numpy.abs(vlisty-y)<precision) 
			(temp_z,) = numpy.where(numpy.abs(vlistz-z)<precision) 
			#If we have gone up to a level where our test doesn't match
			#anything quit
			if numpy.any([not len(obj)==1 for obj in (temp_x,temp_y,temp_z,)]):
				break
			else:
				cur_x,cur_y,cur_z,cur_level=temp_x[0],temp_y[0],temp_z[0],level
		left_edge_index[j] = numpy.array([cur_x,cur_y,cur_z])
		left_edge_level[j] = cur_level
		j+=1
		print "%d of %d" % (j,len(left_edges))
	return left_edge_index,left_edge_level
				
if __name__ == '__main__':
    filename='../../data/ceverino/MW6_D875.a0.941.dat'
    snap = Snapshot()
    rows = snap.read_body(filename)
    pdb.set_trace()