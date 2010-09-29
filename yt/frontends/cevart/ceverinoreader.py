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
			byteorder='>', names='arr_length,id,posx,posy,posz,velx,vely,velz,mass,arr_length2')
		if numpy.all(data['arr_length']==60): 
			return 0,data
		
	#Is it a Stars file? (Can't distinguish between S or Si unless by filename)
	row_size = 4+4+3*8+3*8+8+4+4 # 72
	if numpy.mod(byte_size,row_size) == 0:
		n_dm = byte_size / row_size
		format = 'i4 ,i4 ,f8,f8,f8 ,f8,f8,f8 ,f8,f4 ,i4'
		data = numpy.core.records.fromfile(file, formats=format, shape=n_dm, 
			byteorder='>', names='arr_length,id,posx,posy,posz,velx,vely,velz,mass,age,arr_length2')
		if numpy.all(data['arr_length']==64):
			return 1,data
		
	#Is it a SZ = Stars with Metals info file?
	row_size = 4+4+3*8+3*8+8+4+4+4+4 # 80
	if numpy.mod(byte_size,row_size) == 0:
		n_dm = byte_size / row_size
		format = 'i4 ,i4 ,f8,f8,f8 ,f8,f8,f8 ,f8,f4,f4,f4 ,i4'
		data = numpy.core.records.fromfile(file, formats=format, shape=n_dm, 
			byteorder='>', names='arr_length,id,posx,posy,posz,velx,vely,velz,mass,age,'+
								'metals_snII,metals_snIa,arr_length2')
		if numpy.all(data['arr_length']==72):
			return 2,data
		
	#Is it a Gas file?
	row_size = 4+4+3*4+3*4+4+4+4 # 44
	if numpy.mod(byte_size,row_size) == 0:
		format = 'i4,f4, f4,f4,f4 , f4,f4,f4 ,f4,f4 ,i4'
		data=numpy.core.records.fromfile(file, formats=format, 
						byteorder='>', names='arr_length,cell_size,posx,posy,posz,velx,vely,velz,density,'+
						'temperature,arr_length2')
		if numpy.all(data['arr_length']==36):
			return 3,data
		
	#Is it a GZ = Gas file with metals?
	row_size = 4+4+3*4+3*4+4+4+4+4+4 # 52
	if numpy.mod(byte_size,row_size) == 0:
		format = 'i4, f4, f4,f4,f4, f4,f4,f4, f4, f4, f4,f4, i4'
		data=numpy.core.records.fromfile(file, formats=format, 
						byteorder='>', names='arr_length,cell_size,posx,posy,posz,velx,vely,velz,density,'+
						'temperature,metals_snII,metals_snIa,arr_length2')
		if numpy.all(data['arr_length']==44):
			return 4,data

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
	s = data['cell_size']/1000
	us = numpy.unique(s)
	xl,yl,zl = x-0.5*s,y-0.5*s,z-0.5*s
	xr,yr,zr = x+0.5*s,y+0.5*s,z+0.5*s
	dat=[(s==j)*idx for idx,j in enumerate(us)]
	l=numpy.sum(numpy.array(dat),axis=0)
	#le = numpy.array((xl,yl,zl))
	#re = numpy.array((xr,yr,zr))
	data = append_field(data,'left_edge_x',xl)
	data = append_field(data,'left_edge_y',xl)
	data = append_field(data,'left_edge_z',xl)
	data = append_field(data,'right_edge_x',xl)
	data = append_field(data,'right_edge_y',xl)
	data = append_field(data,'right_edge_z',xl)
	data = append_field(data,'level',l)
	return data
		
if __name__ == '__main__':
    filename='../../data/ceverino/MW6_D875.a0.941.dat'
    snap = Snapshot()
    rows = snap.read_body(filename)
    pdb.set_trace()