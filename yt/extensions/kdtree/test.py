from Forthon import *
from fKDpy import *
import numpy,random

n = 32768


fKD.tags = fzeros((64),'i')
fKD.dist = fzeros((64),'d')
fKD.pos = fzeros((3,n),'d')
fKD.nn = 64
fKD.nparts = n
fKD.sort = True
fKD.rearrange = True
fKD.qv = numpy.array([16./32, 16./32, 16./32])

fp = open('parts.txt','r')
xpos = []
ypos = []
zpos = []
line = fp.readline()
while line:
    line = line.split()
    xpos.append(float(line[0]))
    ypos.append(float(line[1]))
    zpos.append(float(line[2]))
    line= fp.readline()

fp.close()


for k in range(32):
    for j in range(32):
        for i in range(32):
            fKD.pos[0][i + j*32 + k*1024] = float(i)/32 + 1./64 + 0.0001*random.random()
            fKD.pos[1][i + j*32 + k*1024] = float(j)/32 + 1./64 + 0.0001*random.random()
            fKD.pos[2][i + j*32 + k*1024] = float(k)/32 + 1./64 + 0.0001*random.random()

            

#print fKD.pos[0][0],fKD.pos[1][0],fKD.pos[2][0]

create_tree()


find_nn_nearest_neighbors()

#print 'next'

#fKD.qv = numpy.array([0., 0., 0.])

#find_nn_nearest_neighbors()


#print (fKD.tags - 1)
#print fKD.dist

free_tree()
