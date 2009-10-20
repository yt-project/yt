from yt.mods import *
import yt.lagos.Clump as cl
import numpy as na

# load the dataset
pf = load('my_dataset')

# locate the maximum density and create a 5 pc sphere around it
v, c = pf.h.find_max('Density')
sphere = pf.h.sphere(c, 5/pf['pc'], ['Density']) # cache our field

# define the min and max contour values
# use even powers of ten
c_min = 10**na.floor(na.log10(pf['Density'].min()))
c_max = 10**na.floor(na.log10(pf['Density'].max())+1)

# use increments of 0.5 dex
step = 10**(0.5)

# a string, when processed with eval will return True if a clump is bound
bound_function = 'self.data.quantities["IsBound"](truncate=True,include_thermal_energy=True) > 1.0'

# create the main clump
# The None as the second argument indicates that the main clump has no parent.
main_clump = cl.Clump(pf, None, 'Density', function=bound_function)

# begin clump-finding on the main clump
cl.find_clumps(main_clump, c_min, c_max, step)

# output the clump hierarchy, including all parent and child clumps
f = open('clump_hierarchy.txt')
cl.write_clump_hierarchy(main_clump,0,f)
f.close()

# output only the clumps without any parents
f = open('clumps.txt')
cl.write_clumps(main_clump,0,f)
f.close()
