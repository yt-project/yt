# This is a wrapper for find_clump_dataset.py, allowing you to run the clump finder
# over multiple datasets.
# In this example, we are looking for clumps in a sphere of radius 10 pc surrounding 
# the density maximum in each datadump.

import yt.lagos as lagos
from find_clumps_dataset import *

data_dir = "/Users/britton/EnzoRuns/runs_08/cl-2.5/"
#data_dir = "/Volumes/Turducken/EnzoRuns/runs_08/cl-5/"

# Time datadumps.
time_dumps = [34]
time_dump_dir = "DataDir"
time_dump_prefix = "DataDump"

# Redshift datadumps.
redshift_dumps = []
redshift_dump_dir = "RedshiftDir"
redshift_dump_prefix = "RedshiftDump"

field = "Density"
radius = 0.0001
units = "pc"
steps_per_dex = 4.
step = 10**(1./steps_per_dex)

minCells = 64 # not setting anything, only for file prefix

# Prepare list of datasets.
datasets = []

for q in time_dumps:
    datasets.append("%s%s%04d/%s%04d" % (data_dir,time_dump_dir,q,time_dump_prefix,q))
for q in redshift_dumps:
    datasets.append("%s%s%04d/%s%04d" % (data_dir,redshift_dump_dir,q,redshift_dump_prefix,q))

for dataset in datasets:
    print "Finding clumps in %s." % dataset

    prefix = "%s_%.1e%s_spd%d_min%d" % (dataset,radius,units,steps_per_dex,minCells)

    dataset_object = lagos.EnzoStaticOutput(dataset)

    # Look for clumps in a sphere surrounding the density maximum.
    v, c = dataset_object.h.find_max(field)
    sphere = dataset_object.h.sphere(c, radius/dataset_object[units], [field]) # cache our field

    print "Sphere is %s %s." % (radius,units)
    print "Min %s: %e, Max %s: %e." % (field,sphere.data[field].min(),
                                      field,sphere.data[field].max())

    master = find_clumps_dataset(prefix,sphere,field,step)

    del master
    del sphere
    del dataset_object
