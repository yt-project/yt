# This is a wrapper for find_clump_dataset.py, allowing you to run the clump finder
# over multiple datasets.
# In this example, we are looking for clumps in a sphere of radius 10 pc surrounding 
# the density maximum in each datadump.

from find_clumps_dataset import *

data_dir = "/Users/britton/EnzoRuns/real_4/cl-4/"

# Time datadumps.
time_dumps = [q for q in range(25)]
time_dump_dir = "DataDir"
time_dump_prefix = "DataDump"

# Redshift datadumps.
redshift_dumps = [q for q in range(5)]
redshift_dump_dir = "RedshiftDir"
redshift_dump_prefix = "RedshiftDump"

field = "Density"
radius = 10
units = "pc"
step = 10**(1./5.) # 5 steps each order of magnitude

# Prepare list of datasets.
datasets = []

for q in time_dumps:
    datasets.append("%s%s%04d/%s%04d" % (data_dir,time_dump_dir,q,time_dump_prefix,q))
for q in redshift_dumps:
    datasets.append("%s%s%04d/%s%04d" % (data_dir,redshift_dump_dir,q,redshift_dump_prefix,q))

for dataset in datasets:
    print "Finding clumps in %s." % dataset

    dataset_object = lagos.EnzoStaticOutput(dataset)

    # Look for clumps in a sphere surrounding the density maximum.
    v, c = dataset_object.h.find_max(field)
    sphere = dataset_object.h.sphere(c, radius/dataset_object[units], [field]) # cache our field

    find_clumps_dataset(dataset,sphere,field,step)

    del sphere
    del dataset_object
