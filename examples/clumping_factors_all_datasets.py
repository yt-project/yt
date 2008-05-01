# This is a wrapper for clumping_factor.py, allowing you to calculate 
# clumping factor for multiple datasets.
# In this example, we calculate the standard density clumping factor 
# (<rho^2> / <rho>^2), with both mass-weighted and volume-weighted 
# average densities.
# We also get the redshift of each dump and sort the datasets by redshift 
# for output into a file.

from clumping_factor import *

data_dir = ""

# Time datadumps.
time_dumps = [252]
time_dump_dir = "DD"
time_dump_prefix = "DD"

# Redshift datadumps.
redshift_dumps = [(q+8) for q in range(68-8)]
redshift_dump_dir = "RD"
redshift_dump_prefix = "RD"

datasets = []
redshifts = []

c_m = {}
c_v = {}

# Prepare list of datasets.

for q in time_dumps:
    datasets.append("%s%s%04d/%s%04d" % (data_dir,time_dump_dir,q,time_dump_prefix,q))
for q in redshift_dumps:
    datasets.append("%s%s%04d/%s%04d" % (data_dir,redshift_dump_dir,q,redshift_dump_prefix,q))

for dataset in datasets:
    dataset_object = lagos.EnzoStaticOutput(dataset)

    redshifts.append(dataset_object.parameters["CosmologyCurrentRedshift"])

    print "Calculating clumping factors for %s (z = %f)." % (dataset,redshifts[-1])

    # Store clumping factors in dicts with the redshifts as the keys.
    c_m[redshifts[-1]] = clumping_factor(dataset_object,["Density"],[2],"CellMassMsun")
    c_v[redshifts[-1]] = clumping_factor(dataset_object,["Density"],[2],"CellVolume")

    del dataset_object

# Sort list of redshifts.
redshifts.sort()

file = open('clumping_factors.txt','w')
file.write('#z\tc_m\tc_v\n')
for z in redshifts:
    file.write("%.6e\t%.6e\t%.6e\n" % (z,c_m[z],c_v[z]))
file.close()
