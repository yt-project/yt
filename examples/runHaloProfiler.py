from yt.extensions import HaloProfiler

# Instantiate HaloProfiler for this dataset and load parameters from sample_halo_profiles.par.
q = HaloProfiler("/Users/britton/EnzoRuns/cc_run/DD0232/DD0232","sample_halo_profiles.par")

# Make profiles of all halos in dataset.
q.makeProfiles()

# Make projections of all halos in dataset.
q.makeProjections()
