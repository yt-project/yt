"""
Light cone example.
"""

from yt.extensions.lightcone import *

q = LightCone("128Mpc256grid_SFFB.param","lightcone.par")

q.CalculateLightConeSolution()

# If random seed was not provided in the parameter file, it can be given 
# straight to the routine.
q.CalculateLightConeSolution(seed=123456789)

# Save a text file detailing the light cone solution.
q.SaveLightConeSolution()

# Make a density light cone.
# The plot collection is returned so the final image can be
# customized and remade.
# Save the data to an hdf5 file and save images of the individual slices.
pc = q.ProjectLightCone('Density',save_stack=True,save_slice_images=True)

# Make a weighted light cone projection and save the stack to an hdf5 file.
pc = q.ProjectLightCone('Temperature',weight_field='Density',save_stack=True)

# Save the temperature stack to a different file.
q.SaveLightConeStack(file='light_cone_temperature.h5')

# Recycle current light cone solution by creating a new solution 
# that only randomizes the lateral shifts.
# This will allow the projection objects that have already been made 
# to be re-used.
# Just don't use the same random seed as the original.
q.RerandomizeLightConeSolution(987654321,recycle=True)

# Save the recycled solution.
q.SaveLightConeSolution(file='light_cone_recycled.out')

# Change the file prefix so that new light cones will not over-write the old ones.
q.lightConeParameters['OutputPrefix'] = "LightCone_NewSeed"

# Make new projection with the recycled solution.
pc = q.ProjectLightCone('Density')

# Rerandomize the light cone solution with an entirely new solution.
q.RerandomizeLightConeSolution(8675309,recycle=False)
