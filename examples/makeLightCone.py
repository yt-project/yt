"""
Light cone example.
"""

from yt.lagos.lightcone import *

q = LightCone("128Mpc256grid_SFFB.param","lightcone.par")

q.CalculateLightConeSolution()

# Save a text file detailing the light cone solution.
q.SaveLightConeSolution()

# Make a density light cone.
# The plot collection is returned so the final image can be
# customized and remade.
pc = q.ProjectLightCone('Density')

# Make a weighted light cone projection.
pc = q.ProjectLightCone('Temperature',weight_field='Density')

# Save the light cone stack to an hdf5 file.
q.SaveLightConeStack()
