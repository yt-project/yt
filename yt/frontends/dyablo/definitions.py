"""
Constants and definitions for the Dyablo frontend.
"""

import re

# Pattern for hydro output files: {name}_iter{NNNNNNN}.h5
HYDRO_FILE_PATTERN = re.compile(r"(.+)_iter(\d+)\.h5$")

# Pattern for particle output files: {name}_particles_{TYPE}_iter{NNNNNNN}.h5
PARTICLE_FILE_PATTERN = re.compile(r"(.+)_particles_(.+)_iter(\d+)\.h5$")
