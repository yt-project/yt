import numpy as np

import yt
from yt.data_objects.particle_filters import add_particle_filter


# Define filter functions for our particle filters based on stellar age.
# In this dataset particles in the initial conditions are given creation
# times arbitrarily far into the future, so stars with negative ages belong
# in the old stars filter.
def stars_10Myr(pfilter, data):
    age = data.ds.current_time - data["Stars", "creation_time"]
    filter = np.logical_and(age >= 0, age.in_units("Myr") < 10)
    return filter


def stars_100Myr(pfilter, data):
    age = (data.ds.current_time - data["Stars", "creation_time"]).in_units("Myr")
    filter = np.logical_and(age >= 10, age < 100)
    return filter


def stars_old(pfilter, data):
    age = data.ds.current_time - data["Stars", "creation_time"]
    filter = np.logical_or(age < 0, age.in_units("Myr") >= 100)
    return filter


# Create the particle filters
add_particle_filter(
    "stars_young",
    function=stars_10Myr,
    filtered_type="Stars",
    requires=["creation_time"],
)
add_particle_filter(
    "stars_medium",
    function=stars_100Myr,
    filtered_type="Stars",
    requires=["creation_time"],
)
add_particle_filter(
    "stars_old", function=stars_old, filtered_type="Stars", requires=["creation_time"]
)

# Load a dataset and apply the particle filters
filename = "TipsyGalaxy/galaxy.00300"
ds = yt.load(filename)
ds.add_particle_filter("stars_young")
ds.add_particle_filter("stars_medium")
ds.add_particle_filter("stars_old")

# What are the total masses of different ages of star in the whole simulation
# volume?
ad = ds.all_data()
mass_young = ad["stars_young", "particle_mass"].in_units("Msun").sum()
mass_medium = ad["stars_medium", "particle_mass"].in_units("Msun").sum()
mass_old = ad["stars_old", "particle_mass"].in_units("Msun").sum()
print(f"Mass of young stars = {mass_young:g} Msun")
print(f"Mass of medium stars = {mass_medium:g} Msun")
print(f"Mass of old stars = {mass_old:g} Msun")

# Generate 4 projections: gas density, young stars, medium stars, old stars
fields = [
    ("stars_young", "particle_mass"),
    ("stars_medium", "particle_mass"),
    ("stars_old", "particle_mass"),
]

prj1 = yt.ProjectionPlot(ds, "z", ("gas", "density"), center="max", width=(100, "kpc"))
prj1.save()
prj2 = yt.ParticleProjectionPlot(ds, "z", fields, center="max", width=(100, "kpc"))
prj2.save()
