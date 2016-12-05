from yt.funcs import issue_deprecation_warning

issue_deprecation_warning("Particle trajectories are now available from DatasetSeries "
                          "objects as ts.particle_trajectories. The ParticleTrajectories "
                          "analysis module is deprecated.")
from yt.data_objects.particle_trajectories import ParticleTrajectories