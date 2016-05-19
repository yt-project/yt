import os
import yt
from yt.analysis_modules.cosmological_observation.api import \
    LightRay

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
lr = LightRay(ds)

# With a single dataset, a start_position and
# end_position or trajectory must be given.
# Trajectory should be given as (r, theta, phi)
lr.make_light_ray(start_position=[0., 0., 0.],
                  end_position=[1., 1., 1.],
                  solution_filename='lightraysolution.txt',
                  data_filename='lightray.h5',
                  fields=['temperature', 'density'])

# Optionally, we can now overplot this ray on a projection of the source
# dataset
p = yt.ProjectionPlot(ds, 'z', 'density')
p.annotate_ray(lr)
p.save()
