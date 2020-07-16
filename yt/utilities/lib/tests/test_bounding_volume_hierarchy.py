import numpy as np

import yt
from yt.testing import requires_file
from yt.utilities.lib.bounding_volume_hierarchy import BVH, test_ray_trace
from yt.visualization.volume_rendering.api import Camera, Scene


def get_rays(camera):

    normal_vector = camera.unit_vectors[2].d
    W = np.array([8.0, 8.0])
    N = np.array([800, 800])
    dx = W / N

    x_points = np.linspace((-N[0] / 2 + 0.5) * dx[0], (N[0] / 2 - 0.5) * dx[0], N[0])
    y_points = np.linspace((-N[1] / 2 + 0.5) * dx[1], (N[1] / 2 - 0.5) * dx[0], N[1])

    X, Y = np.meshgrid(x_points, y_points)
    result = np.dot(camera.unit_vectors[0:2].T, [X.ravel(), Y.ravel()])
    vec_origins = camera.scene.arr(result.T, "unitary") + camera.position
    return np.array(vec_origins), np.array(normal_vector)


fn = "MOOSE_sample_data/out.e-s010"


@requires_file(fn)
def test_bounding_volume_hierarchy():
    ds = yt.load(fn)
    vertices = ds.index.meshes[0].connectivity_coords
    indices = ds.index.meshes[0].connectivity_indices - 1

    ad = ds.all_data()
    field_data = ad["connect1", "diffused"]

    bvh = BVH(vertices, indices, field_data)

    sc = Scene()
    cam = Camera(sc)
    cam.set_position(np.array([8.0, 8.0, 8.0]))
    cam.focus = np.array([0.0, 0.0, 0.0])
    origins, direction = get_rays(cam)

    image = np.empty(800 * 800, np.float64)
    test_ray_trace(image, origins, direction, bvh)
    image = image.reshape((800, 800))
    return image
