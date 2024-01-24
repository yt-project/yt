import numpy as np

from yt import create_scene
from yt.testing import fake_vr_orientation_test_ds


def test_annotations(tmp_path):
    from matplotlib.image import imread

    tmp_file = tmp_path / "test_scene_annotated.png"
    ds = fake_vr_orientation_test_ds(N=16)
    sc = create_scene(ds)
    sc.annotate_axes()
    sc.annotate_domain(ds)
    sc.render()
    # ensure that there are actually red, green, blue, and white pixels
    # in the image. see Issue #1595
    im = sc._last_render
    for c in ([1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 1], [1, 1, 1, 1]):
        assert np.where((im == c).all(axis=-1))[0].shape[0] > 0
    assert im.shape == sc.camera.resolution + (4,)
    sc[0].tfh.tf.add_layers(10, colormap="cubehelix")
    sc.save_annotated(
        tmp_file,
        text_annotate=[[(0.1, 1.05), "test_string"]],
    )
    image = imread(tmp_file)
    assert image.shape == sc.camera.resolution + (4,)
