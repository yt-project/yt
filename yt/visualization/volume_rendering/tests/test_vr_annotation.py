from yt.testing import fake_vr_orientation_test_ds
from yt import create_scene
import numpy as np
import pytest


@pytest.mark.mpl_image_compare
def test_annotations_answer():
    from matplotlib.image import imread
    from matplotlib.pyplot import imshow, subplots


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
    sc[0].tfh.tf.add_layers(10, colormap="cubehelix")
    sc.save_annotated(
        "test_scene_annotated.png",
        text_annotate=[[(0.1, 1.05), "test_string"]],
    )
    f, ax = subplots(1)
    ax.imshow(imread("test_scene_annotated.png"))
    return f
