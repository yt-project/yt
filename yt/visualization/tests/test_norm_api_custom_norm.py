import matplotlib
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
from nose.plugins.attrib import attr
from packaging.version import Version

from yt.testing import ANSWER_TEST_TAG, fake_random_ds, skip_case
from yt.utilities.answer_testing.framework import GenericImageTest
from yt.visualization.api import SlicePlot

MPL_VERSION = Version(matplotlib.__version__)


@attr(ANSWER_TEST_TAG)
def test_sliceplot_custom_norm():
    if MPL_VERSION < Version("3.4"):
        skip_case(
            reason="in MPL<3.4, SymLogNorm emits a deprecation warning "
            "that cannot be easily filtered"
        )
    # don't import this at top level because it's only available since MPL 3.2
    from matplotlib.colors import TwoSlopeNorm

    norms_to_test = [
        (Normalize(), "linear"),
        (LogNorm(), "log"),
        (TwoSlopeNorm(vcenter=0, vmin=-0.5, vmax=1), "twoslope"),
        (SymLogNorm(linthresh=0.01, vmin=-1, vmax=1), "symlog"),
    ]

    ds = fake_random_ds(16)

    def create_image(filename_prefix):
        field = ("gas", "density")
        for norm, name in norms_to_test:
            p = SlicePlot(ds, "z", field)
            p.set_norm(field, norm=norm)
            p.save(f"{filename_prefix}_{name}")

    test = GenericImageTest(ds, create_image, 12)
    test.prefix = "test_sliceplot_custom_norm"
    test.answer_name = "sliceplot_custom_norm"
    yield test
