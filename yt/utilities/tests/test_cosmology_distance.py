from numpy.testing import assert_almost_equal

from yt.utilities.cosmology import Cosmology


def test_observer_redshift_distance_duality():
    """
    Test angular diameter distance and luminosity distance with observer redshift
    zi > 0 following Hogg (1999) and the distance duality (Etherington) relation.
    """
    co = Cosmology()
    for zi, zf in [(0.0, 1.0), (0.0, 2.0), (0.5, 2.0), (1.0, 3.0)]:
        dm = co.comoving_transverse_distance(zi, zf)
        da = co.angular_diameter_distance(zi, zf)
        dl = co.luminosity_distance(zi, zf)

        # Hogg (1999) eq. 18: D_A = D_M / (1 + z_f)
        expected_da = dm / (1.0 + zf)
        assert_almost_equal(da.to("Mpc").v, expected_da.to("Mpc").v, decimal=10)

        # Etherington reciprocity / distance duality: D_L = ((1 + z_f)/(1 + z_i))^2 * D_A
        expected_dl = da * ((1.0 + zf) / (1.0 + zi)) ** 2
        assert_almost_equal(dl.to("Mpc").v, expected_dl.to("Mpc").v, decimal=10)

        # Angular diameter distance must be strictly positive for z_f > z_i
        assert da.v > 0
        assert dl.v > 0
