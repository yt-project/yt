"""
Unit test the PPVCube analysis module.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.frontends.stream.api import load_uniform_grid
from yt.analysis_modules.ppv_cube.api import PPVCube
import yt.units as u
from yt.utilities.physical_constants import kboltz, mh, clight
import numpy as np
from yt.testing import assert_allclose_units

def setup():
    """Test specific setup."""
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

def test_ppv():

    np.random.seed(seed=0x4d3d3d3)

    dims = (8, 8, 128)
    v_shift = 1.0e7*u.cm/u.s
    sigma_v = 2.0e7*u.cm/u.s
    T_0 = 1.0e8*u.Kelvin
    data = {"density":(np.ones(dims),"g/cm**3"),
            "temperature":(T_0.v*np.ones(dims), "K"),
            "velocity_x":(np.zeros(dims),"cm/s"),
            "velocity_y":(np.zeros(dims),"cm/s"),
            "velocity_z":(np.random.normal(loc=v_shift.v,scale=sigma_v.v,size=dims), "cm/s")}

    ds = load_uniform_grid(data, dims)

    cube = PPVCube(ds, "z", "density", (-300., 300., 1024, "km/s"),
                   dims=8, thermal_broad=True)

    dv = cube.dv
    v_th = np.sqrt(2.*kboltz*T_0/(56.*mh) + 2.*sigma_v**2).in_units("km/s")
    a = cube.data.mean(axis=(0,1)).v
    b = dv*np.exp(-((cube.vmid+v_shift)/v_th)**2)/(np.sqrt(np.pi)*v_th)

    assert_allclose_units(a, b, 1.0e-2)

    E_0 = 6.8*u.keV

    cube.transform_spectral_axis(E_0.v, str(E_0.units))

    dE = -cube.dv
    delta_E = E_0*v_th.in_cgs()/clight
    E_shift = E_0*(1.+v_shift/clight)

    c = dE*np.exp(-((cube.vmid-E_shift)/delta_E)**2)/(np.sqrt(np.pi)*delta_E)

    assert_allclose_units(a, c, 1.0e-2)

def test_ppv_nothermalbroad():

    np.random.seed(seed=0x4d3d3d3)

    dims = (16, 16, 128)
    v_shift = 1.0e6*u.cm/u.s
    sigma_v = 2.0e6*u.cm/u.s
    data = {"density":(np.ones(dims),"g/cm**3"),
            "velocity_x":(np.zeros(dims),"cm/s"),
            "velocity_y":(np.zeros(dims),"cm/s"),
            "velocity_z":(np.random.normal(loc=v_shift.v,scale=sigma_v.v,size=dims), "cm/s")}

    ds = load_uniform_grid(data, dims)

    cube = PPVCube(ds, "z", "density", (-100., 100., 128, "km/s"),
                   dims=16, thermal_broad=False)

    dv = cube.dv
    v_noth = np.sqrt(2)*(sigma_v).in_units("km/s")
    a = cube.data.mean(axis=(0,1)).v
    b = dv*np.exp(-((cube.vmid+v_shift)/v_noth)**2)/(np.sqrt(np.pi)*v_noth)

    assert_allclose_units(a, b, atol=5.0e-3)
