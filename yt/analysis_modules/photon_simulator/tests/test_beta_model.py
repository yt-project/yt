"""
Answer test the photon_simulator analysis module.
"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.analysis_modules.photon_simulator.api import \
    XSpecThermalModel, XSpecAbsorbModel, \
    ThermalPhotonModel, PhotonList
from yt.config import ytcfg
from yt.testing import requires_file
from yt.utilities.answer_testing.framework import \
    requires_module
import numpy as np
from yt.utilities.physical_ratios import \
    K_per_keV, mass_hydrogen_grams
from yt.frontends.stream.api import load_uniform_grid

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

test_dir = ytcfg.get("yt", "test_data_dir")

ETC = test_dir+"/enzo_tiny_cosmology/DD0046/DD0046"
ARF = test_dir+"/xray_data/sxt-s_120210_ts02um_intallpxl.arf"
RMF = test_dir+"/xray_data/ah_sxs_5ev_basefilt_20100712.rmf"

@requires_module("xspec")
@requires_file(ARF)
@requires_file(RMF)
def test_beta_model():

    R = 1.0
    r_c = 0.05
    rho_c = 0.04*mass_hydrogen_grams
    beta = 1.
    T = 6.0
    v_shift = 4.0e7
    nx = 128

    ddims = (nx,nx,nx)

    x, y, z = np.mgrid[-R:R:nx*1j,
                       -R:R:nx*1j,
                       -R:R:nx*1j]

    r = np.sqrt(x**2+y**2+z**2)

    dens = np.zeros(ddims)
    dens[r <= R] = rho_c*(1.+(r[r <= R]/r_c)**2)**(-1.5*beta)
    dens[r > R] = 0.0
    temp = T*K_per_keV*np.ones(ddims)
    bbox = np.array([[-0.5,0.5],[-0.5,0.5],[-0.5,0.5]])

    data = {}
    data["density"] = (dens, "g/cm**3")
    data["temperature"] = (temp, "K")
    data["velocity_x"] = (np.zeros(ddims), "cm/s")
    data["velocity_y"] = (np.zeros(ddims), "cm/s")
    data["velocity_z"] = (v_shift*np.ones(ddims), "cm/s")

    ds = load_uniform_grid(data, ddims, length_unit=(2*R, "Mpc"),
                           nprocs=64, bbox=bbox)

    A = 3000.
    exp_time = 1.0e5
    redshift = 0.05

    apec_model = XSpecThermalModel("bapec", 0.1, 11.5, 40000,
                                   thermal_broad=True)
    abs_model = XSpecAbsorbModel("TBabs", 0.02)

    sphere = ds.sphere("c", (0.5, "Mpc"))

    thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3)
    photons = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                      thermal_model)

    events = photons.project_photons("z", responses=[ARF,RMF],
                                     absorb_model=abs_model)

