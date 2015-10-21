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
from yt.utilities.answer_testing.framework import \
    requires_module
import numpy as np
from yt.utilities.physical_ratios import \
    K_per_keV, mass_hydrogen_grams
from yt.frontends.stream.api import load_uniform_grid
import os
import tempfile
import shutil

def setup():
    from yt.config import ytcfg
    ytcfg["yt", "__withintesting"] = "True"

xray_data_dir = ytcfg.get("yt", "xray_data_dir")

rmfs = ["pn-med.rmf", "acisi_aimpt_cy17.rmf",
        "aciss_aimpt_cy17.rmf", "nustar.rmf",
        "ah_sxs_5ev_basefilt_20100712.rmf"]
arfs = ["pn-med.arf", "acisi_aimpt_cy17.arf",
        "aciss_aimpt_cy17.arf", "nustar_3arcminA.arf",
        "sxt-s_120210_ts02um_intallpxl.arf"]

@requires_module("xspec")
def test_beta_model():
    import xspec
    
    xspec.Fit.statMethod = "cstat"
    xspec.Xset.addModelString("APECTHERMAL","yes")
    xspec.Fit.query = "yes"
    xspec.Fit.method = ["leven","10","0.01"]
    xspec.Fit.delta = 0.01
    xspec.Xset.chatter = 5

    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

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

    apec_model = XSpecThermalModel("bapec", 0.1, 11.5, 20000,
                                   thermal_broad=True)
    abs_model = XSpecAbsorbModel("TBabs", 0.02)

    sphere = ds.sphere("c", (0.5, "Mpc"))

    thermal_model = ThermalPhotonModel(apec_model, Zmet=0.3)
    photons = PhotonList.from_scratch(sphere, redshift, A, exp_time,
                                      thermal_model)

    for a, r in zip(rmfs, arfs):
        arf = xray_data_dir+a
        rmf = xray_data_dir+r
        events = photons.project_photons([1.0,-1.0,0.5], responses=[arf,rmf],
                                         absorb_model=abs_model)
        events.write_spectrum("beta_model_evt.pi", clobber=True)

        s = xspec.Spectrum("beta_model_evt.pi")
        s.ignore("**-0.5")
        s.ignore("7.0-**")
        m = xspec.Model("tbabs*bapec")
        m.bapec.kT = 5.0
        m.bapec.Abundanc = 0.25
        m.bapec.norm = 1.0
        m.bapec.Redshift = 0.05
        m.bapec.Velocity = 100.0
        m.TBabs.nH = 0.015

        m.bapec.Velocity.frozen = False
        m.bapec.Abundanc.frozen = False
        m.bapec.Redshift.frozen = False
        m.TBabs.nH.frozen = False

        xspec.Fit.renorm()
        xspec.Fit.nIterations = 100
        xspec.Fit.perform()

        assert np.abs(m.bapec.Redshift.values[0]-v_shift) < 1.0
        assert np.abs(m.bapec.kT.values[0]-6.0) < 1.0
        assert np.abs(m.bapec.Abundanc.values[0]-0.3) < 1.0
        assert np.abs(m.bapec.Velocity.values[0]-0.0) < 1.0
        assert np.abs(m.TBabs.nH.values[0]-0.02) < 1.0

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
