"""
RAMSES-specific fields



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import glob
import os
import numpy as np

from yt.utilities.physical_constants import \
    boltzmann_constant_cgs, \
    mass_hydrogen_cgs, \
    mh
from yt.utilities.linear_interpolators import \
    BilinearFieldInterpolator
import yt.utilities.fortran_utils as fpu
from yt.fields.field_info_container import \
    FieldInfoContainer
from yt.utilities.physical_ratios import cm_per_km, cm_per_mpc

b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_density"
vel_units = "code_velocity"
pressure_units = "code_pressure"

known_species_masses = dict(
  (sp, mh * v) for sp, v in [
                ("HI", 1.0),
                ("HII", 1.0),
                ("Electron", 1.0),
                ("HeI", 4.0),
                ("HeII", 4.0),
                ("HeIII", 4.0),
                ("H2I", 2.0),
                ("H2II", 2.0),
                ("HM", 1.0),
                ("DI", 2.0),
                ("DII", 2.0),
                ("HDI", 3.0),
    ])

_cool_axes = ("lognH", "logT", "logTeq")
_cool_arrs = ("metal", "cool", "heat", "metal_prime", "cool_prime",
              "heat_prime", "mu", "abundances")
_cool_species = ("Electron_number_density",
                 "HI_number_density",
                 "HII_number_density",
                 "HeI_number_density",
                 "HeII_number_density",
                 "HeIII_number_density")

_X = 0.76 # H fraction, hardcoded
_Y = 0.24 # He fraction, hardcoded

class RAMSESFieldInfo(FieldInfoContainer):
    known_other_fields = (
        ("Density", (rho_units, ["density"], None)),
        ("x-velocity", (vel_units, ["velocity_x"], None)),
        ("y-velocity", (vel_units, ["velocity_y"], None)),
        ("z-velocity", (vel_units, ["velocity_z"], None)),
        ("Pres_IR", (pressure_units, ["pres_IR"], None)),
        ("Pressure", (pressure_units, ["pressure"], None)),
        ("Metallicity", ("", ["metallicity"], None)),
        ("HII",  ("", ["H_p1_fraction"], None)),
        ("HeII", ("", ["He_p1_fraction"], None)),
        ("HeIII",("", ["He_p2_fraction"], None)),
    )
    known_particle_fields = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", (vel_units, [], None)),
        ("particle_velocity_y", (vel_units, [], None)),
        ("particle_velocity_z", (vel_units, [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("particle_identifier", ("", ["particle_index"], None)),
        ("particle_refinement_level", ("", [], None)),
        ("particle_age", ("", [""], None)),
        ("particle_metallicity", ("", [], None)),
    )

    def setup_fluid_fields(self):
        def _temperature(field, data):
            rv = data["gas", "pressure"]/data["gas", "density"]
            rv *= mass_hydrogen_cgs/boltzmann_constant_cgs
            return rv
        self.add_field(("gas", "temperature"), sampling_type="cell",  function=_temperature,
                        units=self.ds.unit_system["temperature"])
        self.create_cooling_fields()
        # See if we need to load the rt fields
        foldername  = os.path.abspath(os.path.dirname(self.ds.parameter_filename))
        rt_flag = any(glob.glob(os.sep.join([foldername, 'info_rt_*.txt'])))
        if rt_flag: # rt run
            self.setup_rt_fields()

    def setup_rt_fields(self):
        def _temp_IR(field, data):
            rv = data["gas", "pres_IR"]/data["gas", "density"]
            rv *= mass_hydrogen_cgs/boltzmann_constant_cgs
            return rv
        self.add_field(("gas", "temp_IR"), sampling_type="cell",
                       function=_temp_IR,
                       units=self.ds.unit_system["temperature"])
        for species in ['H_p1', 'He_p1', 'He_p2']:
            def _species_density(field, data):
                return data['gas', species+'_fraction']*data['gas', 'density']
            self.add_field(('gas', species+'_density'), sampling_type='cell',
                           function=_species_density,
                           units=self.ds.unit_system['density'])
            def _species_mass(field, data):
                return (data['gas', species+'_density']*
                        data['index', 'cell_volume'])
            self.add_field(('gas', species+'_mass'), sampling_type='cell',
                           function=_species_mass,
                           units=self.ds.unit_system['mass'])


    def setup_particle_fields(self, ptype):
        super(RAMSESFieldInfo, self).setup_particle_fields(ptype)
        self.setup_age_fields(ptype)

    def create_cooling_fields(self):
        num = os.path.basename(self.ds.parameter_filename).split("."
                )[0].split("_")[1]
        filename = "%s/cooling_%05i.out" % (
            os.path.dirname(self.ds.parameter_filename), int(num))

        if not os.path.exists(filename): return
        def _create_field(name, interp_object):
            def _func(field, data):
                shape = data["temperature"].shape
                d = {'lognH': np.log10(_X*data["density"]/mh).ravel(),
                     'logT' : np.log10(data["temperature"]).ravel()}
                rv = 10**interp_object(d).reshape(shape)
                # Return array in unit 'per volume' consistently with line below
                return data.ds.arr(rv, 'code_length**-3')
            self.add_field(name = name, sampling_type="cell", function=_func,
                                 units = "code_length**-3")
        avals = {}
        tvals = {}
        with open(filename, "rb") as f:
            n1, n2 = fpu.read_vector(f, 'i')
            n = n1 * n2
            for ax in _cool_axes:
                avals[ax] = fpu.read_vector(f, 'd')
            for tname in _cool_arrs:
                var = fpu.read_vector(f, 'd')
                if var.size == n1*n2:
                    tvals[tname] = var.reshape((n1, n2), order='F')
                else:
                    var = var.reshape((n1, n2, var.size // (n1*n2)), order='F')
                    for i in range(var.shape[-1]):
                        tvals[_cool_species[i]] = var[:,:,i]

        for n in tvals:
            interp = BilinearFieldInterpolator(tvals[n],
                        (avals["lognH"], avals["logT"]),
                        ["lognH", "logT"], truncate = True)
            _create_field(("gas", n), interp)

    def setup_age_fields(self, ptype):
        ds = self.ds
        if (ptype, 'particle_age') in self:
            if ds.cosmological_simulation:
                def _age(field, data):
                    tf = ds.t_frw
                    dtau = ds.dtau
                    tauf = ds.tau_frw
                    tsim = ds.time_simu
                    h100 = ds.hubble_constant
                    nOver2 = ds.n_frw/2
                    t_scale = 1./(h100 * 100 * cm_per_km / cm_per_mpc)
                    t_scale /= ds['unit_t']
                    ages = data[ptype, 'particle_age']
                    dage = 1 + (10*ages/dtau)
                    dage = np.minimum(dage, nOver2 + (dage - nOver2)/10.)
                    iage = np.array(dage,dtype=np.int32)

                    t = ((tf[iage]*(ages - tauf[iage - 1]) /
                          (tauf[iage] - tauf[iage - 1])))
                    t = t + ((tf[iage-1]*(ages-tauf[iage]) /
                              (tauf[iage-1]-tauf[iage])))
                    return ds.arr((tsim - t)*t_scale, 'code_time')
                self.add_field((ptype, 'age'), sampling_type='particle',
                               function=_age, units=self.ds.unit_system['time'])
            else:
                self[ptype, 'particle_age'].units = 'code_time'
                self.alias((ptype, 'age'), (ptype, 'particle_age'))
