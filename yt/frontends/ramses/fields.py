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
        ("Pressure", (pressure_units, ["pressure"], None)),
        ("Metallicity", ("", ["metallicity"], None)),
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
        ("particle_age", ("code_time", ['age'], None)),
        ("particle_metallicity", ("", [], None)),
    )

    def setup_fluid_fields(self):
        def _temperature(field, data):
            rv = data["gas", "pressure"]/data["gas", "density"]
            rv *= mass_hydrogen_cgs/boltzmann_constant_cgs
            return rv
        self.add_field(("gas", "temperature"), function=_temperature,
                        units=self.ds.unit_system["temperature"])
        self.create_cooling_fields()

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
                return rv
            self.add_field(name = name, function=_func,
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
                    var = var.reshape((n1, n2, var.size / (n1*n2)), order='F')
                    for i in range(var.shape[-1]):
                        tvals[_cool_species[i]] = var[:,:,i]
        
        for n in tvals:
            interp = BilinearFieldInterpolator(tvals[n],
                        (avals["lognH"], avals["logT"]),
                        ["lognH", "logT"], truncate = True)
            _create_field(("gas", n), interp)
