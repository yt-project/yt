import os
from typing import Tuple

import numpy as np

from yt.config import ytcfg
from yt.fields.species_fields import (
    add_species_field_by_density,
    add_species_field_by_fraction,
)
from yt.frontends.sph.fields import SPHFieldInfo
from yt.funcs import download_file, mylog

from . import owls_ion_tables as oit


def _get_ion_mass_frac(ion, ftype, itab, data):
    # get element symbol from ion string. ion string will
    # be a member of the tuple _ions (i.e. si13)
    # --------------------------------------------------------
    if ion[0:2].isalpha():
        symbol = ion[0:2].capitalize()
    else:
        symbol = ion[0:1].capitalize()

    # mass fraction for the element
    # --------------------------------------------------------
    m_frac = data[ftype, symbol + "_fraction"]

    # get nH and T for lookup
    # --------------------------------------------------------
    log_nH = np.log10(data["PartType0", "H_number_density"])
    log_T = np.log10(data["PartType0", "Temperature"])

    # get name of owls_ion_file for given ion
    # --------------------------------------------------------
    itab.set_iz(data.ds.current_redshift)

    # find ion balance using log nH and log T
    # --------------------------------------------------------
    i_frac = itab.interp(log_nH, log_T)

    return i_frac, m_frac


class OWLSFieldInfo(SPHFieldInfo):

    _ions: Tuple[str, ...] = (
        "c1",
        "c2",
        "c3",
        "c4",
        "c5",
        "c6",
        "fe2",
        "fe17",
        "h1",
        "he1",
        "he2",
        "mg1",
        "mg2",
        "n2",
        "n3",
        "n4",
        "n5",
        "n6",
        "n7",
        "ne8",
        "ne9",
        "ne10",
        "o1",
        "o6",
        "o7",
        "o8",
        "si2",
        "si3",
        "si4",
        "si13",
    )

    _elements = ("H", "He", "C", "N", "O", "Ne", "Mg", "Si", "Fe")

    _num_neighbors = 48

    _add_elements = ("PartType0", "PartType4")

    _add_ions = "PartType0"

    def __init__(self, ds, field_list, slice_info=None):

        new_particle_fields = (
            ("Hydrogen", ("", ["H_fraction"], None)),
            ("Helium", ("", ["He_fraction"], None)),
            ("Carbon", ("", ["C_fraction"], None)),
            ("Nitrogen", ("", ["N_fraction"], None)),
            ("Oxygen", ("", ["O_fraction"], None)),
            ("Neon", ("", ["Ne_fraction"], None)),
            ("Magnesium", ("", ["Mg_fraction"], None)),
            ("Silicon", ("", ["Si_fraction"], None)),
            ("Iron", ("", ["Fe_fraction"], None)),
        )

        if ds.gen_hsmls:
            new_particle_fields += (("smoothing_length", ("code_length", [], None)),)

        self.known_particle_fields += new_particle_fields

        super().__init__(ds, field_list, slice_info=slice_info)

        # This enables the machinery in yt.fields.species_fields
        self.species_names += list(self._elements)

    def setup_particle_fields(self, ptype):
        """additional particle fields derived from those in snapshot.
        we also need to add the smoothed fields here b/c setup_fluid_fields
        is called before setup_particle_fields."""

        smoothed_suffixes = ("_number_density", "_density", "_mass")

        # we add particle element fields for stars and gas
        # -----------------------------------------------------
        if ptype in self._add_elements:

            # this adds the particle element fields
            # X_density, X_mass, and X_number_density
            # where X is an item of self._elements.
            # X_fraction are defined in snapshot
            # -----------------------------------------------
            for s in self._elements:
                field_names = add_species_field_by_fraction(self, ptype, s)
                if ptype == self.ds._sph_ptypes[0]:
                    for fn in field_names:
                        self.alias(("gas", fn[1]), fn)

        # this needs to be called after the call to
        # add_species_field_by_fraction for some reason ...
        # not sure why yet.
        # -------------------------------------------------------
        if ptype == "PartType0":
            ftype = "gas"
        else:
            ftype = ptype

        super().setup_particle_fields(
            ptype, num_neighbors=self._num_neighbors, ftype=ftype
        )

        # and now we add the smoothed versions for PartType0
        # -----------------------------------------------------
        if ptype == "PartType0":

            # we only add ion fields for gas.  this takes some
            # time as the ion abundances have to be interpolated
            # from cloudy tables (optically thin)
            # -----------------------------------------------------

            # this defines the ion density on particles
            # X_density for all items in self._ions
            # -----------------------------------------------
            self.setup_gas_ion_particle_fields(ptype)

            # this adds the rest of the ion particle fields
            # X_fraction, X_mass, X_number_density
            # -----------------------------------------------
            for ion in self._ions:

                # construct yt name for ion
                # ---------------------------------------------------
                if ion[0:2].isalpha():
                    symbol = ion[0:2].capitalize()
                    roman = int(ion[2:])
                else:
                    symbol = ion[0:1].capitalize()
                    roman = int(ion[1:])

                if (ptype, symbol + "_fraction") not in self.field_aliases:
                    continue

                pstr = f"_p{roman - 1}"
                yt_ion = symbol + pstr

                # add particle field
                # ---------------------------------------------------
                add_species_field_by_density(self, ptype, yt_ion)

            def _h_p1_density(field, data):
                return data[ptype, "H_density"] - data[ptype, "H_p0_density"]

            self.add_field(
                (ptype, "H_p1_density"),
                sampling_type="particle",
                function=_h_p1_density,
                units=self.ds.unit_system["density"],
            )

            add_species_field_by_density(self, ptype, "H_p1")
            for sfx in ["mass", "density", "number_density"]:
                fname = "H_p1_" + sfx
                self.alias(("gas", fname), (ptype, fname))

            def _el_number_density(field, data):
                n_e = data[ptype, "H_p1_number_density"]
                n_e += data[ptype, "He_p1_number_density"]
                n_e += 2.0 * data[ptype, "He_p2_number_density"]
                return n_e

            self.add_field(
                (ptype, "El_number_density"),
                sampling_type="particle",
                function=_el_number_density,
                units=self.ds.unit_system["number_density"],
            )
            self.alias(("gas", "El_number_density"), (ptype, "El_number_density"))

            # alias ion fields
            # -----------------------------------------------
            for ion in self._ions:

                # construct yt name for ion
                # ---------------------------------------------------
                if ion[0:2].isalpha():
                    symbol = ion[0:2].capitalize()
                    roman = int(ion[2:])
                else:
                    symbol = ion[0:1].capitalize()
                    roman = int(ion[1:])

                if (ptype, symbol + "_fraction") not in self.field_aliases:
                    continue

                pstr = f"_p{roman - 1}"
                yt_ion = symbol + pstr

                for sfx in smoothed_suffixes:
                    fname = yt_ion + sfx
                    self.alias(("gas", fname), (ptype, fname))

    def setup_gas_ion_particle_fields(self, ptype):
        """Sets up particle fields for gas ion densities."""

        # loop over all ions and make fields
        # ----------------------------------------------
        for ion in self._ions:

            # construct yt name for ion
            # ---------------------------------------------------
            if ion[0:2].isalpha():
                symbol = ion[0:2].capitalize()
                roman = int(ion[2:])
            else:
                symbol = ion[0:1].capitalize()
                roman = int(ion[1:])

            if (ptype, symbol + "_fraction") not in self.field_aliases:
                continue

            pstr = f"_p{roman - 1}"
            yt_ion = symbol + pstr
            ftype = ptype

            # add ion density and mass field for this species
            # ------------------------------------------------
            fname = yt_ion + "_density"
            dens_func = self._create_ion_density_func(ftype, ion)
            self.add_field(
                (ftype, fname),
                sampling_type="particle",
                function=dens_func,
                units=self.ds.unit_system["density"],
            )
            self._show_field_errors.append((ftype, fname))

            fname = yt_ion + "_mass"
            mass_func = self._create_ion_mass_func(ftype, ion)
            self.add_field(
                (ftype, fname),
                sampling_type="particle",
                function=mass_func,
                units=self.ds.unit_system["mass"],
            )
            self._show_field_errors.append((ftype, fname))

    def _create_ion_density_func(self, ftype, ion):
        """returns a function that calculates the ion density of a particle."""

        def get_owls_ion_density_field(ion, ftype, itab):
            def _func(field, data):
                m_frac, i_frac = _get_ion_mass_frac(ion, ftype, itab, data)
                return data[ftype, "Density"] * m_frac * i_frac

            return _func

        ion_path = self._get_owls_ion_data_dir()
        fname = os.path.join(ion_path, ion + ".hdf5")
        itab = oit.IonTableOWLS(fname)
        return get_owls_ion_density_field(ion, ftype, itab)

    def _create_ion_mass_func(self, ftype, ion):
        """returns a function that calculates the ion mass of a particle"""

        def get_owls_ion_mass_field(ion, ftype, itab):
            def _func(field, data):
                m_frac, i_frac = _get_ion_mass_frac(ion, ftype, itab, data)
                return data[ftype, "particle_mass"] * m_frac * i_frac

            return _func

        ion_path = self._get_owls_ion_data_dir()
        fname = os.path.join(ion_path, ion + ".hdf5")
        itab = oit.IonTableOWLS(fname)
        return get_owls_ion_mass_field(ion, ftype, itab)

    # this function sets up the X_mass, X_density, X_fraction, and
    # X_number_density fields where X is the name of an OWLS element.
    # -------------------------------------------------------------
    def setup_fluid_fields(self):

        return

    # this function returns the owls_ion_data directory. if it doesn't
    # exist it will download the data from http://yt-project.org/data
    # -------------------------------------------------------------
    def _get_owls_ion_data_dir(self):

        txt = "Attempting to download ~ 30 Mb of owls ion data from %s to %s."
        data_file = "owls_ion_data.tar.gz"
        data_url = "http://yt-project.org/data"

        # get test_data_dir from yt config (ytcgf)
        # ----------------------------------------------
        tdir = ytcfg.get("yt", "test_data_dir")

        # set download destination to tdir or ./ if tdir isn't defined
        # ----------------------------------------------
        if tdir == "/does/not/exist":
            data_dir = "./"
        else:
            data_dir = tdir

        # check for owls_ion_data directory in data_dir
        # if not there download the tarball and untar it
        # ----------------------------------------------
        owls_ion_path = os.path.join(data_dir, "owls_ion_data")

        if not os.path.exists(owls_ion_path):
            mylog.info(txt, data_url, data_dir)
            fname = os.path.join(data_dir, data_file)
            download_file(os.path.join(data_url, data_file), fname)

            cmnd = f"cd {data_dir}; tar xf {data_file}"
            os.system(cmnd)

        if not os.path.exists(owls_ion_path):
            raise RuntimeError("Failed to download owls ion data.")

        return owls_ion_path
