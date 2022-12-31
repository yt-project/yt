"""
    Chimera-specific IO functions



"""

import numpy as np
import unyt as un

from yt.frontends.chimera.data_structures import _find_files
from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py


class ChimeraIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "chimera"

    def __init__(self, ds):
        super().__init__(ds)
        self._handle = ds._handle
        self.filename = ds.filename

    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        pass

    def _read_particle_fields(self, chunks, ptf, selector):
        # This gets called after the arrays have been allocated.  It needs to
        # yield ((ptype, field), data) where data is the masked results of
        # reading ptype, field and applying the selector to the data read in.
        # Selector objects have a .select_points(x,y,z) that returns a mask, so
        # you need to do your masking here.
        pass

    def _read_fluid_selection(self, chunks, selector, fields, size):
        rv = {}
        nodal_fields = []
        for field in fields:
            finfo = self.ds.field_info[field]
            nodal_flag = finfo.nodal_flag
            if np.any(nodal_flag):
                num_nodes = 2 ** sum(nodal_flag)
                rv[field] = np.empty((size, num_nodes), dtype="=f8")
                nodal_fields.append(field)
            else:
                rv[field] = np.empty(size, dtype="=f8")
        ind = 0
        for field, mesh, data in self.io_iter(chunks, fields):
            if data is None:
                continue
            else:
                ind += mesh.select(selector, data.flatten(), rv[field], ind)  # caches
        return rv

    def io_iter(self, chunks, fields):
        for n, chunk in enumerate(chunks):
            file = _find_files(self.filename)
            with h5py.File(file[n], "r") as f:
                # Generates mask according to the "ongrid_mask" variable
                m = int(file[n][-5:-3]) - 1
                k = f["fluid"]["entropy"].shape[0]
                mask_0 = f["mesh"]["ongrid_mask"][k * m : k * (m + 1), :]

                if f["mesh"]["array_dimensions"][2] > 1:
                    nrd = f["mesh"]["array_dimensions"][0] - 2
                else:
                    nrd = f["mesh"]["array_dimensions"][0]

                mask = np.repeat(mask_0[:, :, np.newaxis], nrd, axis=2).transpose()
                specials = (
                    "abar",
                    "e_rms_1",
                    "e_rms_2",
                    "e_rms_3",
                    "e_rms_4",
                    "lumin_1",
                    "lumin_2",
                    "lumin_3",
                    "lumin_4",
                    "num_lumin_1",
                    "num_lumin_2",
                    "num_lumin_3",
                    "num_lumin_4",
                    "shock",
                    "nse_c",
                )
                for field in fields:  # Reads data by locating subheading
                    ftype, fname = field
                    a_name_2 = [i.decode("utf-8") for i in f["abundance"]["a_name"]]
                    a_name_dict = {name.strip(): name for name in a_name_2}
                    if fname not in specials:
                        if fname in f["fluid"]:
                            ds = f["fluid"][f"{fname}"]
                        elif fname in f["abundance"]:
                            ds = f["abundance"][f"{fname}"]
                        elif fname in a_name_dict:
                            ind_xn = a_name_2.index(a_name_dict[fname])
                            ds = f["abundance"]["xn_c"][:, :, :, ind_xn]
                        else:
                            mylog.warning("Invalid field name %s", fname)
                        dat_1 = ds[:, :, :].transpose()

                    elif fname == "nse_c":
                        if np.shape(f["abundance"]["nse_c"]) != np.shape(
                            f["fluid"]["rho_c"]
                        ):
                            ds = f["abundance"]["nse_c"][:, :, 1:]
                        else:
                            ds = f["abundance"]["nse_c"]
                        dat_1 = ds[:, :, :].transpose()
                    elif fname == "abar":
                        xn_c = np.array(f["abundance"]["xn_c"])
                        a_nuc_rep_c = np.array(f["abundance"]["a_nuc_rep_c"])
                        a_nuc = np.array(f["abundance"]["a_nuc"])
                        a_nuc_tile = np.tile(
                            a_nuc, (xn_c.shape[0], xn_c.shape[1], xn_c.shape[2], 1)
                        )
                        yn_c = np.empty(xn_c.shape)
                        yn_c[:, :, :, :-1] = xn_c[:, :, :, :-1] / a_nuc_tile[:, :, :, :]
                        yn_c[:, :, :, -1] = xn_c[:, :, :, -1] / a_nuc_rep_c[:, :, :]
                        ytot = np.sum(yn_c, axis=3)
                        atot = np.sum(xn_c, axis=3)
                        abar = np.divide(atot, ytot)
                        dat_1 = abar[:, :, :].transpose()
                    elif fname in ("e_rms_1", "e_rms_2", "e_rms_3", "e_rms_4"):
                        dims = f["mesh"]["array_dimensions"]
                        n_groups = f["radiation"]["raddim"][0]
                        n_species = f["radiation"]["raddim"][1]
                        n_hyperslabs = f["mesh"]["nz_hyperslabs"][()]
                        energy_edge = f["radiation"]["unubi"][()]
                        energy_center = f["radiation"]["unui"][()]
                        d_energy = []
                        for i in range(0, n_groups):
                            d_energy.append(energy_edge[i + 1] - energy_edge[i])
                        d_energy = np.array(d_energy)
                        e3de = energy_center**3 * d_energy
                        e5de = energy_center**5 * d_energy

                        psi0_c = f["radiation"]["psi0_c"][:]
                        row = np.empty(
                            (n_species, int(dims[2] / n_hyperslabs), dims[1], dims[0])
                        )
                        for n in range(0, n_species):
                            numerator = np.sum(psi0_c[:, :, :, n] * e5de, axis=3)
                            denominator = np.sum(psi0_c[:, :, :, n] * e3de, axis=3)
                            row[n][:][:][:] = np.sqrt(
                                numerator / (denominator + 1e-100)
                            )
                            species = int(fname[-1]) - 1
                            dat_1 = row[species, :, :, :].transpose()
                    elif fname in (
                        "lumin_1",
                        "lumin_2",
                        "lumin_3",
                        "lumin_4",
                        "num_lumin_1",
                        "num_lumin_2",
                        "num_lumin_3",
                        "num_lumin_4",
                    ):
                        dims = f["mesh"]["array_dimensions"]
                        n_groups = f["radiation"]["raddim"][0]
                        n_hyperslabs = f["mesh"]["nz_hyperslabs"][()]
                        ergmev = float((1 * un.MeV) / (1 * un.erg))
                        cvel = float(un.c.to("cm/s"))
                        h = float(un.h.to("MeV * s"))
                        ecoef = 4.0 * np.pi * ergmev / (h * cvel) ** 3
                        radius = f["mesh"]["x_ef"][()]
                        agr_e = f["fluid"]["agr_e"][()]
                        cell_area_GRcorrected = 4 * np.pi * radius**2 / agr_e**4
                        psi1_e = f["radiation"]["psi1_e"]
                        energy_edge = f["radiation"]["unubi"][()]
                        energy_center = f["radiation"]["unui"][()]
                        d_energy = []
                        for i in range(0, n_groups):
                            d_energy.append(energy_edge[i + 1] - energy_edge[i])
                        d_energy = np.array(d_energy)
                        species = int(fname[-1]) - 1
                        if fname in ("lumin_1", "lumin_2", "lumin_3", "lumin_4"):
                            eNde = energy_center**3 * d_energy
                        else:
                            eNde = energy_center**2 * d_energy

                        lumin = (
                            np.sum(psi1_e[:, :, :, species] * eNde, axis=3)
                            * np.tile(
                                cell_area_GRcorrected[1 : dims[0] + 1],
                                (int(dims[2] / n_hyperslabs), dims[1], 1),
                            )
                            * (cvel * ecoef * 1e-51)
                        )
                        dat_1 = lumin[:, :, :].transpose()

                    if f["mesh"]["array_dimensions"][2] > 1:
                        data = dat_1[:-2, :, :]  # Clips off ghost zones for 3D
                    else:
                        data = dat_1[:, :, :]
                        data = np.ma.masked_where(mask == 0.0, data)  # Masks
                        data = np.ma.filled(
                            data, fill_value=0.0
                        )  # Replaces masked value with 0

                    yield field, chunk.objs[0], data

        pass

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk without doing any selection,
        # and is only used for caching data that might be used by multiple
        # different selectors later. For instance, this can speed up ghost zone
        # computation.
        pass
