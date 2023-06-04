import os

import numpy as np

from yt.frontends.gadget.data_structures import GadgetHDF5Dataset
from yt.funcs import only_on_root
from yt.utilities.cosmology import Cosmology
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import GizmoFieldInfo


class GizmoDataset(GadgetHDF5Dataset):
    _field_info_class = GizmoFieldInfo

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        need_groups = ["Header"]
        veto_groups = ["FOF", "Group", "Subhalo"]
        valid = True
        valid_fname = filename
        # If passed arg is a directory, look for the .0 file in that dir
        if os.path.isdir(filename):
            valid_files = []
            for f in os.listdir(filename):
                fname = os.path.join(filename, f)
                fext = os.path.splitext(fname)[-1]
                if (
                    (".0" in f)
                    and (fext not in {".ewah", ".kdtree"})
                    and os.path.isfile(fname)
                ):
                    valid_files.append(fname)
            if len(valid_files) == 0:
                valid = False
            elif len(valid_files) > 1:
                valid = False
            else:
                valid_fname = valid_files[0]
        try:
            fh = h5py.File(valid_fname, mode="r")
            valid = all(ng in fh["/"] for ng in need_groups) and not any(
                vg in fh["/"] for vg in veto_groups
            )
            dmetal = "/PartType0/Metallicity"
            if dmetal not in fh or (fh[dmetal].ndim > 1 and fh[dmetal].shape[1] < 11):
                valid = False
            fh.close()
        except Exception:
            valid = False
            pass
        return valid

    def _set_code_unit_attributes(self):
        super()._set_code_unit_attributes()

    def _parse_parameter_file(self):
        hvals = self._get_hvals()

        self.dimensionality = 3
        self.refine_by = 2
        self.parameters["HydroMethod"] = "sph"
        # Set standard values

        # We may have an overridden bounding box.
        if self.domain_left_edge is None and hvals["BoxSize"] != 0:
            self.domain_left_edge = np.zeros(3, "float64")
            self.domain_right_edge = np.ones(3, "float64") * hvals["BoxSize"]

        self.domain_dimensions = np.ones(3, "int32")
        self._periodicity = (True, True, True)

        self.cosmological_simulation = 1

        try:
            self.current_redshift = hvals["Redshift"]
        except KeyError:
            # Probably not a cosmological dataset, we should just set
            # z = 0 and let the user know
            self.current_redshift = 0.0
            only_on_root(mylog.info, "Redshift is not set in Header. Assuming z=0.")

        try:
            # The current version of GIZMO has updated the names of the Omegas
            # to include an _, added baryons and radiation and added the
            # ComovingIntegrationOn field. ComovingIntegrationOn is always set,
            # but the Omega's are only included if ComovingIntegrationOn is true
            if "ComovingIntegrationOn" in hvals:
                self.cosmological_simulation = hvals["ComovingIntegrationOn"]
                self.omega_lambda = hvals["Omega_Lambda"]
                self.omega_matter = hvals["Omega_Matter"]
                self.omega_baryon = hvals["Omega_Baryon"]
                self.omega_radiation = hvals["Omega_Radiation"]
                self.hubble_constant = hvals["HubbleParam"]
            else:
                # Should still support old GIZMO versions too
                only_on_root(
                    mylog.info,
                    "ComovingIntegrationOn does not exist, falling back to OmegaLambda",
                )
                self.omega_lambda = hvals["OmegaLambda"]
                self.omega_matter = hvals["Omega0"]
                self.hubble_constant = hvals["HubbleParam"]
                self.cosmological_simulation = self.omega_lambda == 0.0
        except KeyError:
            # If these are not set it is definitely not a cosmological dataset.
            self.omega_lambda = 0.0
            self.omega_matter = 1.0  # Just in case somebody asks for it.
            self.cosmological_simulation = 0
            # Hubble is set below for Omega Lambda = 0.

        # Since new versions of GIZMO output the ComovingIntegrationOn flag
        # we will rely on its presence to determine if this is a cosmological
        # dataset.
        if not self.cosmological_simulation:
            only_on_root(
                mylog.info,
                "ComovingIntegrationOn != 1 or (not found "
                "and OmegaLambda is 0.0), so we are turning off Cosmology.",
            )
            self.hubble_constant = 1.0  # So that scaling comes out correct
            self.cosmological_simulation = 0
            self.current_redshift = 0.0
            # This may not be correct.
            self.current_time = hvals["Time"]
        else:
            # Now we calculate our time based on the cosmology, because in
            # ComovingIntegration hvals["Time"] will in fact be the expansion
            # factor, not the actual integration time, so we re-calculate
            # global time from our Cosmology.
            cosmo = Cosmology(
                hubble_constant=self.hubble_constant,
                omega_matter=self.omega_matter,
                omega_lambda=self.omega_lambda,
            )
            self.current_time = cosmo.lookback_time(self.current_redshift, 1e6)
            only_on_root(
                mylog.info,
                "Calculating time from %0.3e to be %0.3e seconds",
                hvals["Time"],
                self.current_time,
            )
        self.parameters = hvals

        prefix = os.path.abspath(
            os.path.join(
                os.path.dirname(self.parameter_filename),
                os.path.basename(self.parameter_filename).split(".", 1)[0],
            )
        )

        if hvals["NumFiles"] > 1:
            for t in (
                f"{prefix}.%(num)s{self._suffix}",
                f"{prefix}.gad.%(num)s{self._suffix}",
            ):
                if os.path.isfile(t % {"num": 0}):
                    self.filename_template = t
                    break
            else:
                raise RuntimeError("Could not determine correct data file template.")
        else:
            self.filename_template = self.parameter_filename

        self.file_count = hvals["NumFiles"]
