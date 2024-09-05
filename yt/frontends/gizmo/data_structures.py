import os

import numpy as np

from yt.frontends.gadget.data_structures import GadgetHDF5Dataset
from yt.utilities.cosmology import Cosmology
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import GizmoFieldInfo


class GizmoDataset(GadgetHDF5Dataset):
    _load_requirements = ["h5py"]
    _field_info_class = GizmoFieldInfo

    @classmethod
    def _is_valid(cls, filename: str, *args, **kwargs) -> bool:
        if cls._missing_load_requirements():
            return False

        need_groups = ["Header"]
        veto_groups = ["Config", "Constants", "FOF", "Group", "Subhalo"]
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
            # From Apr 2021, 7f1f06f, public gizmo includes a header variable
            # GIZMO_version, which is set to the year of the most recent commit
            # We should prefer this to checking the metallicity, which might
            # not exist
            if "GIZMO_version" not in fh["/Header"].attrs:
                dmetal = "/PartType0/Metallicity"
                if dmetal not in fh or (
                    fh[dmetal].ndim > 1 and fh[dmetal].shape[1] < 11
                ):
                    valid = False
            fh.close()
        except Exception:
            valid = False
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

        self.current_redshift = hvals.get("Redshift", 0.0)
        if "Redshift" not in hvals:
            mylog.info("Redshift is not set in Header. Assuming z=0.")

        if "ComovingIntegrationOn" in hvals:
            # In 1d8479, Nov 2020, public GIZMO updated the names of the Omegas
            # to include an _, added baryons and radiation and added the
            # ComovingIntegrationOn field. ComovingIntegrationOn is always set,
            # but the Omega's are only included if ComovingIntegrationOn is true
            mylog.debug("Reading cosmological parameters using post-1d8479 format")
            self.cosmological_simulation = hvals["ComovingIntegrationOn"]
            if self.cosmological_simulation:
                self.omega_lambda = hvals["Omega_Lambda"]
                self.omega_matter = hvals["Omega_Matter"]
                self.omega_baryon = hvals["Omega_Baryon"]
                self.omega_radiation = hvals["Omega_Radiation"]
            self.hubble_constant = hvals["HubbleParam"]
        elif "OmegaLambda" in hvals:
            # Should still support GIZMO versions prior to 1d8479 too
            mylog.info(
                "ComovingIntegrationOn does not exist, falling back to OmegaLambda",
            )
            self.omega_lambda = hvals["OmegaLambda"]
            self.omega_matter = hvals["Omega0"]
            self.hubble_constant = hvals["HubbleParam"]
            self.cosmological_simulation = self.omega_lambda != 0.0
        else:
            # If these are not set it is definitely not a cosmological dataset.
            mylog.debug("No cosmological information found, assuming defaults")
            self.omega_lambda = 0.0
            self.omega_matter = 0.0  # Just in case somebody asks for it.
            self.cosmological_simulation = 0
            # Hubble is set below for Omega Lambda = 0.

        if not self.cosmological_simulation:
            mylog.info(
                "ComovingIntegrationOn != 1 or (not found "
                "and OmegaLambda is 0.0), so we are turning off Cosmology.",
            )
            self.hubble_constant = 1.0  # So that scaling comes out correct
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
            mylog.info(
                "Calculating time from %0.3e to be %0.3e seconds",
                hvals["Time"],
                self.current_time,
            )
        self.parameters = hvals

        prefix = os.path.join(self.directory, self.basename.split(".", 1)[0])

        if hvals["NumFiles"] > 1:
            self.filename_template = f"{prefix}.%(num)s{self._suffix}"
        else:
            self.filename_template = self.parameter_filename

        self.file_count = hvals["NumFiles"]
