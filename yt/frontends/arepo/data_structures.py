import numpy as np

from yt.frontends.gadget.api import GadgetHDF5Dataset
from yt.funcs import mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import ArepoFieldInfo


class ArepoHDF5Dataset(GadgetHDF5Dataset):
    _field_info_class = ArepoFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="arepo_hdf5",
        unit_base=None,
        smoothing_factor=2.0,
        index_order=None,
        index_filename=None,
        kernel_name=None,
        bounding_box=None,
        units_override=None,
        unit_system="cgs",
        default_species_fields=None,
    ):
        super().__init__(
            filename,
            dataset_type=dataset_type,
            unit_base=unit_base,
            index_order=index_order,
            index_filename=index_filename,
            kernel_name=kernel_name,
            bounding_box=bounding_box,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
        )
        # The "smoothing_factor" is a user-configurable parameter which
        # is multiplied by the radius of the sphere with a volume equal
        # to that of the Voronoi cell to create smoothing lengths.
        self.smoothing_factor = smoothing_factor
        self.gamma = 5.0 / 3.0
        self.gamma_cr = self.parameters.get("GammaCR", 4.0 / 3.0)

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        need_groups = ["Header", "Config"]
        veto_groups = ["FOF", "Group", "Subhalo"]
        valid = True
        try:
            fh = h5py.File(filename, mode="r")
            valid = (
                all(ng in fh["/"] for ng in need_groups)
                and not any(vg in fh["/"] for vg in veto_groups)
                and (
                    "VORONOI" in fh["/Config"].attrs.keys()
                    or "AMR" in fh["/Config"].attrs.keys()
                )
                # Datasets with GFM_ fields present are AREPO
                or any(field.startswith("GFM_") for field in fh["/PartType0"])
            )
            fh.close()
        except Exception:
            valid = False
        return valid

    def _get_uvals(self):
        handle = h5py.File(self.parameter_filename, mode="r")
        uvals = {}
        missing = [False] * 3
        for i, unit in enumerate(
            ["UnitLength_in_cm", "UnitMass_in_g", "UnitVelocity_in_cm_per_s"]
        ):
            if unit in handle["/Header"].attrs:
                uvals[unit] = handle["/Header"].attrs[unit]
                if unit == "UnitLength_in_cm":
                    # We assume this is comoving, because in the absence of comoving
                    # integration the redshift will be zero.
                    uvals["cmcm"] = 1.0 / uvals[unit]
            else:
                missing[i] = True
        handle.close()
        if all(missing):
            uvals = None
        return uvals

    def _set_code_unit_attributes(self):
        arepo_unit_base = self._get_uvals()
        # This rather convoluted logic is required to ensure that
        # units which are present in the Arepo dataset will be used
        # no matter what but that the user gets warned
        if arepo_unit_base is not None:
            if self._unit_base is None:
                self._unit_base = arepo_unit_base
            else:
                for unit in arepo_unit_base:
                    if unit == "cmcm":
                        continue
                    short_unit = unit.split("_")[0][4:].lower()
                    if short_unit in self._unit_base:
                        which_unit = short_unit
                        self._unit_base.pop(short_unit, None)
                    elif unit in self._unit_base:
                        which_unit = unit
                    else:
                        which_unit = None
                    if which_unit is not None:
                        msg = f"Overwriting '{which_unit}' in unit_base with what we found in the dataset."
                        mylog.warning(msg)
                    self._unit_base[unit] = arepo_unit_base[unit]
                if "cmcm" in arepo_unit_base:
                    self._unit_base["cmcm"] = arepo_unit_base["cmcm"]
        super()._set_code_unit_attributes()
        munit = np.sqrt(self.mass_unit / (self.time_unit**2 * self.length_unit)).to(
            "gauss"
        )
        if self.cosmological_simulation:
            self.magnetic_unit = self.quan(munit.value, f"{munit.units}/a**2")
        else:
            self.magnetic_unit = munit
