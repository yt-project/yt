from yt.frontends.gadget.api import GadgetHDF5Dataset
from yt.funcs import mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import \
    ArepoFieldInfo

import numpy as np


class ArepoHDF5Dataset(GadgetHDF5Dataset):
    _field_info_class = ArepoFieldInfo

    def __init__(self, filename, dataset_type="arepo_hdf5",
                 unit_base=None,
                 smoothing_factor=2.0,
                 index_order=None,
                 index_filename=None,
                 kernel_name=None,
                 bounding_box=None,
                 units_override=None,
                 unit_system="cgs"):
        super(ArepoHDF5Dataset, self).__init__(
            filename, dataset_type=dataset_type, unit_base=unit_base,
            index_order=index_order, index_filename=index_filename,
            kernel_name=kernel_name, bounding_box=bounding_box,
            units_override=units_override, unit_system=unit_system)
        # The "smoothing_factor" is a user-configurable parameter which
        # is multiplied by the radius of the sphere with a volume equal
        # to that of the Voronoi cell to create smoothing lengths.
        self.smoothing_factor = smoothing_factor
        self.gamma = 5./3.

    @classmethod
    def _is_valid(self, *args, **kwargs):
        need_groups = ['Header', 'Config']
        veto_groups = ['FOF', 'Group', 'Subhalo']
        valid = True
        try:
            fh = h5py.File(args[0], mode='r')
            valid = all(ng in fh["/"] for ng in need_groups) and \
                    not any(vg in fh["/"] for vg in veto_groups) and \
                    ("VORONOI" in fh["/Config"].attrs.keys() or
                     "AMR" in fh["/Config"].attrs.keys())
            fh.close()
        except:
            valid = False
            pass
        return valid

    def _get_uvals(self):
        handle = h5py.File(self.parameter_filename, mode="r")
        uvals = {}
        missing = False
        for unit in ["UnitLength_in_cm", "UnitMass_in_g", 
                     "UnitVelocity_in_cm_per_s"]:
            if unit in handle["/Header"].attrs:
                uvals[unit] = handle["/Header"].attrs[unit]
            else:
                mylog.warning("Arepo header is missing %s!" % unit)
                missing = True
        handle.close()
        if missing:
            uvals = None
        return uvals

    def _set_code_unit_attributes(self):
        self._unit_base = self._get_uvals()
        super(ArepoHDF5Dataset, self)._set_code_unit_attributes()
        munit = np.sqrt(self.mass_unit /
                        (self.time_unit**2 * self.length_unit)).to("gauss")
        if self.cosmological_simulation:
            self.magnetic_unit = self.quan(munit.value, "%s/a**2" % munit.units)
        else:
            self.magnetic_unit = munit
