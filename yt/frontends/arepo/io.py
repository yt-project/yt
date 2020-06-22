from yt.frontends.gadget.api import IOHandlerGadgetHDF5
import numpy as np
from yt.utilities.on_demand_imports import _h5py as h5py

class IOHandlerArepoHDF5(IOHandlerGadgetHDF5):
    _dataset_type = "arepo_hdf5"

    def _get_smoothing_length(self, data_file, position_dtype, position_shape):
        ptype = self.ds._sph_ptypes[0]
        ind = int(ptype[-1])
        si, ei = data_file.start, data_file.end
        with h5py.File(data_file.filename, "r") as f:
            pcount = f["/Header"].attrs["NumPart_ThisFile"][ind].astype("int")
            pcount = np.clip(pcount - si, 0, ei - si)
            # Arepo cells do not have "smoothing lengths" by definition, so
            # we compute one here by finding the radius of the sphere
            # corresponding to the volume of the Voroni cell and multiplying
            # by a user-configurable smoothing factor.
            hsml = f[ptype]["Masses"][si:ei,...]/f[ptype]["Density"][si:ei,...]
            hsml *= 3.0/(4.0*np.pi)
            hsml **= (1./3.)
            hsml *= self.ds.smoothing_factor
            dt = hsml.dtype.newbyteorder("N") # Native
            if position_dtype is not None and dt < position_dtype:
                dt = position_dtype
            return hsml.astype(dt)

    def _identify_fields(self, data_file):
        fields, _units = super(IOHandlerArepoHDF5, 
                               self)._identify_fields(data_file)
        fields.append(("PartType0", "smoothing_length"))
        return fields, _units

