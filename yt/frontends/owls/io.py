from yt.frontends.gadget.io import IOHandlerGadgetHDF5


class IOHandlerOWLS(IOHandlerGadgetHDF5):
    _dataset_type = "OWLS"
