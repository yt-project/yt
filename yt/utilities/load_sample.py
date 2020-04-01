"""
sample data manager for yt

This utility will check to see if sample data exists on disc.
If not, it will download it.

"""
import os
import yt.utilities.sample_data as sd
from yt.funcs import mylog
from yt.convenience import load
from yt.utilities.on_demand_imports import _pooch as pch

def load_sample(name=None, specific_file=None):
    """
    Load sample data with yt. Simple wrapper around yt.load to include fetching
    data with pooch.

    Parameters
    ----------
    name : str or None
        The name of the sample data to load. This is generally the name of the
        folder of the dataset. For IsolatedGalaxy, the name would be
        `IsolatedGalaxy`.  If `None` is supplied, the return value
        will be a list of all known datasets (by name).

    specific_file : str, optional
        optional argument -- the name of the file to load that is located
        within sample dataset of `name`. For the dataset `enzo_cosmology_plus`,
        which has a number of timesteps available, one may wish to choose
        DD0003. The file specifically would be
        `enzo_cosmology_plus/DD0003/DD0003`, and the argument passed to this
        variable would be `DD0003/DD0003`

    """
    fido = sd.Fido()
    if name is None:
        keys = []
        for key in fido._registry:
            for ext in sd._extensions_to_strip:
                if key.endswith(ext): key = key[:-len(ext)]
            keys.append(key)
        return keys

    base_path = fido.fido.path
    fileext, name, extension = _validate_sampledata_name(name)

    if extension == "h5":
        fname = fetch_noncompressed_file(fileext, fido)
    else:
        # we are going to assume most files that exist on the hub are
        # compressed in .tar folders. Some may not.
        fname = fetch_compressed_file(fileext, fido)

    # The `folder_path` variable is used here to notify the user where the
    # files have been unpacked to. However, we can't assume this is reliable
    # because in some cases the common path will overlap with the `load_name`
    # variable of the file.
    folder_path = os.path.commonprefix(fname)
    mylog.info("Files located at %s", folder_path)

    # Location of the file to load automatically, registered in the Fido class
    info = fido[fileext]
    file_lookup = info['load_name']
    optional_args = info['load_kwargs']

    if specific_file is None:
        # right now work on loading only untarred files. build out h5 later
        mylog.info("Default to loading %s for %s dataset", file_lookup, name)
        loaded_file = os.path.join(base_path, "%s.untar" %fileext,
                                   name, file_lookup)
    else:
        mylog.info("Loading %s for %s dataset", specific_file, name)
        loaded_file = os.path.join(base_path, "%s.untar" %fileext, name,
                                   specific_file)

    return load(loaded_file, **optional_args)

def _validate_sampledata_name(name):
    """
    format name of sample data passed to function, accepts a named string
    argument and parses it to determine the sample data name, what type of
    extension it has, or other relevant information.

    returns
    -------
    fileext : str
        The name of the sample data, with the file extension
        example: "IsolatedGalaxy.tar.gz"
    basename : str
        The name of the sample data, without the file extension
        example: "IsolatedGalaxy"
    extension : str
        name of extension of remote sample data
        example: "h5" or "tar"
    """

    if not isinstance(name, str):
        mylog.error("The argument {} passed to ".format(name) + \
                    "load_sample() is not a string.")

    # now get the extension if it exists
    base, ext = os.path.splitext(name)
    if ext == '':
        # Right now we are assuming that any name passed without an explicit
        # extension is packed in a tarball. This logic can be modified later to
        # be more flexible.
        fileext = "%s.tar.gz" %name
        basename = name
        extension = "tar"
    elif ext == ".gz":
        fileext = name
        basename = os.path.splitext(base)[0]
        extension = "tar"
    elif ext in [".h5", ".hdf5"]:
        fileext = name
        basename = base
        extension = "h5"
    else:
        mylog.info(
            """extension of %s for dataset %s is unexpected. the `load_data`
            function  may not work as expected""",
            ext, name )
        extension = ext
        fileext = name
        basename = base
    return fileext, basename, extension


def fetch_compressed_file(name, fido):
    """
    Load a large compressed file from the data registry
    """
    fname = fido.fido.fetch(name, processor=pch.pooch.Untar())
    return fname

def fetch_noncompressed_file(name, fido):
    """
    Load an uncompressed file from the data registry
    """
    fname = fido.fido.fetch(name)
    return fname

