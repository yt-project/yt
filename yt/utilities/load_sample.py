"""
sample data manager for yt

This utility will check to see if sample data exists on disc.
If not, it will download it.

"""

import yt.utilities.sample_data as sd
from yt.funcs import mylog
from yt.convenience import load
from yt.utilities.on_demand_imports import _pooch as pch

def load_sample(name, file=None):
    """
    Load sample data with yt. Simple wrapper around yt.load to include fetching
    data with pooch.

    Parameters
    ----------
    name : str
        The name of the sample data to load. This is generally the name of the
        folder of the dataset. For IsolatedGalaxy, the name would be
        `IsolaatedGalaxy`

    file : str, optionaal
        optional argument -- the name of the file to load that is located
        within sample dataset of `name`. For the file
        IsolatedGalaxy/galaxy0030/galaxy0030, this variable would be
        `galaxy0030/galaxy0030`

    """
    fido = sd.Fido()

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
    folder_path = commonprefix(fname)
    mylog.info("Files located at %s", folder_path)

    # Location of the file to load automatically, registered in the Fido class
    info = fido[fileext]
    file_lookup = info['load_name']
    optional_args = info['load_kwargs']

    if file is None:
        # right now work on loading only untarred files. build out h5 later
        loaded_file = "%s/%s.untar/%s/%s" %(base_path,fileext,name,file_lookup)
    else:
        loaded_file = "%s/%s.untar/%s/%s" %(base_path,fileext,name,file)

    load(loaded_file, **optional_args)

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

    assert isinstance(name, str)

    # now get the extension if it exists
    (base,ext) = splitext(name)
    if ext == '':
        fileext = "%s.tar.gz" %name
        basename = name
        extension = "tar"
    elif ext == ".gz":
        fileext = name
        basename = splitext(ext)[-1].lower()
        extension = "tar"
    elif ext == ".h5" or ext == ".hdf5":
        fileext = name
        basename = base
        extension = "h5"
    else:
        mylog.info(
            "extension of %s for dataset %s is unexpected. pooch may not work",
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

