"""
sample data manager for yt

This utility will check to see if sample data exists on disc.
If not, it will download it.

"""

import logging
from yt.utilities.sample_data import Fido
from yt.funcs import mylog
from yt.convenience import load
from pooch import Untar
import os.path.commonprefix as commonprefix

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

    base_path = Fido.fido.env
    fileext, name, extension = _validate_sampledata_name

    if extension == "h5":
        fname = fetch_noncompressed_file(fileext)
    else:
        # we are going to assume most files that exist on the hub are
        # compressed in .tar folders. Some may not.
        fname = fetch_compressed_file(fileext)

    folder_path = commonprefix(fname)
    mylog.info("Files located at %s", folder_path)

    # Location of the file to load automatically, registered in the Fido class
    file_lookup = Fido.fido.filenames[fileext][loadname]

    if file is None:
        loaded_file = "%s/%s/%s" %(base_path,name,file_lookup)
    else:
        loaded_file = "%s/%s/%s" %(base_path,name,file)

    load(loaded_file)

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
    assert(name, str)

    # now get the extension if it exists
    (base,ext) = os.path.splitext(name)
    if ext == '':
        fileext = "%s.tar.gz" %name
        basename = name
        extension = "tar"
    elif ext == ".gz":
        fileext = name
        basename = os.path.splitext(ext)[-1].lower()
        extension = "tar"
    elif ext == ".h5" or elif ext == ".hdf5":
        fileext = name
        basename = base
        extension = "h5"
    else:
        mylog.info("extension of %s for dataset %s is unexpected. pooch may not
                work", %(ext, name))
        extension = ext
        fileext = name
        basename = base
    return fileext, basename, extension


def fetch_compressed_file(name):
    """
    Load a large compressed file from the data registry
    """
    fname = Fido.fido.fetch(name, processor=Untar())
    return fname

def fetch_noncompressed_file(name):
    """
    Load an uncompressed file from the data registry
    """
    fname = Fido.fido.fetch(name)
    return fname

