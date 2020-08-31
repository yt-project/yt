"""
Title: sample_data.py
Purpose: Contains functions used for automatic downloading and loading
            of sample data that is not already present locally.
"""
import json
import os

import pkg_resources

from yt.config import ytcfg
from yt.funcs import mylog
from yt.utilities.on_demand_imports import _pooch as pooch

## The format of the data registry json:
##
## {
##   'dataset_archive_name.tar.gz': {
##          'hash': '...',
##          'url': '...',
##          'load_kwargs': {},
##          'load_name': 'supplied_to_load'
##   }
## }

_extensions_to_strip = (".tgz", ".tar.gz", ".gz")


class PoochHandle:
    r"""
    Container for a pooch object used to fetch remote data that isn't
    already stored locally.
    """

    def __init__(self, filename="sample_data_registry.json", cache_dir=None):
        self.filename = filename
        self._registry = json.load(pkg_resources.resource_stream("yt", self.filename))
        if cache_dir is None:
            if os.path.isdir(ytcfg.get("yt", "test_data_dir")):
                cache_dir = ytcfg.get("yt", "test_data_dir")
            else:
                cache_dir = pooch.pooch.os_cache("yt")
        self.pooch_obj = pooch.pooch.create(
            path=cache_dir,
            registry={_: self._registry[_]["hash"] for _ in self._registry},
            urls={_: self._registry[_]["url"] for _ in self._registry},
            env="YT_DATA_DIR",
            base_url="https://yt-project.org/data/",
        )
        # Load the external registry file. It contains data file names,
        # hashes used for validation, and the url for the data file

    def __getitem__(self, item):
        if item in self._registry:
            return self._registry[item]
        for ext in _extensions_to_strip:
            if item + ext in self._registry:
                return self._registry[item + ext]
        raise KeyError(item)

    @staticmethod
    def _validate_sample_fname(name):
        """
        format name of sample data passed to function, accepts a named string
        argument and parses it to determine the sample data name, what type of
        extension it has, or other relevant information.

        Returns
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
            mylog.error(
                "The argument %s passed to load_sample() is not a string.", name
            )

        # now get the extension if it exists
        base, ext = os.path.splitext(name)
        if ext == "":
            # Right now we are assuming that any name passed without an explicit
            # extension is packed in a tarball. This logic can be modified later to
            # be more flexible.
            fileext = f"{name}.tar.gz"
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
                ext,
                name,
            )
            extension = ext
            fileext = name
            basename = base
        return fileext, basename, extension
