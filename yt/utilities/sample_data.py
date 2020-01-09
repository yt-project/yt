"""
Title: sample_data.py
Purpose: Contains functions used for automatic downloading and loading
            of sample data that is not already present locally.
"""
import pooch

from yt.config import ytcfg
import pkg_resources
import json

## The format of the data registry json:
##
## {
##   'dataset_archive_name.tar.gz': {
##          'sha256': '...',
##          'url': '...',
##          'kwargs': {},
##          'item_name': 'supplied_to_load'
##   }
## }

class Fido:
    r"""
    Container for a pooch object used to fetch remote data that isn't
    already stored locally.
    """
    def __init__(self, filename = "sample_data_registry.json", cache_dir = None):
        self.filename = filename
        self._retrieve_registry()
        self._registry = json.load(pkg_resources.resource_stream("yt", self.filename))
        if cache_dir is None:
            if os.path.isdir(ytcfg.get("test_data_dir")):
                cache_dir = ytcfg.get("test_data_dir")
            else:
                cache_dir = pooch.os_cache("yt")
        self.fido = pooch.create(
            path=cache_dir,
            registry={_: self._registry[_]['sha256'] for _ in self._registry},
            urls={_: self._registry[_]['url'] for _ in self._registry},
            env="YT_DATA_DIR",
        )
        # Load the external registry file. It contains data file names,
        # hashes used for validation, and the url for the data file

    def __getitem__(self, item):
        if item in self._registry:
            return self._registry[item]
        if item + ".tar.gz" in self._registry:
            return self._registry[item + ".tar.gz"]
        raise KeyError(item)
