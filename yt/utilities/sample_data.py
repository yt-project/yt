"""
Title: sample_data.py
Purpose: Contains functions used for automatic downloading and loading
            of sample data that is not already present locally.
"""
import pooch

from yt.config import ytcfg


class Fido:
    r"""
    Container for a pooch object used to fetch remote data that isn't
    already stored locally.
    """
    def __init__(self):
        self.fido = pooch.create(
            path=pooch.os_cache("yt"),
            registry=None,
            env="YT_DATA_DIR",
        )
        # Load the external registry file. It contains data file names,
        # hashes used for validation, and the url for the data file
        self.fido.load_registry("sample_data_registry.txt")
