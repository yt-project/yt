"""
This is a collection of helper functions to yt.load_sample
"""
import json
import re
from functools import lru_cache
from itertools import chain
from pathlib import Path
from typing import Optional, Union

import pkg_resources

from yt.funcs import mylog
from yt.utilities.on_demand_imports import (
    _pandas as pd,
    _pooch as pooch,
    _requests as requests,
)

num_exp = re.compile(r"\d*(\.\d*)?")
byte_unit_exp = re.compile(r"[KMG]?B")


def _parse_byte_size(s: str):
    """
    Convert a string size specification to integer byte size.

    This function should be insensitive to case and whitespace.
    Examples
    --------
    >>> _parse_byte_size("11B")
    11
    >>> _parse_byte_size("10 kB")
    10000
    >>> _parse_byte_size("67.4GB")
    67400000000
    >>> _parse_byte_size("0.1 MB")
    100000
    >>> _parse_byte_size(".1 MB")
    100000
    >>> _parse_byte_size(None)
    <NA>
    >>> from numpy import nan
    >>> _parse_byte_size(nan)
    <NA>
    """
    try:
        s = s.upper()
    except AttributeError:
        # input is not a string (likely a np.nan)
        return pd.NA
    val = float(re.search(num_exp, s).group())
    unit = re.search(byte_unit_exp, s).group()
    exponent = {"B": 0, "KB": 3, "MB": 6, "GB": 9}[unit]
    return int(val * 10 ** exponent)


@lru_cache(maxsize=128)
def get_data_registry_table():
    """
    Load the sample data registry as a pandas.Dataframe instance.

    This function is considered experimental and is exposed for exploratory purposed.
    The output format is subject to change.

    The output of this function is cached so it will only generate one request per session.
    """

    # it would be nicer to have an actual api on the yt website server,
    # but this will do for now
    api_url = "https://raw.githubusercontent.com/yt-project/website/master/data/datafiles.json"  # noqa

    response = requests.get(api_url)

    if not response.ok:
        raise RuntimeError(
            "Could not retrieve registry data. Please check your network setup."
        )

    website_json = response.json()
    # this dict follows this schema: {frontend_name: {flat dataframe-like}}

    columns = ["code", "filename", "size", "url", "description"]
    website_table = pd.concat(pd.DataFrame(d) for d in website_json.values())[columns]

    # add a int-type byte size column
    # note that we cast to pandas specific type "Int64" because we expect missing values
    # see https://pandas.pydata.org/pandas-docs/stable/user_guide/missing_data.html#integer-dtypes-and-missing-data
    website_table["byte_size"] = (
        website_table["size"].apply(_parse_byte_size).astype("Int64")
    )

    # normalize urls to match the local json
    website_table["url"] = website_table["url"].apply(
        lambda u: u.replace("http:", "https:")
    )

    # load local data
    pooch_json = json.load(
        pkg_resources.resource_stream("yt", "sample_data_registry.json")
    )
    pooch_table = pd.DataFrame(pooch_json.values())

    # merge tables
    unified_table = website_table.merge(pooch_table, on="url", how="outer")

    # PR 3089
    # ideally we should be able to do this, but it's not possible
    # at the time of writing because fhe "filename" is incomplete
    # see the companion comment in load_sample
    # unified_table.set_index("filename", inplace=True)
    # unified_table.index.rename("id", inplace=True)
    return unified_table


def _get_test_data_dir_path():
    from yt.config import ytcfg

    return Path(ytcfg.get("yt", "test_data_dir"))


def lookup_on_disk_data(fn):
    """
    Look for data file/dir on disk.

    Returns
    -------

    pathlib.Path to a file/dir matching fn if any

    Raises
    ------
    FileNotFoundError
    """

    path = Path(fn).expanduser()

    if path.exists():
        return path

    alt_path = _get_test_data_dir_path() / fn
    err_msg = f"No such file or directory: '{fn}'."
    if not alt_path.parent.is_dir():
        raise FileNotFoundError(err_msg)

    if not alt_path.exists():
        err_msg += f"\n(Also tried '{alt_path}')."
        raise FileNotFoundError(err_msg)

    return alt_path


@lru_cache(maxsize=128)
def _get_pooch_instance():
    data_registry = get_data_registry_table()
    cache_storage = _get_test_data_dir_path() / "yt_download_cache"

    sample_data_registry = json.load(
        pkg_resources.resource_stream("yt", "sample_data_registry.json")
    )
    registry = {k: v["hash"] for k, v in sample_data_registry.items()}
    return pooch.create(
        path=cache_storage, base_url="https://yt-project.org/data/", registry=registry
    )


def _download_sample_data_file(
    filename, progressbar=True, timeout: Optional[int] = None
):
    """
    Download a file by url.

    Returns
    -------
    storage_filename : location of the downloaded file

    """
    downloader = pooch.HTTPDownloader(progressbar=progressbar, timeout=timeout)

    poochie = _get_pooch_instance()
    poochie.fetch(filename, downloader=downloader)
    return Path.joinpath(poochie.path, filename)
