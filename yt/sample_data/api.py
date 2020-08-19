import json
import re
from itertools import chain
from pathlib import Path

import pandas as pd
import pkg_resources
import pooch
import requests

from yt.funcs import mylog

num_exp = re.compile(r"[0-9]+")
byte_unit_exp = re.compile(r"[KMG]?B")


def _parse_byte_size(s):
    """
    Convert a string size specification to integer byte size.

    This function should be insensitive to case and whitespace.
    Examples
    --------
    >>> parse_byte_size("11B")
    11
    >>> parse_byte_size("10 kB")
    10000
    >>> parse_byte_size("67.4GB")
    67400000000
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


def _get_archive_filename(url):
    archive_filename = url.split("/")[-1]
    return archive_filename


sample_data_registry = json.load(
    pkg_resources.resource_stream("yt", "sample_data_registry.json")
)


def get_data_registry_table():
    """
    Load the sample data registry as a pandas.Dataframe instance.

    Returns
    -------
    return_type="dataframe" (default)
    todo: describe table format

    """

    # it would be nicer to have an actual api on the yt website server,
    # but this will do for now
    api_url = "https://raw.githubusercontent.com/yt-project/website/master/data/datafiles.json"  # noqa

    response = requests.get(api_url)

    if not response.ok:
        raise RuntimeError(
            "Could not retrieve registry data. Please check your network setup."
        )

    # todo: we should cache this the first time it's loaded
    website_json = json.loads(response.text)
    # this dict follows this schema: {frontend_name: {flat dataframe-like}}

    columns = ["code", "filename", "size", "url"]
    website_table = pd.concat(pd.DataFrame(d) for d in website_json.values())[columns]

    # add a int-type byte size column
    website_table["byte size"] = website_table["size"].apply(_parse_byte_size)

    # normalize urls to match the local json
    website_table["url"] = website_table["url"].apply(
        lambda u: u.replace("http", "https")
    )

    # load local data
    pooch_json = json.load(
        pkg_resources.resource_stream("yt", "sample_data_registry.json")
    )
    pooch_table = pd.DataFrame(pooch_json.values())

    # merge tables
    unified_table = website_table.merge(pooch_table, on="url", how="outer")

    # unified_table.set_index("filename", inplace=True)
    # unified_table.index.rename("id", inplace=True)
    return unified_table


def _get_test_data_dir_path():
    from yt.config import ytcfg

    test_data_dir = Path(ytcfg.get("yt", "test_data_dir"))
    return test_data_dir


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
    # this should be used in yt.loaders.load

    path = Path(fn).expanduser()

    if path.exists():
        return path

    err_msg = f"No such file or directory: '{fn}'."

    data_dir = _get_test_data_dir_path()
    if data_dir.is_dir():
        alt_path = Path.joinpath(data_dir, fn)
        if alt_path.exists():
            return alt_path
        err_msg += f"\n(Also tried '{alt_path}')."

    raise FileNotFoundError(err_msg)


# create a pooch.Pooch instance
# pardon the poor variable naming: it is common practice to give those dog names...
data_registry = get_data_registry_table()
cache_storage = _get_test_data_dir_path() / "yt_download_cache"
pooch_registry = {k: v["hash"] for k, v in sample_data_registry.items()}
base_url = "https://yt-project.org/data/"
poochie = pooch.create(path=cache_storage, base_url=base_url, registry=pooch_registry)


def _download_sample_data_file(filename, progressbar=True):
    """
    Download a file by url.

    Returns
    -------
    storage_filename : location of the downloaded file

    """

    if progressbar:
        try:
            import tqdm  # noqa
        except ImportError:
            # note that pooch cannot use our vendored version of tqdm
            mylog.error(
                "progress bar support requires tqdm. "
                "Install with `pip install -u tqdm`."
            )
            progressbar = False

    downloader = pooch.HTTPDownloader(progressbar=progressbar)

    poochie.fetch(filename, downloader=downloader)
    storage_filename = Path.joinpath(poochie.path, filename)
    return storage_filename
