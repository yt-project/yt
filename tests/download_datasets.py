"""
Script to download small datasets on Travis for testing with real datasets

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2018, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

import logging
import os
import shutil
import tarfile

import yaml

from yt.config import ytcfg
from yt.funcs import download_file
from yt.utilities.exceptions import YTConfigNotSet

logging.basicConfig(level=logging.INFO)
log = logging.getLogger("yt-download-dataset")

def get_readable_byte_format(size):
    byte_units = {0: 'B', 1: 'KB', 2: 'MB', 3: 'GB', 4: 'TB'}
    power = 2 ** 10
    n = 0
    while size >= power and n < 4:
        size /= power
        n += 1
    return "{0:.2f} {1}".format(size, byte_units[n])

def extract_archive(src_dir, dest_dir, filename):
    src_file = os.path.join(src_dir, filename)
    shutil.copy2(src_file, dest_dir)
    old_wd = os.getcwd()
    os.chdir(dest_dir)
    with tarfile.open(filename, "r:*") as tar:
        tar.extractall()
        size = 0
        for tarinfo in tar:
            size += tarinfo.size
    os.remove(filename)
    os.chdir(old_wd)
    size = get_readable_byte_format(size)
    response = ("Successfully extracted file: %s of %s to destination dir: "
                "%s\n" % (filename, size, dest_dir))
    return response

def dowload_dataset_util(test_data_cache_dir, test_data_dir, filename):
    test_data_url = ytcfg.get("yt", "test_data_url")
    test_data_url = test_data_url.format(filename=filename)

    # Extracted dataset present in test data directory
    if os.path.exists(os.path.join(test_data_dir, filename)):
        log.info("File: %s already present, not downloading again.\n" % filename)
        return

    # Dataset present in cache but not extracted
    if os.path.exists(os.path.join(test_data_cache_dir, filename)):
        # extract and return
        response = extract_archive(test_data_cache_dir, test_data_dir, filename)
        log.info(response)
        return

    download_filename_path = os.path.join(test_data_cache_dir, filename)
    log.info("Downloading data file using url %s and saving at location: %s"
             % (test_data_url, download_filename_path))
    download_file(url=test_data_url, filename=download_filename_path)
    response = extract_archive(test_data_cache_dir, test_data_dir, filename)
    log.info(response)

def download_datasets():
    test_data_dir = ytcfg.get("yt", "test_data_dir")
    test_data_cache_dir = ytcfg.get("yt", "test_data_cache_dir")

    if test_data_dir == '/does/not/exist':
        raise YTConfigNotSet(var="test_data_dir",
                             extra_msg="Abort downloading datasets.")
    if not test_data_cache_dir:
        raise YTConfigNotSet(var="test_data_cache_dir",
                             extra_msg="Abort downloading datasets.")

    if not os.path.exists(test_data_dir):
        os.makedirs(test_data_dir)
    if not os.path.exists(test_data_cache_dir):
        os.makedirs(test_data_cache_dir)

    # List of datasets to be downloaded and cached
    with open('tests/download_datasets.yaml', 'r') as file:
        download_datasets = yaml.load(file)
        download_datasets = download_datasets['datasets'].split(' ')

    datasets_present = os.listdir(test_data_cache_dir)

    # Cleanup cache directory of unwanted datasets (not listed in yaml)
    for dataset in datasets_present:
        if dataset not in download_datasets:
            os.remove(os.path.join(test_data_cache_dir, dataset))

    # Download and extract datasets
    for dataset in download_datasets:
        log.info("Processing dataset: " + dataset)
        dowload_dataset_util(test_data_cache_dir, test_data_dir, dataset)

if __name__ == "__main__":
    download_datasets()
