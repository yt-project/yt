import os
import shutil
import tempfile

import pytest
import yaml

import yt
from yt.config import ytcfg
from yt.utilities.answer_testing import utils
from yt.utilities.logger import ytLogger


def pytest_addoption(parser):
    """
    Lets options be passed to test functions.
    """
    parser.addoption(
        "--with-answer-testing", action="store_true",
    )
    parser.addoption(
        "--answer-store", action="store_true",
    )
    parser.addoption(
        "--answer-big-data", action="store_true",
    )
    parser.addoption(
        "--answer-raw-arrays", action="store_true",
    )
    parser.addoption(
        "--raw-answer-store", action="store_true",
    )
    parser.addoption(
        "--force-overwrite", action="store_true",
    )
    parser.addoption(
        "--no-hash", action="store_true",
    )
    parser.addoption("--local-dir", default=None, help="Where answers are saved.")
    # Tell pytest about the local-dir option in the ini files. This
    # option is used for creating the answer directory on CI
    parser.addini("local-dir", default="../answer-store", help="answer directory.")
    parser.addini(
        "test_data_dir",
        default=ytcfg.get("yt", "test_data_dir"),
        help="Directory where data for tests is stored.",
    )
    parser.addini(
        "answer_file_list",
        default="./tests/tests_pytest.yaml",
        help="Contains answer file names",
    )


def pytest_configure(config):
    r"""
    Reads in the tests/tests.yaml file. This file contains a list of
    each answer test's answer file (including the changeset number).
    """
    yt._called_from_pytest = True
    # Read the list of answer test classes and their associated answer file
    with open(config.getini("answer_file_list"), "r") as f:
        pytest.answer_files = yaml.safe_load(f)
    # Register custom marks for answer tests and big data
    config.addinivalue_line("markers", "answer_test: Run the answer tests.")
    config.addinivalue_line(
        "markers", "big_data: Run answer tests that require" " large data files."
    )


def pytest_collection_modifyitems(config, items):
    r"""
    Decide which tests to skip based on command-line options.
    """
    # Set up the skip marks
    skip_answer = pytest.mark.skip(reason="--with-answer-testing not set.")
    skip_unit = pytest.mark.skip(reason="Running answer tests, so skipping unit tests.")
    skip_big = pytest.mark.skip(reason="--answer-big-data not set.")
    # Loop over every collected test function
    for item in items:
        # If it's an answer test and the appropriate CL option hasn't
        # been set, skip it
        if "answer_test" in item.keywords and not config.getoption(
            "--with-answer-testing"
        ):
            item.add_marker(skip_answer)
        # If it's an answer test that requires big data and the CL
        # option hasn't been set, skip it
        if "big_data" in item.keywords and not config.getoption("--with-answer-testing") and not config.getoption("--answer-big-data"):
            item.add_marker(skip_big)
        if "answer_test" not in item.keywords and config.getoption(
            "--with-answer-testing"
        ):
            item.add_marker(skip_unit)


def _param_list(request):
    r"""
    Saves the non-ds, non-fixture function arguments for saving to
    the answer file.
    """
    # pytest treats parameterized arguments as fixtures, so there's no
    # clean way to separate them out from other other fixtures (that I
    # know of), so we do it explicitly
    blacklist = ["hashing", "answer_file", "request", "answer_compare", "temp_dir"]
    test_params = {}
    for key, val in request.node.funcargs.items():
        # if key not in blacklist and not key.startswith('ds'):
        if key not in blacklist:
            test_params[key] = str(val)
    # Convert python-specific data objects (such as tuples) to a more
    # io-friendly format (in order to not have python-specific anchors
    # in the answer yaml file)
    test_params = utils._streamline_for_io(test_params)
    return test_params


def _get_answer_files(request):
    """
    Gets the path to where the hashed and raw answers are saved.
    """
    try:
        answer_file, raw_answer_file = pytest.answer_files[request.cls.__name__]
    except KeyError:
        ytLogger.error(f"Answer files for `{request.cls.__name__}` not found.")
        sys.exit()
    except ValueError:
        ytLogger.error(
            "Missing either hashed file name or raw file name for `{request.cls.__name__}`in `{request.config.getini('answer_file_list')}`."
        )
    # Add the local-dir aspect of the path. If there's a command line value,
    # have that override the ini file value
    clLocalDir = request.config.getoption("--local-dir")
    iniLocalDir = request.config.getini("local-dir")
    if clLocalDir is not None:
        answer_file = os.path.join(os.path.expanduser(clLocalDir), answer_file)
        raw_answer_file = os.path.join(os.path.expanduser(clLocalDir), raw_answer_file)
    else:
        answer_file = os.path.join(os.path.expanduser(iniLocalDir), answer_file)
        raw_answer_file = os.path.join(os.path.expanduser(iniLocalDir), raw_answer_file)
    # Make sure we don't overwrite unless we mean to
    overwrite = request.config.getoption("--force-overwrite")
    storing = request.config.getoption("--answer-store")
    raw_storing = request.config.getoption("--raw-answer-store")
    raw = request.config.getoption("--answer-raw-arrays")
    if os.path.exists(answer_file) and storing and not overwrite:
        raise FileExistsError(
            "Use `--force-overwrite` to overwrite an existing answer file."
        )
    if os.path.exists(raw_answer_file) and raw_storing and raw and not overwrite:
        raise FileExistsError(
            "Use `--force-overwrite` to overwrite an existing raw answer file."
        )
    # If we do mean to overwrite, do so here by deleting the original file
    if os.path.exists(answer_file) and storing and overwrite:
        os.remove(answer_file)
    if os.path.exists(raw_answer_file) and raw_storing and raw and overwrite:
        os.remove(raw_answer_file)
    print(os.path.abspath(answer_file))
    return answer_file, raw_answer_file


@pytest.fixture(scope="function")
def hashing(request):
    r"""
    Handles initialization, generation, and saving of answer test
    result hashes.
    """
    no_hash = request.config.getoption("--no-hash")
    store_hash = request.config.getoption("--answer-store")
    raw = request.config.getoption("--answer-raw-arrays")
    raw_store = request.config.getoption("--raw-answer-store")
    # This check is so that, when checking if the answer file exists in _get_answer_files,
    # we don't continuously fail. With this check, _get_answer_files is called once
    # per class, despite this having function scope
    if request.cls.answer_file is None:
        request.cls.answer_file, request.cls.raw_answer_file = _get_answer_files(
            request
        )
    request.cls.hashes = {}
    # Load the saved answers if we're comparing. We don't do this for the raw
    # answers because those are huge
    if not no_hash and not store_hash and request.cls.saved_hashes is None:
        with open(request.cls.answer_file, "r") as fd:
            request.cls.saved_hashes = yaml.safe_load(fd)
    yield
    # Get arguments and their values passed to the test (e.g., axis, field, etc.)
    params = _param_list(request)
    # Hash the test results. Don't save to request.cls.hashes so we still have
    # raw data, in case we want to work with that
    hashes = utils._hash_results(request.cls.hashes)
    # Add the other test parameters
    hashes.update(params)
    # Add the function name as the "master" key to the hashes dict
    hashes = {request.node.name: hashes}
    # Save hashes
    if not no_hash and store_hash:
        utils._save_result(hashes, request.cls.answer_file)
    # Compare hashes
    elif not no_hash and not store_hash:
        utils._compare_result(hashes, request.cls.saved_hashes)
    # Save raw data
    if raw and raw_store:
        utils._save_raw_arrays(
            request.cls.hashes, request.cls.raw_answer_file, request.node.name
        )
    # Compare raw data. This is done one test at a time because the
    # arrays can get quite large and storing everything in memory would
    # be bad
    if raw and not raw_store:
        utils._compare_raw_arrays(
            request.cls.hashes, request.cls.raw_answer_file, request.node.name
        )


@pytest.fixture(scope="function")
def temp_dir():
    r"""
    Creates a temporary directory needed by certain tests.
    """
    curdir = os.getcwd()
    if int(os.environ.get("GENERATE_YTDATA", 0)):
        tmpdir = os.getcwd()
    else:
        tmpdir = tempfile.mkdtemp()
    os.chdir(tmpdir)
    yield tmpdir
    os.chdir(curdir)
    if tmpdir != curdir:
        shutil.rmtree(tmpdir)


@pytest.fixture(scope="class")
def ds(request):
    try:
        if isinstance(request.param, list):
            # data_dir_load can take the cls, args, and kwargs. These optional
            # arguments are given in a dictionary as the second element of the
            # list
            ds_fn = request.param[0]
            opts = request.param[1]
            cls = None if "cls" not in opts.keys() else opts["cls"]
            args = None if "args" not in opts.keys() else opts["args"]
            kwargs = None if "kwargs" not in opts.keys() else opts["kwargs"]
            dataset = utils.data_dir_load(ds_fn, cls=cls, args=args, kwargs=kwargs)
        else:
            dataset = utils.data_dir_load(request.param)
        return dataset
    except FileNotFoundError:
        return pytest.skip(f"Data file: `{request.param}` not found.")


@pytest.fixture(scope="class")
def f(request):
    """
    Fixture for returning the field. Needed because indirect=True is
    used for loading the datasets.
    """
    return request.param


@pytest.fixture(scope="class")
def d(request):
    """
    Fixture for returning the ds_obj. Needed because indirect=True is
    used for loading the datasets.
    """
    return request.param


@pytest.fixture(scope="class")
def a(request):
    """
    Fixture for returning the axis. Needed because indirect=True is
    used for loading the datasets.
    """
    return request.param


@pytest.fixture(scope="class")
def w(request):
    """
    Fixture for returning the weight_field. Needed because
    indirect=True is used for loading the datasets.
    """
    return request.param


@pytest.fixture(scope="class")
def ds_repr(request):
    """
    Fixture for returning the string representation of a dataset.
    Needed because indirect=True is used for loading the datasets.
    """
    return request.param


@pytest.fixture(scope="class")
def N(request):
    """
    Fixture for returning the number of particles in a dataset.
    Needed because indirect=True is used for loading the datasets.
    """
    return request.param


@pytest.fixture(scope="class")
def big_data(request):
    """
    There isn't a way to access command-line options given to pytest from
    an actual test, so we need a fixture that retrieves and returns its
    value. In this case, the --answer-big-data option.
    """
    if request.config.getoption("--answer-big-data"):
        return True
    else:
        return False
