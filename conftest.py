"""
Title:   conftest.py
Purpose: Contains hooks and fixtures for yt testing
Notes:
    1.) https://docs.pytest.org/en/latest/example/simple.html
    2.) https://docs.pytest.org/en/latest/historical-notes.html#string-conditions
    3.) https://docs.pytest.org/en/latest/skipping.html#skipif
    4.) https://docs.pytest.org/en/latest/reference.html
"""
import os
import shutil
import tempfile

import pytest
import yaml

import yt
from yt.config import ytcfg
from yt.utilities.answer_testing import utils

# Global variables can be added to the pytest namespace
# if this is not done (i.e., just using answer_files = {}), then
# the dictionary is empty come pytest configuration time and
# setting each class' answer_file attribute
pytest.answer_files = {}

# List of answer files
answer_file_list = 'tests/tests_pytest.yaml'
answer_dir = os.path.join(ytcfg.get('yt', 'test_data_dir'), 'answers')
array_dir = os.path.join(answer_dir, 'raw_arrays')


def pytest_addoption(parser):
    """
    Lets options be passed to test functions.
    """
    parser.addoption(
        "--with-answer-testing",
        action="store_true",
    )
    parser.addoption(
        "--answer-store",
        action="store_true",
    )
    parser.addoption(
        "--answer-big-data",
        action="store_true",
    )
    parser.addoption(
        "--answer-raw-arrays",
        action="store_true",
    )
    parser.addoption(
        "--raw-answer-store",
        action="store_true",
    )
    parser.addoption(
        "--force-overwrite",
        action="store_true",
    )
    parser.addoption(
        "--no-hash",
        action="store_true",
    )


def pytest_configure(config):
    r"""
    Reads in the tests/tests.yaml file. This file contains a list of
    each answer test's answer file (including the changeset number).
    """
    yt._called_from_pytest = True
    # Make sure that the answers dir exists. If not, try to make it
    if not os.path.isdir(answer_dir):
        os.mkdir(answer_dir)
    if not os.path.isdir(array_dir):
        os.mkdir(array_dir)
    # Read the list of answer test classes and their associated answer
    # file
    with open(answer_file_list, 'r') as f:
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
        if "big_data" in item.keywords and not config.getoption("--answer-big-data"):
            item.add_marker(skip_big)


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
def answer_file(request):
    r"""
    Assigns the name of the appropriate answer file as an attribute of
    the calling answer test class.

    The answer file is the file that either already contains previously
    generated answers or else is going to contain newly generated
    answers. If an answer file already exists, then it cannot be used
    to hold newly generated answers; a fresh file must be used instead.

    Each class that performs answer tests (e.g., TestEnzo) has a
    corresponding answer file. These answer files are stored in
    test_data_dir/answers. The names of the files and their
    corresponding answer class are in ./tests/tests.yaml.

    Parameters:
    -----------
        request : pytest.FixtureRequest
            Provides access to the requesting test context. For
            example, if an answer class uses this fixture, such as
            TestEnzo, then request provides access to all of the
            methods and attributes of the TestEnzo class, since
            that class is the user of this fixture (the calling
            context).

    Example:
    --------
        # This fixture should be used whenever a new answer class is
        # defined
        >>> @pytest.mark.usefixtures('answer_file')
        >>> class TestNewFrontend:
        >>>     def test1(self):
                    ...
    """
    # See if the class being tested is one we know about
    try:
        answer_file, raw_answer_file  = pytest.answer_files[request.cls.__name__]
        answer_file = os.path.join(answer_dir, answer_file)
        raw_answer_file = os.path.join(array_dir, raw_answer_file)
    except KeyError:
        raise KeyError("Answer file: `{}` not found in list: `{}`.".format(
            request.config.__name__, answer_file_list))
    except ValueError:
        raise ValueError("Either no hashed answer file name or no raw "
            "answer file name given for: `{}` in: `{}`.".format(
            request.cls.__name__, answer_file_list))
    no_hash = request.config.getoption("--no-hash")
    answer_store = request.config.getoption("--answer-store")
    answer_raw_arrays = request.config.getoption("--answer-raw-arrays")
    raw_answer_store = request.config.getoption("--raw-answer-store")
    force_overwrite = request.config.getoption("--force-overwrite")
    # If we're saving answers, make sure answer file doesn't already exist
    # unless force-overwrite is enabled. Only need to check if that type of
    # answer file is under consideration (e.g., we don't need to check the
    # hash answer files if the user wants to only save raw answer data)
    if not no_hash:
        if answer_store and os.path.isfile(answer_file):
            if not force_overwrite:
                raise FileExistsError("Attempting to overwrite existing answer " +
                "file: `{}`. If you're sure about this, use the `--force-overwrite` " +
                "option.".format(answer_file))
        # If we're comparing new and old answers, first make sure that the
        # file containing the old answers exists so we don't waste time running
        # all the tests only to find we have nothing to compare to
        elif not answer_store and not os.path.isfile(answer_file):
            raise FileNotFoundError("Cannot find `{}` containing gold " +
                "answers.".format(answer_file))
    # Now do the same thing as above but for the raw arrays
    if answer_raw_arrays:
        if raw_answer_store and os.path.isfile(raw_answer_file):
            if not force_overwrite:
                raise FileExistsError("Attempting to overwrite existing answer " +
                "file: `{}`. If you're sure about this, use the `--force-overwrite` " +
                "option.".format(raw_answer_file))
        elif not raw_answer_store and not os.path.isfile(raw_answer_file):
            raise FileNotFoundError("Cannot find `{}` containing gold " +
                "answers.".format(raw_answer_file))
    # If the file exists and we are using the force-overwrite option, then
    # we need to delete the existing file, otherwise it will get appended
    # to and this will mess everything up, including with group name conflicts
    # in h5py when saving raw arrays
    if not no_hash and answer_store and os.path.isfile(answer_file) and force_overwrite:
        os.remove(answer_file)
    if answer_raw_arrays and raw_answer_store and os.path.isfile(raw_answer_file) and force_overwrite:
        os.remove(raw_answer_file)
    request.cls.answer_file = answer_file
    request.cls.raw_answer_file = raw_answer_file


def _param_list(request):
    r"""
    Saves the non-ds, non-fixture function arguments for saving to
    the answer file.
    """
    # pytest treats parameterized arguments as fixtures, so there's no
    # clean way to separate them out from other other fixtures (that I
    # know of), so we do it explicitly
    blacklist = ['hashing', 'answer_file', 'request']
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


@pytest.fixture(scope="function")
def hashing(request):
    r"""
    Handles initialization, generation, and saving of answer test
    result hashes.

    Answer tests that require generated data to be saved to disk
    have that data hashed. This fixture creates an empty dictionary
    as an attribute of the answer class of which the current answer
    test is a method.

    Once the test has been run and the raw data has been saved to this
    hashes dictionary, this fixture hashes the raw data and prepares
    an entry to the answer file containing the test name as well as the
    test parameter names and values to accompany the hash(es).

    These entries are then either compared to an existing entry or
    saved in the new answer file.

    Parameters:
    -----------
        request : pytest.FixtureRequest
            Provides access to the requesting test context. For
            example, if an answer class uses this fixture, such as
            TestEnzo, then request provides access to all of the
            methods and attributes of the TestEnzo class, since
            that class is the user of this fixture (the calling
            context).

    Example:
    --------
        # If every method of an answer class saves data then the
        # fixture can be applied to each method like so:
        >>> @pytest.mark.usefixtures('hashing')
        >>> class TestNewFrontend:
        >>>     def test1(self):
        >>>         ...
        >>>     def test2(self):
        >>>         ...
        # If only certain methods save data, then it must be applied
        # directly to those methods, like so:
        >>> class TestNewFrontend:
        >>>     @pytest.mark.usefixtures('hashing')
        >>>     def test1(self):
        >>>         ...
        >>>     def test2(self):
        >>>         ...
    """
    # Set up hashes dictionary
    if request.cls is not None:
        request.cls.hashes = {}
    else:
        assert False
    # Yield to the caller in order to actually perform the tests
    yield
    # Get param list
    params = _param_list(request)
    # Hash the test results
    hashes = utils._hash_results(request.cls.hashes)
    # Add the other test parameters
    hashes.update(params)
    # Add the function name as the "master" key to the hashes dict
    hashes = {request.node.name: hashes}
    # Either save or compare
    if not request.config.getoption("--no-hash"):
        utils._handle_hashes(request.cls.answer_file, hashes,
            request.config.getoption('--answer-store'))
    if request.config.getoption('--answer-raw-arrays'):
        utils._handle_raw_arrays(request.cls.raw_answer_file,
            request.cls.hashes, request.config.getoption('--raw-answer-store'),
            request.node.name)

@pytest.fixture(scope='class')
def ds(request):
    if isinstance(request.param, list):
        dataset = utils.data_dir_load(request.param[0], kwargs=request.param[1])
    else:
        dataset = utils.data_dir_load(request.param)
    if dataset:
        return dataset
    else:
        pytest.skip(f"Data file: `{request.param}` not found.")

@pytest.fixture(scope='class')
def f(request):
    """
    Fixture for returning the field. Needed because indirect=True is
    used for loading the datasets.
    """
    return request.param

@pytest.fixture(scope='class')
def d(request):
    """
    Fixture for returning the ds_obj. Needed because indirect=True is
    used for loading the datasets.
    """
    return request.param

@pytest.fixture(scope='class')
def a(request):
    """
    Fixture for returning the axis. Needed because indirect=True is
    used for loading the datasets.
    """
    return request.param

@pytest.fixture(scope='class')
def w(request):
    """
    Fixture for returning the weight_field. Needed because
    indirect=True is used for loading the datasets.
    """
    return request.param
