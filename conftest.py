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

from yt.config import ytcfg
from yt.utilities.answer_testing import utils


# Global variables can be added to the pytest namespace
pytest.answer_files = {}

# List of answer files
answer_file_list = 'tests/tests.yaml'
answer_dir = os.path.join(ytcfg.get('yt', 'test_data_dir'), 'answers')


#============================================
#             pytest_addoption
#============================================
def pytest_addoption(parser):
    """
    Lets options be passed to test functions.
    """
    parser.addoption(
        "--with-answer-testing",
        action="store_true",
        default=False,
    )
    parser.addoption(
        "--answer-store",
        action="store_true",
        default=False,
    )
    parser.addoption(
        "--answer-big-data",
        action="store_true",
        default=False,
    )


#============================================
#             pytest_configure
#============================================
def pytest_configure(config):
    """
    Reads in the tests/tests.yaml file. This file contains a list of
    each answer test's answer file (including the changeset number).
    """
    # Make sure that the answers dir exists. If not, try to make it
    if not os.path.isdir(answer_dir):
        os.mkdir(answer_dir)
    # Read the list of answer test classes and their associated answer
    # file
    with open(answer_file_list, 'r') as f:
        pytest.answer_files = yaml.safe_load(f)
    # Register custom marks for answer tests and big data
    config.addinivalue_line('markers', 'answer_test: Run the answer tests.')
    config.addinivalue_line('markers', 'big_data: Run answer tests that require'
        ' large data files.')


#============================================
#     pytest_collection_modifyitems
#============================================
def pytest_collection_modifyitems(config, items):
    """
    Decide which tests to skip based on command-line options.
    """
    # Set up the skip marks
    skip_answer = pytest.mark.skip(reason="--with-answer-testing not set.")
    skip_big = pytest.mark.skip(reason="--answer-big-data not set.")
    # Loop over every collected test function
    for item in items:
        # If it's an answer test and the appropriate CL option hasn't
        # been set, skip it
        if "answer_test" in item.keywords and not config.getoption("--with-answer-testing"):
            item.add_marker(skip_answer)
        # If it's an answer test that requires big data and the CL
        # option hasn't been set, skip it
        if "big_data" in item.keywords and not config.getoption("--answer-big-data"):
            item.add_marker(skip_big)


#============================================
#                   temp_dir
#============================================
@pytest.fixture(scope='function')
def temp_dir():
    """
    Creates a temporary directory needed by certain tests.
    """
    curdir = os.getcwd()
    if int(os.environ.get('GENERATE_YTDATA', 0)):
        tmpdir = os.getcwd()
    else:
        tmpdir = tempfile.mkdtemp()
    os.chdir(tmpdir)
    yield tmpdir
    os.chdir(curdir)
    if tmpdir != curdir:
        shutil.rmtree(tmpdir)


#============================================
#                answer_file
#============================================
@pytest.fixture(scope='class')
def answer_file(request):
    if request.cls.__name__ in pytest.answer_files:
        answer_file = pytest.answer_files[request.cls.__name__]
        # Make sure we're not overwriting an existing answer set
        if os.path.isfile(os.path.join(answer_dir, answer_file)):
            if request.config.getoption('--answer-store'):
                raise FileExistsError("Error, attempting to overwrite "
                    "answer file {}. Either specify a new version or "
                    "set the `--force-override-answers` option".format(
                        answer_file
                    )
                )
    else:
        assert False
    request.cls.answer_file = answer_file


#============================================
#                param_list
#============================================
def param_list(request):
    """
    Saves the non-ds, non-fixture function arguments for saving to
    the answer file.
    """
    test_params = {}
    func = request.node.function
    # co_varnames is all of the variable names local to the function
    # starting with self, then the passed args, then the vars defined
    # in the function body. This excludes fixture names
    args = func.__code__.co_varnames[1:func.__code__.co_argcount]
    # funcargs includes the names and values of all arguments, including
    # fixtures, so we use args to weed out the fixtures. Need to have
    # special treatment of the data files loaded in fixtures for the
    # frontends
    for key, val in request.node.funcargs.items():
        if key in args and not key.startswith('ds_'):
            test_params[key] = val
    # Convert python-specific data objects (such as tuples) to a more
    # io-friendly format (in order to not have python-specific anchors
    # in the answer yaml file)
    test_params = utils.streamline_for_io(test_params)
    return test_params


#============================================
#                  hashing
#============================================
@pytest.fixture(scope='function')
def hashing(request):
    """
    This fixture reduces answer test boilerplate by handling the
    initialization of the hashes, the actual hashing of the arrays
    returned by the tests, and performing the writing/comparison.
    It also handles saving the values and names of the test parameters.
    """
    # Set up hashes dictionary
    if request.cls is not None:
        request.cls.hashes = {}
    else:
        assert False
    # Yield to the caller in order to actually perform the tests
    yield
    # Get param list
    params = param_list(request)
    # Hash the test results
    hashes = utils.hash_results(request.cls.hashes)
    # Add the other test parameters
    hashes.update(params)
    # Add the function name as the "master" key to the hashes dict
    hashes = {request.node.name : hashes}
    # Either save or compare
    utils.handle_hashes(answer_dir, request.cls.answer_file, hashes,
        request.config.getoption('--answer-store'))
