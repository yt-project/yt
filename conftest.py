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

import numpy as np
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
        "--answer-store",
        action="store_true",
        default=False,
    )
    parser.addoption(
        "--answer-big-data",
        action="store_true",
        default=False,
    )
    parser.addoption(
        "--with-answer-testing",
        action="store_true",
        default=False,
    )
    parser.addoption(
        "--with-unit-testing",
        action="store_true",
        default=False,
    )
    parser.addoption(
        "--force-override-answers",
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
        pytest.answer_files = yaml.load(f, Loader=yaml.Loader)
    # Register custom marks for answer tests and big data
    config.addinivalue_line('markers', 'answer_test: Run the answer tests.')
    config.addinivalue_line('markers', 'big_data: Run answer tests that require'
        ' large data files.')


#============================================
#     pytest_collection_modifyitems
#============================================
def pytest_collection_modifyitems(config, items):
    """
    Sets up the marks.
    
    pytest version 5.x.x no longer provides access to pytest.config
    (which was originally used for skipping answer tests when the CL
    option wasn't set), but pytest 5.x.x is what works with the most
    recent pytest-xdist, which enables parallizing by class. This
    function is called once, right after test collection.
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
#                answer_store
#============================================
@pytest.fixture(scope='class')
def cli_testing_opts(request):
    """
    Gets the value of the parameters passed at the command line.
    """
    # When using the class scope for a fixture, request has an
    # attribute called cls, which is the class that's using the
    # fixture. In this case, we can add the options from the parser as
    # attributes to the class using the fixture. i.e., we can now do
    # self.answer_store from the class using this fixture
    if request.cls is not None:
        request.cls.answer_store = request.config.getoption("--answer-store")
        request.cls.save_dir = answer_dir 


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
    if request.cls is not None:
        if request.cls.__name__ in pytest.answer_files:
            request.cls.answer_file = pytest.answer_files[request.cls.__name__]
            # Make sure we're not overwriting an existing answer set
            if os.path.isfile(os.path.join(answer_dir, request.cls.answer_file)):
                if request.config.getoption('--answer-store'):
                    raise FileExistsError("Error, attempting to overwrite "
                        "answer file {}. Either specify a new version or "
                        "set the `--force-override-answers` option".format(
                            request.cls.answer_file
                        )
                    )
        else:
            assert False


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
    test_params = {}
    # Yield to the caller in order to actually perform the tests
    yield
    # Get the test parameters
    func = request.node.function
    # co_varnames is all of the variable names local to the function
    # it starts with self, then the passed args, then the vars defined
    # in the function body. This excludes fixture names
    args = func.__code__.co_varnames[1:func.__code__.co_argcount]
    # funcargs includes the names and values of all arguments, including
    # fixtures, so we use args to weed out the fixtures. Need to have
    # special treatment of the data files loaded in fixtures for the
    # frontends
    for key, val in request.node.funcargs.items():
        if key in args and not key.startswith('ds_'):
            test_params[key] = val
    request.cls.hashes.update(test_params)
    # Hash the arrays 
    hashes = utils.array_to_hash(request.cls.hashes)
    # Finalize by adding the function name as a key
    hashes = {request.node.name : hashes}
    # Either save or compare
    utils.handle_hashes(request.cls.save_dir, request.cls.answer_file, hashes,
        request.cls.answer_store)
