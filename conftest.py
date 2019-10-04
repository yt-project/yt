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
from collections import OrderedDict
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
        pytest.answer_files = yaml.load(f)


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
#+===========================================
@pytest.fixture(scope='class')
def answer_file(request):
    if request.cls is not None:
        if request.cls.__name__ in pytest.answer_files:
            request.cls.answer_file = pytest.answer_files[request.cls.__name__]
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
    """
    # Set up hashes dictionary
    if request.cls is not None:
        request.cls.hashes = OrderedDict()
    else:
        assert False
    # Yield to the caller in order to actually perform the tests
    yield
    # Hash the arrays 
    import pdb; pdb.set_trace()
    hashes = utils.array_to_hash(request.cls.hashes)
    # Finalize by adding the function name as a key
    hashes = {request.function.__name__ : hashes}
    # Either save or compare
    utils.handle_hashes(request.cls.save_dir, request.cls.answer_file, hashes,
        request.cls.answer_store)
