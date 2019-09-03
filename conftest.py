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

from yt.config import ytcfg


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
        test_dir = ytcfg.get('yt', 'test_data_dir')
        request.cls.save_dir = os.path.join(test_dir, 'answers')


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
