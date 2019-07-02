"""
Title:   conftest.py
Author:  Jared Coughlin
Date:    6/21/19
Purpose: Contains hooks and fixtures for yt testing
Notes:   https://docs.pytest.org/en/latest/example/simple.html
"""
import pytest


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
        request.cls.run_big_data = request.config.getoption("--answer-big-data")
