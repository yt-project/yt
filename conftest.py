import os
import shutil
import sys
import tempfile
from importlib.metadata import version
from importlib.util import find_spec
from pathlib import Path

import pytest
import yaml
from packaging.version import Version

from yt.config import ytcfg
from yt.utilities.answer_testing.testing_utilities import (
    _compare_raw_arrays,
    _hash_results,
    _save_raw_arrays,
    _save_result,
    _streamline_for_io,
    data_dir_load,
)

NUMPY_VERSION = Version(version("numpy"))
PILLOW_VERSION = Version(version("Pillow"))
MATPLOTLIB_VERSION = Version(version("matplotlib"))

# setuptools does not ship with the standard lib starting in Python 3.12, so we need to
# be resilient if it's not available at runtime
if find_spec("setuptools") is not None:
    SETUPTOOLS_VERSION = Version(version("setuptools"))
else:
    SETUPTOOLS_VERSION = None

if find_spec("pandas") is not None:
    PANDAS_VERSION = Version(version("pandas"))
else:
    PANDAS_VERSION = None


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
    parser.addoption("--local-dir", default=None, help="Where answers are saved.")
    # Tell pytest about the local-dir option in the ini files. This
    # option is used for creating the answer directory on CI
    parser.addini(
        "local-dir",
        default=str(Path(__file__).parent / "answer-store"),
        help="answer directory.",
    )
    parser.addini(
        "test_data_dir",
        default=ytcfg.get("yt", "test_data_dir"),
        help="Directory where data for tests is stored.",
    )


def pytest_configure(config):
    r"""
    Reads in the tests/tests.yaml file. This file contains a list of
    each answer test's answer file (including the changeset number).
    """
    # Register custom marks for answer tests and big data
    config.addinivalue_line("markers", "answer_test: Run the answer tests.")
    config.addinivalue_line(
        "markers", "big_data: Run answer tests that require large data files."
    )
    for value in (
        # treat most warnings as errors
        "error",
        # >>> warnings emitted by testing frameworks, or in testing contexts
        # we still have some yield-based tests, awaiting for transition into pytest
        "ignore::pytest.PytestCollectionWarning",
        # matplotlib warnings related to the Agg backend which is used in CI, not much we can do about it
        "ignore:Matplotlib is currently using agg, which is a non-GUI backend, so cannot show the figure.:UserWarning",
        r"ignore:tight_layout.+falling back to Agg renderer:UserWarning",
        #
        # >>> warnings from wrong values passed to numpy
        # these should normally be curated out of the test suite but they are too numerous
        # to deal with in a reasonable time at the moment.
        "ignore:invalid value encountered in log10:RuntimeWarning",
        "ignore:divide by zero encountered in log10:RuntimeWarning",
        #
        # >>> there are many places in yt (most notably at the frontend level)
        # where we open files but never explicitly close them
        # Although this is in general bad practice, it can be intentional and
        # justified in contexts where reading speeds should be optimized.
        # It is not clear at the time of writing how to approach this,
        # so I'm going to ignore this class of warnings altogether for now.
        "ignore:unclosed file.*:ResourceWarning",
    ):
        config.addinivalue_line("filterwarnings", value)

    if SETUPTOOLS_VERSION is not None and SETUPTOOLS_VERSION >= Version("67.3.0"):
        # may be triggered by multiple dependencies
        # see https://github.com/glue-viz/glue/issues/2364
        # see https://github.com/matplotlib/matplotlib/issues/25244
        config.addinivalue_line(
            "filterwarnings",
            r"ignore:(Deprecated call to `pkg_resources\.declare_namespace\('.*'\)`\.\n)?"
            r"Implementing implicit namespace packages \(as specified in PEP 420\) "
            r"is preferred to `pkg_resources\.declare_namespace`\.:DeprecationWarning",
        )

    if SETUPTOOLS_VERSION is not None and SETUPTOOLS_VERSION >= Version("67.5.0"):
        # may be triggered by multiple dependencies
        # see https://github.com/glue-viz/glue/issues/2364
        # see https://github.com/matplotlib/matplotlib/issues/25244
        config.addinivalue_line(
            "filterwarnings",
            "ignore:pkg_resources is deprecated as an API:DeprecationWarning",
        )

    if NUMPY_VERSION >= Version("1.25"):
        if find_spec("h5py") is not None and (
            Version(version("h5py")) < Version("3.9")
        ):
            # https://github.com/h5py/h5py/pull/2242
            config.addinivalue_line(
                "filterwarnings",
                "ignore:`product` is deprecated as of NumPy 1.25.0:DeprecationWarning",
            )

    if PILLOW_VERSION >= Version("11.3.0") and MATPLOTLIB_VERSION <= Version("3.10.3"):
        # patched upstream: https://github.com/matplotlib/matplotlib/pull/30221
        config.addinivalue_line(
            "filterwarnings",
            r"ignore:'mode' parameter is deprecated:DeprecationWarning",
        )

    if PANDAS_VERSION is not None and PANDAS_VERSION >= Version("2.2.0"):
        config.addinivalue_line(
            "filterwarnings",
            r"ignore:\s*Pyarrow will become a required dependency of pandas:DeprecationWarning",
        )

    if sys.version_info >= (3, 12):
        # already patched (but not released) upstream:
        # https://github.com/dateutil/dateutil/pull/1285
        config.addinivalue_line(
            "filterwarnings",
            r"ignore:datetime\.datetime\.utcfromtimestamp\(\) is deprecated:DeprecationWarning",
        )

        if find_spec("ratarmount"):
            # On Python 3.12+, there is a deprecation warning when calling os.fork()
            # in a multi-threaded process. We use this mechanism to mount archives.
            config.addinivalue_line(
                "filterwarnings",
                r"ignore:This process \(pid=\d+\) is multi-threaded, use of fork\(\) "
                r"may lead to deadlocks in the child\."
                ":DeprecationWarning",
            )

    if find_spec("datatree"):
        # the cf_radial dependency arm-pyart<=1.9.2 installs the now deprecated
        # xarray-datatree package (which imports as datatree), which triggers
        # a bunch of runtimewarnings when importing xarray.
        # https://github.com/yt-project/yt/pull/5042#issuecomment-2457797694
        config.addinivalue_line(
            "filterwarnings",
            "ignore:" r"Engine.*loading failed.*" ":RuntimeWarning",
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
        if (
            "big_data" in item.keywords
            and not config.getoption("--with-answer-testing")
            and not config.getoption("--answer-big-data")
        ):
            item.add_marker(skip_big)
        if "answer_test" not in item.keywords and config.getoption(
            "--with-answer-testing"
        ):
            item.add_marker(skip_unit)


def pytest_itemcollected(item):
    # Customize pytest-mpl decorator to add sensible defaults

    mpl_marker = item.get_closest_marker("mpl_image_compare")
    if mpl_marker is not None:
        # in a future version, pytest-mpl may gain an option for doing this:
        # https://github.com/matplotlib/pytest-mpl/pull/181
        mpl_marker.kwargs.setdefault("tolerance", 0.5)


def _param_list(request):
    r"""
    Saves the non-ds, non-fixture function arguments for saving to
    the answer file.
    """
    # pytest treats parameterized arguments as fixtures, so there's no
    # clean way to separate them out from other other fixtures (that I
    # know of), so we do it explicitly
    blacklist = [
        "hashing",
        "answer_file",
        "request",
        "answer_compare",
        "temp_dir",
        "orbit_traj",
        "etc_traj",
    ]
    test_params = {}
    for key, val in request.node.funcargs.items():
        if key not in blacklist:
            # For plotwindow, the callback arg is a tuple and the second
            # element contains a memory address, so we need to drop it.
            # The first element is the callback name, which is all that's
            # needed
            if key == "callback":
                val = val[0]
            test_params[key] = str(val)
    # Convert python-specific data objects (such as tuples) to a more
    # io-friendly format (in order to not have python-specific anchors
    # in the answer yaml file)
    test_params = _streamline_for_io(test_params)
    return test_params


def _get_answer_files(request):
    """
    Gets the path to where the hashed and raw answers are saved.
    """
    answer_file = f"{request.cls.__name__}_{request.cls.answer_version}.yaml"
    raw_answer_file = f"{request.cls.__name__}_{request.cls.answer_version}.h5"
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
    # This check is so that, when checking if the answer file exists in
    # _get_answer_files, we don't continuously fail. With this check,
    # _get_answer_files is called once per class, despite this having function
    # scope
    if request.cls.answer_file is None:
        request.cls.answer_file, request.cls.raw_answer_file = _get_answer_files(
            request
        )
    if not no_hash and not store_hash and request.cls.saved_hashes is None:
        try:
            with open(request.cls.answer_file) as fd:
                request.cls.saved_hashes = yaml.safe_load(fd)
        except FileNotFoundError:
            module_filename = f"{request.function.__module__.replace('.', os.sep)}.py"
            with open(f"generate_test_{os.getpid()}.txt", "a") as fp:
                fp.write(f"{module_filename}::{request.cls.__name__}\n")
            pytest.fail(msg="Answer file not found.", pytrace=False)
    request.cls.hashes = {}
    # Load the saved answers if we're comparing. We don't do this for the raw
    # answers because those are huge
    yield
    # Get arguments and their values passed to the test (e.g., axis, field, etc.)
    params = _param_list(request)
    # Hash the test results. Don't save to request.cls.hashes so we still have
    # raw data, in case we want to work with that
    hashes = _hash_results(request.cls.hashes)
    # Add the other test parameters
    hashes.update(params)
    # Add the function name as the "master" key to the hashes dict
    hashes = {request.node.name: hashes}
    # Save hashes
    if not no_hash and store_hash:
        _save_result(hashes, request.cls.answer_file)
    # Compare hashes
    elif not no_hash and not store_hash:
        try:
            for test_name, test_hash in hashes.items():
                assert test_name in request.cls.saved_hashes
                assert test_hash == request.cls.saved_hashes[test_name]
        except AssertionError:
            pytest.fail(f"Comparison failure: {request.node.name}", pytrace=False)
    # Save raw data
    if raw and raw_store:
        _save_raw_arrays(
            request.cls.hashes, request.cls.raw_answer_file, request.node.name
        )
    # Compare raw data. This is done one test at a time because the
    # arrays can get quite large and storing everything in memory would
    # be bad
    if raw and not raw_store:
        _compare_raw_arrays(
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
    # data_dir_load can take the cls, args, and kwargs. These optional
    # arguments, if present,  are given in a dictionary as the second
    # element of the list
    if isinstance(request.param, str):
        ds_fn = request.param
        opts = {}
    else:
        ds_fn, opts = request.param
    try:
        return data_dir_load(
            ds_fn, cls=opts.get("cls"), args=opts.get("args"), kwargs=opts.get("kwargs")
        )
    except FileNotFoundError:
        return pytest.skip(f"Data file: `{request.param}` not found.")


@pytest.fixture(scope="class")
def field(request):
    """
    Fixture for returning the field. Needed because indirect=True is
    used for loading the datasets.
    """
    return request.param


@pytest.fixture(scope="class")
def dobj(request):
    """
    Fixture for returning the ds_obj. Needed because indirect=True is
    used for loading the datasets.
    """
    return request.param


@pytest.fixture(scope="class")
def axis(request):
    """
    Fixture for returning the axis. Needed because indirect=True is
    used for loading the datasets.
    """
    return request.param


@pytest.fixture(scope="class")
def weight(request):
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
def Npart(request):
    """
    Fixture for returning the number of particles in a dataset.
    Needed because indirect=True is used for loading the datasets.
    """
    return request.param
