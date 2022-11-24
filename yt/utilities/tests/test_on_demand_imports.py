import pytest

from yt.utilities.on_demand_imports import OnDemand, safe_import


def test_access_available_module():
    class os_imports(OnDemand):
        @safe_import
        def path(self):
            from os import path

            return path

    _os = os_imports()

    _os.path.join("eggs", "saussage")


def test_access_unavailable_module():
    class Bacon_imports(OnDemand):
        @safe_import
        def spam(self):
            from Bacon import spam

            return spam

    _bacon = Bacon_imports()
    with pytest.raises(
        ImportError,
        match=r"No module named 'Bacon'",
    ) as excinfo:
        _bacon.spam()

    # yt should add information to the original error message
    # but this done slightly differently in Python>=3.11
    # (using exception notes), so we can't just match the error message
    # directly. Instead this implements a Python-version agnostic check
    # that the user-visible error message is what we expect.
    complete_error_message = excinfo.exconly()
    assert complete_error_message == (
        "ModuleNotFoundError: No module named 'Bacon'\n"
        "Something went wrong while trying to lazy-import Bacon. "
        "Please make sure that Bacon is properly installed.\n"
        "If the problem persists, please file an issue at "
        "https://github.com/yt-project/yt/issues/new"
    )


def test_class_invalidation():
    with pytest.raises(
        TypeError, match="class .*'s name needs to be suffixed '_imports'"
    ):

        class Bacon(OnDemand):
            pass


def test_base_class_instanciation():
    with pytest.raises(
        TypeError, match="The OnDemand base class cannot be instanciated."
    ):
        OnDemand()
