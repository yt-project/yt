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
        match=r"This functionality requires the Bacon package to be installed\.",
    ):
        _bacon.spam()


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
