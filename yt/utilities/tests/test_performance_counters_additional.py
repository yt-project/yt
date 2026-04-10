import builtins
from unittest.mock import patch

from numpy.testing import assert_allclose

import yt.utilities.performance_counters as pc


def test_performance_counters_disabled_is_noop_and_returns_original_function():
    with patch.object(pc.ytcfg, "get", side_effect=lambda *args: False), patch.object(
        pc.atexit, "register"
    ) as register_mock:
        counters = pc.PerformanceCounters()

        def identity(value):
            return value

        wrapped = counters.call_func(identity)
        counters("solve")

        assert wrapped is identity
        assert register_mock.call_count == 0
        assert list(counters.counters) == []


def test_performance_counters_record_elapsed_time_and_print_status():
    class _FakeDatetime:
        values = iter([10, 20])

        @classmethod
        def now(cls):
            return next(cls.values)

    with patch.object(
        pc.ytcfg, "get", side_effect=lambda *args: True
    ), patch.object(pc.atexit, "register") as register_mock, patch.object(
        pc.time, "time", side_effect=[10.0, 13.5]
    ), patch.object(
        pc, "dt", _FakeDatetime
    ):
        counters = pc.PerformanceCounters()
        counters("solve")
        counters("solve")

        assert register_mock.call_count == 1
        assert register_mock.call_args[0][0] == counters.print_stats
        assert counters.counting["solve"] is False
        assert counters.starttime["solve"] == 10
        assert counters.endtime["solve"] == 20
        assert_allclose(counters.counters["solve"], 3.5)

        counters.counters["running"] = 99.0
        counters.counting["running"] = True
        counters.starttime["running"] = 30

        with patch.object(pc.mylog, "info") as info_mock:
            counters.print_stats()

        assert info_mock.call_args_list[0].args == ("Current counter status:\n",)
        rendered = info_mock.call_args_list[1].args[1]
        assert "solve" in rendered
        assert "3.500e+00" in rendered
        assert "running" in rendered
        assert "still running" in rendered


def test_time_function_wrapper_tracks_calls_when_enabled():
    class _FakeDatetime:
        values = iter([100, 200])

        @classmethod
        def now(cls):
            return next(cls.values)

    recorded = []

    with patch.object(
        pc.ytcfg, "get", side_effect=lambda *args: True
    ), patch.object(pc.atexit, "register"), patch.object(
        pc.time, "time", side_effect=[1.0, 3.0]
    ), patch.object(
        pc, "dt", _FakeDatetime
    ):
        counters = pc.PerformanceCounters()

        def measure(value):
            recorded.append(value)
            return value * 2

        wrapped = counters.call_func(measure)
        result = wrapped(4)

        assert result is None
        assert recorded == [4]
        assert_allclose(counters.counters["measure"], 2.0)


def test_profile_function_returns_original_when_cprofile_is_unavailable():
    controller = pc.ProfilingController()
    original_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "cProfile":
            raise ImportError("missing cProfile")
        return original_import(name, *args, **kwargs)

    with patch("builtins.__import__", side_effect=fake_import):
        def sentinel():
            return 5

        wrapped = controller.profile_function("sentinel")(sentinel)

    assert wrapped is sentinel


def test_profile_function_and_write_out_cover_parallel_and_non_parallel_names():
    controller = pc.ProfilingController()
    called = []

    @controller.profile_function("profiled")
    def profiled(value):
        called.append(value)
        return value + 1

    result = profiled(7)

    assert result is None
    assert called == [7]
    assert "profiled" in controller.profilers

    dumped = []

    class _DummyProfiler:
        def __init__(self, name):
            self.name = name

        def dump_stats(self, filename):
            dumped.append((self.name, filename))

    controller.profilers = {"beta": _DummyProfiler("beta"), "alpha": _DummyProfiler("alpha")}

    def fake_get_nonparallel(*args):
        values = {
            ("yt", "internals", "parallel"): False,
        }
        return values[args]

    with patch.object(pc.ytcfg, "get", side_effect=fake_get_nonparallel), patch.object(
        pc.mylog, "info"
    ) as info_mock:
        controller.write_out("profile")

    assert dumped == [
        ("alpha", "profile_alpha.cprof"),
        ("beta", "profile_beta.cprof"),
    ]
    assert info_mock.call_count == 2

    dumped.clear()
    controller.profilers = {"run": _DummyProfiler("run")}

    def fake_get_parallel(*args):
        values = {
            ("yt", "internals", "parallel"): True,
            ("yt", "internals", "global_parallel_rank"): 2,
            ("yt", "internals", "global_parallel_size"): 8,
        }
        return values[args]

    with patch.object(pc.ytcfg, "get", side_effect=fake_get_parallel), patch.object(
        pc.mylog, "info"
    ):
        controller.write_out("profile")

    assert dumped == [("run", "profile_002_008_run.cprof")]
