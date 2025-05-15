import atexit
import time
from bisect import insort
from collections import defaultdict
from datetime import datetime as dt
from functools import wraps

from yt.config import ytcfg
from yt.funcs import mylog


class PerformanceCounters:
    _shared_state = {}  # type: ignore

    def __new__(cls, *args, **kwargs):
        self = object.__new__(cls, *args, **kwargs)
        self.__dict__ = cls._shared_state
        return self

    def __init__(self):
        self.counters = defaultdict(lambda: 0.0)
        self.counting = defaultdict(lambda: False)
        self.starttime = defaultdict(lambda: 0)
        self.endtime = defaultdict(lambda: 0)
        self._on = ytcfg.get("yt", "time_functions")
        self.exit()

    def __call__(self, name):
        if not self._on:
            return
        if self.counting[name]:
            self.counters[name] = time.time() - self.counters[name]
            self.counting[name] = False
            self.endtime[name] = dt.now()
        else:
            self.counters[name] = time.time()
            self.counting[name] = True
            self.starttime[name] = dt.now()

    def call_func(self, func):
        if not self._on:
            return func

        @wraps(func)
        def func_wrapper(*args, **kwargs):
            self(func.__name__)
            func(*args, **kwargs)
            self(func.__name__)

        return func_wrapper

    def print_stats(self):
        mylog.info("Current counter status:\n")
        times = []
        for i in self.counters:
            insort(times, [self.starttime[i], i, 1])  # 1 for 'on'
            if not self.counting[i]:
                insort(times, [self.endtime[i], i, 0])  # 0 for 'off'
        shifts = {}
        order = []
        endtimes = {}
        shift = 0
        multi = 5
        for i in times:
            # a starting entry
            if i[2] == 1:
                shifts[i[1]] = shift
                order.append(i[1])
                shift += 1
            if i[2] == 0:
                shift -= 1
                endtimes[i[1]] = self.counters[i[1]]
        line_fragments: list[str] = []
        for i in order:
            line_fragments.append(
                f"{' ' * shifts[i] * multi}{shifts[i]} : {i} : "
                f"{'still running' if self.counting[i] else f'{self.counters[i]:0.3e}'}\n"
            )
        mylog.info("\n%s", "".join(line_fragments))

    def exit(self):
        if self._on:
            atexit.register(self.print_stats)


yt_counters = PerformanceCounters()
time_function = yt_counters.call_func


class ProfilingController:
    def __init__(self):
        self.profilers = {}

    def profile_function(self, function_name):
        def wrapper(func):
            try:
                import cProfile
            except ImportError:
                return func
            my_prof = cProfile.Profile()
            self.profilers[function_name] = my_prof

            @wraps(func)
            def run_in_profiler(*args, **kwargs):
                my_prof.enable()
                func(*args, **kwargs)
                my_prof.disable()

            return run_in_profiler

        return wrapper

    def write_out(self, filename_prefix):
        pfn = str(filename_prefix)
        if ytcfg.get("yt", "internals", "parallel"):
            global_parallel_rank = ytcfg.get("yt", "internals", "global_parallel_rank")
            global_parallel_size = ytcfg.get("yt", "internals", "global_parallel_size")
            pfn += f"_{global_parallel_rank:03}_{global_parallel_size:03}"

        for n, p in sorted(self.profilers.items()):
            fn = f"{pfn}_{n}.cprof"
            mylog.info("Dumping %s into %s", n, fn)
            p.dump_stats(fn)
