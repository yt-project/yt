from yt.config import ytcfg
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _dask as dask

# a couple of warning messages used in several places
here_be_dragons = (
    "You have chosen to attempt multi-threading. Many frontends are not"
    " thread-safe and many yt functions are built for single-threading so "
    "you will likely encounter bugs (including silent failures and hangs)."
    "Good luck!"
)

thread_override = (
    "multi-threading is prohibited by default. Supply allow_threads=True "
    "to start_client() to override if you really want to (you probably "
    "dont want to)."
)


class ClientContainer:
    """
    a container to help with managing dask clients for yt. Enforces
    single-threading for new clients and allows us to easily track if a client
    is open.

    Intentionally does NOT spin up a dask client on instantiation so that
    a user is responsible for starting their own client. With no client,
    the wrapper functions in dask_helper assume one single-threaded processor.
    """

    def __init__(self):
        self.client = None

    @property
    def _dask_enabled(self):
        # always check config in case user has changed it dynamically
        return ytcfg.get("yt", "internals", "dask_enabled")

    def start_client(self, *args, allow_threads=False, **kwargs):
        """spins up a new dask.distributed.Client and stores in self.client

        Parameters
        ----------
        *args :
            all *args passed to dask.distributed.Client
        allow_threads : bool
            if False (default), prevents multithreading. Switch to True at your
            own risk (multithreading is not supported but might work for some
            operations if you are lucky).
        **kwargs :
            all **kwargs are passed to dask.distributed.Client with some
            exceptions.

            The following **kwargs have defaults set here that differ from the
            defaults of dask.distributed.Client to enforce singlethreaded, single
            processor execution by default. These kwargs and their defaults are:

            threads_per_worker = 1
            n_workers = 1
            proccesses = False

            see dask.distributed.Client for descriptions of arguments.

        Example Usage
        -------------

        To start a default client:
        from yt.utilities.parallel_tools.dask_helper import dask_client

        dask_client.start_client()
        print(dask_client.client)
        print(dask_client.client.status)
        print(dask_client.client.dashboard_link)

        ds = yt.load_sample("snapshot_033")
        reg = ds.region(ds.domain_center,ds.domain_right_edge*0.25, ds.domain_right_edge*0.75)
        si = reg[('PartType4','Silicon')]
        """

        if self._dask_enabled is False:
            raise RuntimeError(
                "To use dask with yt, please set the dask_enabled "
                "config option to True."
            )

        # some sanitizing and defaults different from dask's normal defaults
        # for dask.distributed.Client:
        tpw = self._sanitize_threading(allow_threads, **kwargs)
        nw = kwargs.pop("n_workers", 1)
        pr = kwargs.pop("processes", False)

        # spin up the Client:
        self.client = dask.distributed.Client(
            *args, threads_per_worker=tpw, n_workers=nw, processes=pr, **kwargs
        )

    def _sanitize_threading(self, allow_threads, **kwargs):
        tpw = kwargs.pop("threads_per_worker", 1)
        if tpw > 1:
            if allow_threads:
                mylog.warning(here_be_dragons)
            else:
                mylog.warning(thread_override)
                tpw = 1
        return tpw

    def close(self):
        if self.client:
            self.client.close()
            self.client = None


dask_client = ClientContainer()


def compute(*args, allow_threads=False, **kwargs):
    """direct wrapper of dask.compute. If there is not an active `dask_client`,
    will enforce single-threaded computation. Primarily used by yt internally,
    but may be public uses once dask objects can be returned.

    Parameters
    ----------
    *args : passed to dask.compute()
    allow_threads : bool
        if False (default), prevents multithreading. Switch to True at your
        own risk (multithreading is not supported but might work for some
        operations if you are lucky).
    **kwargs : passed to dask.compute()

    If no dask_client is active, the `scheduler` kwarg is forced to be "singe-threaded"

    Returns
    -------
    output from dask.compute()

    """

    if not ytcfg.get("yt", "internals", "dask_enabled"):
        return args
    elif dask_client.client is None:
        sch = kwargs.pop("scheduler", "single-threaded")
        safe_sch = ["single-threaded", "processes"]
        if sch not in safe_sch:
            if allow_threads:
                mylog.warning(here_be_dragons)
            else:
                mylog.warning(thread_override)
                sch = "single-threaded"
        return dask.compute(*args, scheduler=sch, **kwargs)
    else:
        # sanitization already taken care of by dask_client.
        return dask.compute(*args, **kwargs)


def is_delayed(obj):
    """checks if obj is a dask array (or a unyt-dask array)"""
    return isinstance(obj, dask.array.Array)


def dask_delay_wrapper(func):
    def call_func(*args, **kwargs):
        # note: the dask_enabled switch is a permanent config option
        if ytcfg.get("yt", "internals", "dask_enabled"):
            # if dask is enabled, we delay
            return dask.delayed(func)(*args, **kwargs)
        else:
            # if dask is not enabled, immediately execute and return result
            return func(*args, **kwargs)

    return call_func
