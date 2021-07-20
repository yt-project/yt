import cProfile
import itertools
import pstats
import sys
import time
from datetime import datetime
from subprocess import PIPE, Popen

import numpy as np
from nose.tools import nottest


def assert_less_equal(x, y):
    size_match = True
    try:
        xshape = (1,)
        yshape = (1,)
        if isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
            if isinstance(x, np.ndarray):
                xshape = x.shape
            if isinstance(y, np.ndarray):
                yshape = y.shape
            size_match = xshape == yshape
            assert (x <= y).all()
        else:
            assert x <= y
    except:
        if not size_match:
            raise AssertionError(
                "Shape mismatch\n\n"
                + f"x.shape: {str(x.shape)}\ny.shape: {str(y.shape)}\n"
            )
        raise AssertionError(
            "Variables are not less-equal ordered\n\n" + f"x: {str(x)}\ny: {str(y)}\n"
        )


def call_subprocess(np, func, args, kwargs):
    # Create string with arguments & kwargs
    args_str = ""
    for a in args:
        args_str += f"{a},"
    for k, v in kwargs.items():
        args_str += f"{k}={v},"
    if args_str.endswith(","):
        args_str = args_str[:-1]
    cmd = [
        "mpirun",
        "-n",
        str(np),
        sys.executable,
        "-c",
        f"'from {func.__module__} import {func.__name__}; {func.__name__}({args_str})'",
    ]
    cmd = " ".join(cmd)
    print(f"Running the following command:\n{cmd}")
    p = Popen(cmd, stdin=PIPE, stdout=PIPE, stderr=PIPE, shell=True)
    output, err = p.communicate()
    exit_code = p.returncode
    print(output.decode("utf-8"))
    if exit_code != 0:
        print(err.decode("utf-8"))
        raise Exception("Error on spawned process. See output.")
        return None
    return output.decode("utf-8")


def iter_dict(dicts):
    try:
        return (
            dict(itertools.izip(dicts, x))
            for x in itertools.product(*dicts.itervalues())
        )
    except AttributeError:
        # python 3
        return (dict(zip(dicts, x)) for x in itertools.product(*dicts.values()))


def parametrize(**pargs):
    for k in pargs.keys():
        if not isinstance(pargs[k], (tuple, list)):
            pargs[k] = (pargs[k],)

    def dec(func):
        def pfunc(kwargs0):
            # Wrapper so that name encodes parameters
            def wrapped(*args, **kwargs):
                kwargs.update(**kwargs0)
                return func(*args, **kwargs)

            wrapped.__name__ = func.__name__
            for k, v in kwargs0.items():
                wrapped.__name__ += f"_{k}{v}"
            return wrapped

        def func_param(*args, **kwargs):
            out = []
            for ipargs in iter_dict(pargs):
                out.append(pfunc(ipargs)(*args, **kwargs))
            return out

        func_param.__name__ = func.__name__

        return func_param

    return dec


np.random.seed(100)
pts2 = np.random.rand(100, 2).astype("float64")
pts3 = np.random.rand(100, 3).astype("float64")
rand_state = np.random.get_state()
left_neighbors_x = [[], [0], [1], [2], [], [], [4, 5], [5]]  # None  # None  # None
left_neighbors_y = [
    [],  # None
    [],  # None
    [],  # None
    [],  # None
    [0, 1],
    [4],
    [1, 2, 3],
    [6],
]
left_neighbors_x_periodic = [[3], [0], [1], [2], [6], [6, 7], [4, 5], [5]]
left_neighbors_y_periodic = [[5], [5, 7], [7], [7], [0, 1], [4], [1, 2, 3], [6]]


@nottest
def make_points_neighbors(periodic=False):
    ndim = 2
    npts = 50
    leafsize = 10
    np.random.set_state(rand_state)
    pts = np.random.rand(npts, ndim).astype("float64")
    left_edge = np.zeros(ndim, "float64")
    right_edge = np.ones(ndim, "float64")
    if periodic:
        lx = left_neighbors_x_periodic
        ly = left_neighbors_y_periodic
    else:
        lx = left_neighbors_x
        ly = left_neighbors_y
    num_leaves = len(lx)
    ln = [lx, ly]
    rn = [[[] for i in range(num_leaves)] for _ in range(ndim)]
    for d in range(ndim):
        for i in range(num_leaves):
            for j in ln[d][i]:
                rn[d][j].append(i)
        for i in range(num_leaves):
            rn[d][i] = list(set(rn[d][i]))
    return pts, left_edge, right_edge, leafsize, ln, rn


@nottest
def make_points(npts, ndim, leafsize=10, distrib="rand", seed=100):
    ndim = int(ndim)
    npts = int(npts)
    leafsize = int(leafsize)
    np.random.seed(seed)
    LE = 0.0
    RE = 1.0
    left_edge = LE * np.ones(ndim, "float64")
    right_edge = RE * np.ones(ndim, "float64")
    if npts <= 0:
        npts = 100
        leafsize = 10
        if ndim == 2:
            pts = pts2
        elif ndim == 3:
            pts = pts3
        else:
            pts = np.random.rand(npts, ndim).astype("float64")
    else:
        if distrib == "rand":
            pts = np.random.rand(npts, ndim).astype("float64")
        elif distrib == "uniform":
            pts = np.random.uniform(low=LE, high=RE, size=(npts, ndim))
        elif distrib in ("gaussian", "normal"):
            pts = np.random.normal(
                loc=(LE + RE) / 2.0, scale=(RE - LE) / 4.0, size=(npts, ndim)
            )
            np.clip(pts, LE, RE)
        else:
            raise ValueError(f"Invalid 'distrib': {distrib}")
    return pts, left_edge, right_edge, leafsize


@nottest
def run_test(
    npts,
    ndim,
    nproc=0,
    distrib="rand",
    periodic=False,
    leafsize=10,
    profile=False,
    suppress_final_output=False,
    **kwargs,
):
    r"""Run a routine with a designated number of points & dimensions on a
    selected number of processors.

    Args:
        npart (int): Number of particles.
        nproc (int): Number of processors.
        ndim (int): Number of dimensions.
        distrib (str, optional): Distribution that should be used when
            generating points. Defaults to 'rand'.
        periodic (bool, optional): If True, the domain is assumed to be
            periodic. Defaults to False.
        leafsize (int, optional): Maximum number of points that should be in
            an leaf. Defaults to 10.
        profile (bool, optional): If True cProfile is used. Defaults to False.
        suppress_final_output (bool, optional): If True, the final output
            from spawned MPI processes is suppressed. This is mainly for
            timing purposes. Defaults to False.

    """
    from yt.utilities.lib.cykdtree import make_tree

    unique_str = datetime.today().strftime("%Y%j%H%M%S")
    pts, left_edge, right_edge, leafsize = make_points(
        npts, ndim, leafsize=leafsize, distrib=distrib
    )
    # Set keywords for multiprocessing version
    if nproc > 1:
        kwargs["suppress_final_output"] = suppress_final_output
        if profile:
            kwargs["profile"] = f"{unique_str}_mpi_profile.dat"
    # Run
    if profile:
        pr = cProfile.Profile()
        t0 = time.time()
        pr.enable()
    make_tree(
        pts,
        nproc=nproc,
        left_edge=left_edge,
        right_edge=right_edge,
        periodic=periodic,
        leafsize=leafsize,
        **kwargs,
    )
    if profile:
        pr.disable()
        t1 = time.time()
        ps = pstats.Stats(pr)
        ps.add(kwargs["profile"])
        if isinstance(profile, str):
            ps.dump_stats(profile)
            print(f"Stats saved to {profile}")
        else:
            sort_key = "tottime"
            ps.sort_stats(sort_key).print_stats(25)
            # ps.sort_stats(sort_key).print_callers(5)
            print(f"{t1 - t0} s according to 'time'")
        return ps
