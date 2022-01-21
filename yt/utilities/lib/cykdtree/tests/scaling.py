r"""Routines for tracking the scaling of the triangulation routines."""
import cProfile
import os
import pstats
import time

import numpy as np

from yt.utilities.lib.cykdtree.tests import run_test


def stats_run(
    npart,
    nproc,
    ndim,
    periodic=False,
    overwrite=False,
    display=False,
    suppress_final_output=False,
):
    r"""Get timing stats using :package:`cProfile`.

    Args:
        npart (int): Number of particles.
        nproc (int): Number of processors.
        ndim (int): Number of dimensions.
        periodic (bool, optional): If True, the domain is assumed to be
            periodic. Defaults to False.
        overwrite (bool, optional): If True, the existing file for this
            set of input parameters if overwritten. Defaults to False.
        suppress_final_output (bool, optional): If True, the final output
            from spawned MPI processes is suppressed. This is mainly for
            timing purposes. Defaults to False.
        display (bool, optional): If True, display the profile results.
            Defaults to False.

    """
    perstr = ""
    outstr = ""
    if periodic:
        perstr = "_periodic"
    if suppress_final_output:
        outstr = "_noout"
    fname_stat = f"stat_{npart}part_{nproc}proc_{ndim}dim{perstr}{outstr}.txt"
    if overwrite or not os.path.isfile(fname_stat):
        cProfile.run(
            "from yt.utilities.lib.cykdtree.tests import run_test; "
            + f"run_test({npart}, {ndim}, nproc={nproc}, "
            + f"periodic={periodic}, "
            + f"suppress_final_output={suppress_final_output})",
            fname_stat,
        )
    if display:
        p = pstats.Stats(fname_stat)
        p.sort_stats("time").print_stats(10)
        return p
    return fname_stat


def time_run(
    npart, nproc, ndim, nrep=1, periodic=False, leafsize=10, suppress_final_output=False
):
    r"""Get running times using :package:`time`.

    Args:
        npart (int): Number of particles.
        nproc (int): Number of processors.
        ndim (int): Number of dimensions.
        nrep (int, optional): Number of times the run should be performed to
            get an average. Defaults to 1.
        periodic (bool, optional): If True, the domain is assumed to be
            periodic. Defaults to False.
        leafsize (int, optional): The maximum number of points that should be
            in any leaf in the tree. Defaults to 10.
        suppress_final_output (bool, optional): If True, the final output
            from spawned MPI processes is suppressed. This is mainly for
            timing purposes. Defaults to False.

    """
    times = np.empty(nrep, "float")
    for i in range(nrep):
        t1 = time.time()
        run_test(
            npart,
            ndim,
            nproc=nproc,
            periodic=periodic,
            leafsize=leafsize,
            suppress_final_output=suppress_final_output,
        )
        t2 = time.time()
        times[i] = t2 - t1
    return np.mean(times), np.std(times)


def strong_scaling(
    npart=1e6,
    nrep=1,
    periodic=False,
    leafsize=10,
    overwrite=True,
    suppress_final_output=False,
):
    r"""Plot the scaling with number of processors for a particular function.

    Args:
        npart (int, optional): Number of particles. Defaults to 1e6.
        nrep (int, optional): Number of times the run should be performed to
            get an average. Defaults to 1.
        periodic (bool, optional): If True, the domain is assumed to be
            periodic. Defaults to False.
        leafsize (int, optional): The maximum number of points that should be
            in any leaf in the tree. Defaults to 10.
        overwrite (bool, optional): If True, the existing file for this
            set of input parameters if overwritten. Defaults to False.
        suppress_final_output (bool, optional): If True, the final output
            from spawned MPI processes is suppressed. This is mainly for
            timing purposes. Defaults to False.

    """
    import matplotlib.pyplot as plt

    npart = int(npart)
    perstr = ""
    outstr = ""
    if periodic:
        perstr = "_periodic"
    if suppress_final_output:
        outstr = "_noout"
    fname_plot = "plot_strong_scaling_nproc_{}part{}_{}leafsize{}.png".format(
        npart, perstr, leafsize, outstr
    )
    nproc_list = [1, 2, 4, 8]  # , 16]
    ndim_list = [2, 3, 4]
    clr_list = ["b", "r", "g", "m"]
    times = np.empty((len(nproc_list), len(ndim_list), 2), "float")
    for j, nproc in enumerate(nproc_list):
        for i, ndim in enumerate(ndim_list):
            times[j, i, 0], times[j, i, 1] = time_run(
                npart,
                nproc,
                ndim,
                nrep=nrep,
                periodic=periodic,
                leafsize=leafsize,
                suppress_final_output=suppress_final_output,
            )
            print(f"Finished {ndim}D on {nproc}.")
    fig, axs = plt.subplots(1, 1)
    for i in range(len(ndim_list)):
        ndim = ndim_list[i]
        clr = clr_list[i]
        axs.errorbar(
            nproc_list,
            times[:, i, 0],
            yerr=times[:, i, 1],
            fmt=clr,
            label=f"ndim = {ndim}",
        )
    axs.set_xlabel("# of Processors")
    axs.set_ylabel("Time (s)")
    axs.legend()
    fig.savefig(fname_plot)
    print("    " + fname_plot)


def weak_scaling(
    npart=1e4,
    nrep=1,
    periodic=False,
    leafsize=10,
    overwrite=True,
    suppress_final_output=False,
):
    r"""Plot the scaling with number of processors with a constant number of
    particles per processor for a particular function.

    Args:
        npart (int, optional): Number of particles per processor. Defaults to
            1e4.
        nrep (int, optional): Number of times the run should be performed to
            get an average. Defaults to 1.
        periodic (bool, optional): If True, the domain is assumed to be
            periodic. Defaults to False.
        leafsize (int, optional): The maximum number of points that should be
            in any leaf in the tree. Defaults to 10.
        overwrite (bool, optional): If True, the existing file for this
            set of input parameters if overwritten. Defaults to False.
        suppress_final_output (bool, optional): If True, the final output
            from spawned MPI processes is suppressed. This is mainly for
            timing purposes. Defaults to False.

    """
    import matplotlib.pyplot as plt

    npart = int(npart)
    perstr = ""
    outstr = ""
    if periodic:
        perstr = "_periodic"
    if suppress_final_output:
        outstr = "_noout"
    fname_plot = "plot_weak_scaling_nproc_{}part{}_{}leafsize{}.png".format(
        npart, perstr, leafsize, outstr
    )
    nproc_list = [1, 2, 4, 8, 16]
    ndim_list = [2, 3]
    clr_list = ["b", "r", "g", "m"]
    times = np.empty((len(nproc_list), len(ndim_list), 2), "float")
    for j, nproc in enumerate(nproc_list):
        for i, ndim in enumerate(ndim_list):
            times[j, i, 0], times[j, i, 1] = time_run(
                npart * nproc,
                nproc,
                ndim,
                nrep=nrep,
                periodic=periodic,
                leafsize=leafsize,
                suppress_final_output=suppress_final_output,
            )
    fig, axs = plt.subplots(1, 1)
    for i in range(len(ndim_list)):
        ndim = ndim_list[i]
        clr = clr_list[i]
        axs.errorbar(
            nproc_list,
            times[:, i, 0],
            yerr=times[:, i, 1],
            fmt=clr,
            label=f"ndim = {ndim}",
        )
    axs.set_xlabel("# of Processors")
    axs.set_ylabel("Time (s)")
    axs.legend()
    fig.savefig(fname_plot)
    print("    " + fname_plot)
