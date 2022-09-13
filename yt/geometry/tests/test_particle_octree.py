import os

import numpy as np

import yt.units.dimensions as dimensions
from yt.geometry.oct_container import _ORDER_MAX
from yt.geometry.particle_oct_container import ParticleBitmap, ParticleOctreeContainer
from yt.geometry.selection_routines import RegionSelector
from yt.testing import assert_array_equal, assert_equal, assert_true
from yt.units.unit_registry import UnitRegistry
from yt.units.yt_array import YTArray
from yt.utilities.lib.geometry_utils import (
    get_hilbert_indices,
    get_hilbert_points,
    get_morton_indices,
    get_morton_points,
)

NPART = 32**3
DLE = np.array([0.0, 0.0, 0.0])
DRE = np.array([10.0, 10.0, 10.0])
DW = DRE - DLE
PER = np.array([0, 0, 0], "bool")
dx = DW / (2**_ORDER_MAX)


def test_add_particles_random():
    np.random.seed(int(0x4D3D3D3))
    pos = np.random.normal(0.5, scale=0.05, size=(NPART, 3)) * (DRE - DLE) + DLE
    # Now convert to integers
    for i in range(3):
        np.clip(pos[:, i], DLE[i], DRE[i], pos[:, i])
    # Convert to integers
    pos = np.floor((pos - DLE) / dx).astype("uint64")
    morton = get_morton_indices(pos)
    morton.sort()
    for ndom in [1, 2, 4, 8]:
        octree = ParticleOctreeContainer((1, 1, 1), DLE, DRE)
        octree.n_ref = 32
        for split in np.array_split(morton, ndom):
            octree.add(split)
        octree.finalize()
        # This visits every oct.
        tc = octree.recursively_count()
        total_count = np.zeros(len(tc), dtype="int32")
        for i in sorted(tc):
            total_count[i] = tc[i]
        assert_equal(octree.nocts, total_count.sum())
        # This visits every cell -- including those covered by octs.
        # for dom in range(ndom):
        #    level_count += octree.count_levels(total_count.size-1, dom, mask)
        assert_equal(total_count, [1, 8, 64, 64, 256, 536, 1856, 1672])


class FakeDS:
    domain_left_edge = None
    domain_right_edge = None
    domain_width = None
    unit_registry = UnitRegistry()
    unit_registry.add("code_length", 1.0, dimensions.length)
    periodicity = (False, False, False)


class FakeRegion:
    def __init__(self, nfiles, periodic=False):
        self.ds = FakeDS()
        self.ds.domain_left_edge = YTArray(
            [0.0, 0.0, 0.0], "code_length", registry=self.ds.unit_registry
        )
        self.ds.domain_right_edge = YTArray(
            [nfiles, nfiles, nfiles],
            "code_length",
            registry=self.ds.unit_registry,
            dtype="float64",
        )
        self.ds.domain_width = self.ds.domain_right_edge - self.ds.domain_left_edge
        self.ds.periodicity = (periodic, periodic, periodic)
        self.nfiles = nfiles

    def set_edges(self, file_id, dx=0.1):
        self.left_edge = YTArray(
            [file_id + dx, 0.0, 0.0], "code_length", registry=self.ds.unit_registry
        )
        self.right_edge = YTArray(
            [file_id + 1 - dx, self.nfiles, self.nfiles],
            "code_length",
            registry=self.ds.unit_registry,
        )


class FakeBoxRegion:
    def __init__(self, nfiles, left_edge, right_edge):
        self.ds = FakeDS()
        self.ds.domain_left_edge = YTArray(
            left_edge, "code_length", registry=self.ds.unit_registry
        )
        self.ds.domain_right_edge = YTArray(
            right_edge, "code_length", registry=self.ds.unit_registry
        )
        self.ds.domain_width = self.ds.domain_right_edge - self.ds.domain_left_edge
        self.nfiles = nfiles

    def set_edges(self, center, width):
        self.left_edge = self.ds.domain_left_edge + self.ds.domain_width * (
            center - width / 2
        )
        self.right_edge = self.ds.domain_left_edge + self.ds.domain_width * (
            center + width / 2
        )


def FakeBitmap(
    npart,
    nfiles,
    order1,
    order2,
    left_edge=None,
    right_edge=None,
    periodicity=None,
    decomp="sliced",
    buff=0.1,
    distrib="uniform",
    fname=None,
):
    if left_edge is None:
        left_edge = np.array([0.0, 0.0, 0.0])
    if right_edge is None:
        right_edge = np.array([1.0, 1.0, 1.0])
    if periodicity is None:
        periodicity = np.array([0, 0, 0], "bool")
    reg = ParticleBitmap(
        left_edge, right_edge, periodicity, 12345, nfiles, order1, order2
    )
    # Load from file if it exists
    if isinstance(fname, str) and os.path.isfile(fname):
        reg.load_bitmasks(fname)
    else:
        # Create positions for each file
        posgen = yield_fake_decomp(
            decomp, npart, nfiles, left_edge, right_edge, buff=buff, distrib=distrib
        )
        # Coarse index
        max_npart = 0
        for i, (pos, hsml) in enumerate(posgen):
            max_npart = max(max_npart, pos.shape[0])
            reg._coarse_index_data_file(pos, hsml, i)
            reg._set_coarse_index_data_file(i)
        if i != (nfiles - 1):
            raise RuntimeError(
                f"There are positions for {i + 1} files, but there should be {nfiles}."
            )
        # Refined index
        mask = reg.masks.sum(axis=1).astype("uint8")
        sub_mi1 = np.zeros(max_npart, "uint64")
        sub_mi2 = np.zeros(max_npart, "uint64")
        posgen = yield_fake_decomp(
            decomp, npart, nfiles, left_edge, right_edge, buff=buff, distrib=distrib
        )
        coll = None
        for i, (pos, hsml) in enumerate(posgen):
            nsub_mi, coll = reg._refined_index_data_file(
                coll,
                pos,
                hsml,
                mask,
                sub_mi1,
                sub_mi2,
                i,
                0,
                count_threshold=1,
                mask_threshold=2,
            )
            reg.bitmasks.append(i, coll)
        # Save if file name provided
        if isinstance(fname, str):
            reg.save_bitmasks(fname)
    return reg


def test_bitmap_no_collisions():
    # Test init for slabs of points in x
    left_edge = np.array([0.0, 0.0, 0.0])
    right_edge = np.array([1.0, 1.0, 1.0])
    periodicity = np.array([0, 0, 0], "bool")
    npart = 100
    nfiles = 2
    file_hash = 12345
    order1 = 2
    order2 = 2
    reg = ParticleBitmap(
        left_edge, right_edge, periodicity, file_hash, nfiles, order1, order2
    )
    # Coarse index
    posgen = yield_fake_decomp("sliced", npart, nfiles, left_edge, right_edge)
    max_npart = 0
    for i, (pos, hsml) in enumerate(posgen):
        reg._coarse_index_data_file(pos, hsml, i)
        max_npart = max(max_npart, pos.shape[0])
        reg._set_coarse_index_data_file(i)
        assert_equal(reg.count_total(i), np.sum(reg.masks[:, i]))
    mask = reg.masks.sum(axis=1).astype("uint8")
    ncoll = np.sum(mask > 1)
    nc, nm = reg.find_collisions_coarse()
    assert_equal(nc, 0, "%d coarse collisions" % nc)
    assert_equal(ncoll, nc, "%d in mask, %d in bitmap" % (ncoll, nc))
    # Refined index
    sub_mi1 = np.zeros(max_npart, "uint64")
    sub_mi2 = np.zeros(max_npart, "uint64")
    posgen = yield_fake_decomp("sliced", npart, nfiles, left_edge, right_edge)
    coll = None
    for i, (pos, hsml) in enumerate(posgen):
        nsub_mi, coll = reg._refined_index_data_file(
            coll,
            pos,
            hsml,
            mask,
            sub_mi1,
            sub_mi2,
            i,
            0,
            count_threshold=1,
            mask_threshold=2,
        )
        reg.bitmasks.append(i, coll)
        assert_equal(reg.count_refined(i), 0)
    nr, nm = reg.find_collisions_refined()
    assert_equal(nr, 0, "%d collisions" % nr)


def test_bitmap_collisions():
    # Test init for slabs of points in x
    left_edge = np.array([0.0, 0.0, 0.0])
    right_edge = np.array([1.0, 1.0, 1.0])
    periodicity = np.array([0, 0, 0], "bool")
    nfiles = 2
    file_hash = 12345
    order1 = 2
    order2 = 2
    reg = ParticleBitmap(
        left_edge, right_edge, periodicity, file_hash, nfiles, order1, order2
    )
    # Use same points for all files to force collisions
    pos = cell_centers(order1 + order2, left_edge, right_edge)
    hsml = None
    # Coarse index
    max_npart = 0
    for i in range(nfiles):
        reg._coarse_index_data_file(pos, hsml, i)
        max_npart = max(max_npart, pos.shape[0])
        reg._set_coarse_index_data_file(i)
        assert_equal(reg.count_total(i), np.sum(reg.masks[:, i]))
    mask = reg.masks.sum(axis=1).astype("uint8")
    ncoll = np.sum(mask > 1)
    nc, nm = reg.find_collisions_coarse()
    assert_equal(ncoll, nc, "%d in mask, %d in bitmap" % (ncoll, nc))
    assert_equal(nc, 2 ** (3 * order1), "%d coarse collisions" % nc)
    # Refined index
    sub_mi1 = np.zeros(max_npart, "uint64")
    sub_mi2 = np.zeros(max_npart, "uint64")
    for i in range(nfiles):
        nsub_mi, coll = reg._refined_index_data_file(
            None,
            pos,
            hsml,
            mask,
            sub_mi1,
            sub_mi2,
            i,
            0,
            count_threshold=1,
            mask_threshold=2,
        )
        reg.bitmasks.append(i, coll)
        assert_equal(reg.count_refined(i), ncoll)
    nr, nm = reg.find_collisions_refined()
    assert_equal(nr, 2 ** (3 * (order1 + order2)), "%d collisions" % nr)


def test_bitmap_save_load():
    # Test init for slabs of points in x
    left_edge = np.array([0.0, 0.0, 0.0])
    right_edge = np.array([1.0, 1.0, 1.0])
    periodicity = np.array([0, 0, 0], "bool")
    npart = NPART
    file_hash = 12345
    nfiles = 32
    order1 = 2
    order2 = 2
    fname_fmt = "temp_bitmasks{}.dat"
    i = 0
    fname = fname_fmt.format(i)
    while os.path.isfile(fname):
        i += 1
        fname = fname_fmt.format(i)
    # Create bitmap and save to file
    reg0 = FakeBitmap(npart, nfiles, order1, order2, left_edge, right_edge, periodicity)
    reg0.save_bitmasks(fname)
    # Attempt to load bitmap
    reg1 = ParticleBitmap(
        left_edge, right_edge, periodicity, file_hash, nfiles, order1, order2
    )
    reg1.load_bitmasks(fname)
    assert_true(reg0.iseq_bitmask(reg1))
    # Remove file
    os.remove(fname)


def test_bitmap_select():
    np.random.seed(int(0x4D3D3D3))
    dx = 0.1
    for periodic in [False, True]:
        for nfiles in [2, 15, 31, 32, 33]:
            # Now we create particles
            # Note: we set order1 to log2(nfiles) here for testing purposes to
            # ensure no collisions
            order1 = int(np.ceil(np.log2(nfiles)))  # Ensures zero collisions
            order2 = 2  # No overlap for N = nfiles
            exact_division = nfiles == (1 << order1)
            div = float(nfiles) / float(1 << order1)
            reg = FakeBitmap(
                nfiles**3,
                nfiles,
                order1,
                order2,
                decomp="grid",
                left_edge=np.array([0.0, 0.0, 0.0]),
                right_edge=np.array([nfiles, nfiles, nfiles]),
                periodicity=np.array([periodic, periodic, periodic]),
            )
            # Loop over regions selecting single files
            fr = FakeRegion(nfiles, periodic=periodic)
            for i in range(nfiles):
                fr.set_edges(i, dx)
                selector = RegionSelector(fr)
                (df, gf), (dmask, gmask) = reg.identify_data_files(selector, ngz=1)
                if exact_division:
                    assert_equal(len(df), 1, f"selector {i}, number of files")
                    assert_equal(df[0], i, f"selector {i}, file selected")
                    if periodic and (nfiles != 2):
                        ans_gf = sorted([(i - 1) % nfiles, (i + 1) % nfiles])
                    elif i == 0:
                        ans_gf = [i + 1]
                    elif i == (nfiles - 1):
                        ans_gf = [i - 1]
                    else:
                        ans_gf = [i - 1, i + 1]
                    assert_equal(
                        len(gf),
                        len(ans_gf),
                        f"selector {i}, number of ghost files",
                    )
                    for i in range(len(gf)):
                        assert_equal(gf[i], ans_gf[i], f"selector {i}, ghost files")

                else:
                    lf_frac = np.floor(float(fr.left_edge[0]) / div) * div
                    rf_frac = np.floor(float(fr.right_edge[0]) / div) * div
                    # Selected files
                    lf = int(
                        np.floor(lf_frac)
                        if ((lf_frac % 0.5) == 0)
                        else np.round(lf_frac)
                    )
                    rf = int(
                        np.floor(rf_frac)
                        if ((rf_frac % 0.5) == 0)
                        else np.round(rf_frac)
                    )
                    if (rf + 0.5) >= (rf_frac + div):
                        rf -= 1
                    if (lf + 0.5) <= (lf_frac - div):
                        lf += 1
                    df_ans = np.arange(max(lf, 0), min(rf + 1, nfiles))
                    assert_array_equal(df, df_ans, f"selector {i}, file array")
                    # Ghost zones selected files
                    lf_ghost = int(
                        np.floor(lf_frac - div)
                        if (((lf_frac - div) % 0.5) == 0)
                        else np.round(lf_frac - div)
                    )
                    rf_ghost = int(
                        np.floor(rf_frac + div)
                        if (((rf_frac + div) % 0.5) == 0)
                        else np.round(rf_frac + div)
                    )
                    if not periodic:
                        lf_ghost = max(lf_ghost, 0)
                        rf_ghost = min(rf_ghost, nfiles - 1)
                    if (rf_ghost + 0.5) >= (rf_frac + 2 * div):
                        rf_ghost -= 1
                    gf_ans = []
                    if lf_ghost < lf:
                        gf_ans.append(lf_ghost % nfiles)
                    if rf_ghost > rf:
                        gf_ans.append(rf_ghost % nfiles)
                    gf_ans = np.array(sorted(gf_ans))
                    assert_array_equal(gf, gf_ans, f"selector {i}, ghost file array")


def cell_centers(order, left_edge, right_edge):
    ndim = left_edge.size
    ncells = 2**order
    dx = (right_edge - left_edge) / (2 * ncells)
    d = [
        np.linspace(left_edge[i] + dx[i], right_edge[i] - dx[i], ncells)
        for i in range(ndim)
    ]
    dd = np.meshgrid(*d)
    return np.vstack([x.flatten() for x in dd]).T


def fake_decomp_random(npart, nfiles, ifile, DLE, DRE, buff=0.0):
    np.random.seed(int(0x4D3D3D3) + ifile)
    nPF = int(npart / nfiles)
    nR = npart % nfiles
    if ifile == 0:
        nPF += nR
    pos = np.empty((nPF, 3), "float64")
    for i in range(3):
        pos[:, i] = np.random.uniform(DLE[i], DRE[i], nPF)
    return pos


def fake_decomp_sliced(npart, nfiles, ifile, DLE, DRE, buff=0.0):
    np.random.seed(int(0x4D3D3D3) + ifile)
    DW = DRE - DLE
    div = DW / nfiles
    nPF = int(npart / nfiles)
    nR = npart % nfiles
    inp = nPF
    if ifile == 0:
        inp += nR
    iLE = DLE[0] + ifile * div[0]
    iRE = iLE + div[0]
    if ifile != 0:
        iLE -= buff * div[0]
    if ifile != (nfiles - 1):
        iRE += buff * div[0]
    pos = np.empty((inp, 3), dtype="float")
    pos[:, 0] = np.random.uniform(iLE, iRE, inp)
    for i in range(1, 3):
        pos[:, i] = np.random.uniform(DLE[i], DRE[i], inp)
    return pos


def makeall_decomp_hilbert_gaussian(
    npart,
    nfiles,
    DLE,
    DRE,
    buff=0.0,
    order=6,
    verbose=False,
    fname_base=None,
    nchunk=10,
    width=None,
    center=None,
    frac_random=0.1,
):
    import pickle

    np.random.seed(int(0x4D3D3D3))
    DW = DRE - DLE
    if fname_base is None:
        fname_base = f"hilbert{order}_gaussian_np{npart}_nf{nfiles}_"
    if width is None:
        width = 0.1 * DW
    if center is None:
        center = DLE + 0.5 * DW

    def load_pos(file_id):
        filename = fname_base + f"file{file_id}"
        if os.path.isfile(filename):
            fd = open(filename, "rb")
            positions = pickle.load(fd)
            fd.close()
        else:
            positions = np.empty((0, 3), dtype="float64")
        return positions

    def save_pos(file_id, positions):
        filename = fname_base + f"file{file_id}"
        fd = open(filename, "wb")
        pickle.dump(positions, fd)
        fd.close()

    npart_rnd = int(frac_random * npart)
    npart_gau = npart - npart_rnd
    dim_hilbert = 1 << order
    nH = dim_hilbert**3
    if nH < nfiles:
        raise ValueError("Fewer hilbert cells than files.")
    nHPF = nH / nfiles
    rHPF = nH % nfiles
    for ichunk in range(nchunk):
        inp = npart_gau / nchunk
        if ichunk == 0:
            inp += npart_gau % nchunk
        pos = np.empty((inp, 3), dtype="float64")
        ind = np.empty((inp, 3), dtype="int64")
        for k in range(3):
            pos[:, k] = np.clip(
                np.random.normal(center[k], width[k], inp),
                DLE[k],
                DRE[k] - (1.0e-9) * DW[k],
            )
            ind[:, k] = (pos[:, k] - DLE[k]) / (DW[k] / dim_hilbert)
        harr = get_hilbert_indices(order, ind)
        farr = (harr - rHPF) / nHPF
        for ifile in range(nfiles):
            ipos = load_pos(ifile)
            if ifile == 0:
                idx = farr <= ifile  # Put remainders in first file
            else:
                idx = farr == ifile
            ipos = np.concatenate((ipos, pos[idx, :]), axis=0)
            save_pos(ifile, ipos)
    # Random
    for ifile in range(nfiles):
        ipos = load_pos(ifile)
        ipos_rnd = fake_decomp_hilbert_uniform(
            npart_rnd, nfiles, ifile, DLE, DRE, buff=buff, order=order, verbose=verbose
        )
        ipos = np.concatenate((ipos, ipos_rnd), axis=0)
        save_pos(ifile, ipos)


def fake_decomp_hilbert_gaussian(
    npart, nfiles, ifile, DLE, DRE, buff=0.0, order=6, verbose=False, fname=None
):
    np.random.seed(int(0x4D3D3D3))
    DW = DRE - DLE
    dim_hilbert = 1 << order
    nH = dim_hilbert**3
    if nH < nfiles:
        raise Exception("Fewer hilbert cells than files.")
    nHPF = nH / nfiles
    rHPF = nH % nfiles
    hdiv = DW / dim_hilbert
    if ifile == 0:
        hlist = np.arange(0, nHPF + rHPF, dtype="int64")
    else:
        hlist = np.arange(ifile * nHPF + rHPF, (ifile + 1) * nHPF + rHPF, dtype="int64")
    hpos = get_hilbert_points(order, hlist)
    iLE = np.empty((len(hlist), 3), dtype="float")
    iRE = np.empty((len(hlist), 3), dtype="float")
    count = np.zeros(3, dtype="int64")
    pos = np.empty((npart, 3), dtype="float")
    for k in range(3):
        iLE[:, k] = DLE[k] + hdiv[k] * hpos[:, k]
        iRE[:, k] = iLE[:, k] + hdiv[k]
        iLE[hpos[:, k] != 0, k] -= buff * hdiv[k]
        iRE[hpos[:, k] != (dim_hilbert - 1), k] += buff * hdiv[k]
        gpos = np.clip(
            np.random.normal(DLE[k] + DW[k] / 2.0, DW[k] / 10.0, npart), DLE[k], DRE[k]
        )
        for ipos in gpos:
            for i in range(len(hlist)):
                if iLE[i, k] <= ipos < iRE[i, k]:
                    pos[count[k], k] = ipos
                    count[k] += 1
                    break
    return pos[: count.min(), :]


def fake_decomp_hilbert_uniform(
    npart, nfiles, ifile, DLE, DRE, buff=0.0, order=6, verbose=False
):
    np.random.seed(int(0x4D3D3D3) + ifile)
    DW = DRE - DLE
    dim_hilbert = 1 << order
    nH = dim_hilbert**3
    if nH < nfiles:
        raise Exception("Fewer hilbert cells than files.")
    nHPF = nH / nfiles
    rHPF = nH % nfiles
    nPH = npart / nH
    nRH = npart % nH
    hind = np.arange(nH, dtype="int64")
    hpos = get_hilbert_points(order, hind)
    hdiv = DW / dim_hilbert
    if ifile == 0:
        hlist = range(0, nHPF + rHPF)
        nptot = nPH * len(hlist) + nRH
    else:
        hlist = range(ifile * nHPF + rHPF, (ifile + 1) * nHPF + rHPF)
        nptot = nPH * len(hlist)
    pos = np.empty((nptot, 3), dtype="float")
    pc = 0
    for i in hlist:
        iLE = DLE + hdiv * hpos[i, :]
        iRE = iLE + hdiv
        for k in range(3):  # Don't add buffer past domain bounds
            if hpos[i, k] != 0:
                iLE[k] -= buff * hdiv[k]
            if hpos[i, k] != (dim_hilbert - 1):
                iRE[k] += buff * hdiv[k]
        inp = nPH
        if (ifile == 0) and (i == 0):
            inp += nRH
        for k in range(3):
            pos[pc : (pc + inp), k] = np.random.uniform(iLE[k], iRE[k], inp)
        pc += inp
    return pos


def fake_decomp_morton(
    npart, nfiles, ifile, DLE, DRE, buff=0.0, order=6, verbose=False
):
    np.random.seed(int(0x4D3D3D3) + ifile)
    DW = DRE - DLE
    dim_morton = 1 << order
    nH = dim_morton**3
    if nH < nfiles:
        raise Exception("Fewer morton cells than files.")
    nHPF = nH / nfiles
    rHPF = nH % nfiles
    nPH = npart / nH
    nRH = npart % nH
    hind = np.arange(nH, dtype="uint64")
    hpos = get_morton_points(hind)
    hdiv = DW / dim_morton
    if ifile == 0:
        hlist = range(0, nHPF + rHPF)
        nptot = nPH * len(hlist) + nRH
    else:
        hlist = range(ifile * nHPF + rHPF, (ifile + 1) * nHPF + rHPF)
        nptot = nPH * len(hlist)
    pos = np.empty((nptot, 3), dtype="float")
    pc = 0
    for i in hlist:
        iLE = DLE + hdiv * hpos[i, :]
        iRE = iLE + hdiv
        for k in range(3):  # Don't add buffer past domain bounds
            if hpos[i, k] != 0:
                iLE[k] -= buff * hdiv[k]
            if hpos[i, k] != (dim_morton - 1):
                iRE[k] += buff * hdiv[k]
        inp = nPH
        if (ifile == 0) and (i == 0):
            inp += nRH
        for k in range(3):
            pos[pc : (pc + inp), k] = np.random.uniform(iLE[k], iRE[k], inp)
        pc += inp
    return pos


def fake_decomp_grid(npart, nfiles, ifile, DLE, DRE, buff=0.0, verbose=False):
    # TODO: handle 'remainder' particles
    np.random.seed(int(0x4D3D3D3) + ifile)
    DW = DRE - DLE
    nYZ = int(np.sqrt(npart / nfiles))
    div = DW / nYZ
    Y, Z = np.mgrid[
        DLE[1] + 0.1 * div[1] : DRE[1] - 0.1 * div[1] : nYZ * 1j,
        DLE[2] + 0.1 * div[2] : DRE[2] - 0.1 * div[2] : nYZ * 1j,
    ]
    X = 0.5 * div[0] * np.ones(Y.shape, dtype="float64") + div[0] * ifile
    pos = np.array([X.ravel(), Y.ravel(), Z.ravel()], dtype="float64").transpose()
    return pos


def yield_fake_decomp(decomp, npart, nfiles, DLE, DRE, **kws):
    hsml = None
    for ifile in range(nfiles):
        yield fake_decomp(decomp, npart, nfiles, ifile, DLE, DRE, **kws), hsml


def fake_decomp(
    decomp, npart, nfiles, ifile, DLE, DRE, distrib="uniform", fname=None, **kws
):
    import pickle

    if fname is None and distrib == "gaussian":
        fname = f"{decomp}6_{distrib}_np{npart}_nf{nfiles}_file{ifile}"
    if fname is not None and os.path.isfile(fname):
        fd = open(fname, "rb")
        pos = pickle.load(fd)
        fd.close()
        return pos
    if decomp.startswith("zoom_"):
        zoom_factor = 5
        decomp_zoom = decomp.split("zoom_")[-1]
        zoom_npart = npart / 2
        zoom_rem = npart % 2
        pos1 = fake_decomp(
            decomp_zoom,
            zoom_npart + zoom_rem,
            nfiles,
            ifile,
            DLE,
            DRE,
            distrib=distrib,
            **kws,
        )
        DLE_zoom = DLE + 0.5 * DW * (1.0 - 1.0 / float(zoom_factor))
        DRE_zoom = DLE_zoom + DW / zoom_factor
        pos2 = fake_decomp(
            decomp_zoom,
            zoom_npart,
            nfiles,
            ifile,
            DLE_zoom,
            DRE_zoom,
            distrib=distrib,
            **kws,
        )
        pos = np.concatenate((pos1, pos2), axis=0)
    elif "_" in decomp:
        decomp_list = decomp.split("_")
        decomp_np = npart / len(decomp_list)
        decomp_nr = npart % len(decomp_list)
        pos = np.empty((0, 3), dtype="float")
        for i, idecomp in enumerate(decomp_list):
            inp = decomp_np
            if i == 0:
                inp += decomp_nr
            ipos = fake_decomp(
                idecomp, inp, nfiles, ifile, DLE, DRE, distrib=distrib, **kws
            )
            pos = np.concatenate((pos, ipos), axis=0)
    # A perfect grid, no overlap between files
    elif decomp == "grid":
        pos = fake_decomp_grid(npart, nfiles, ifile, DLE, DRE, **kws)
    # Completely random data set
    elif decomp == "random":
        if distrib == "uniform":
            pos = fake_decomp_random(npart, nfiles, ifile, DLE, DRE, **kws)
        else:
            raise ValueError(
                f"Unsupported value {distrib} for input parameter 'distrib'"
            )
    # Each file contains a slab (part of x domain, all of y/z domain)
    elif decomp == "sliced":
        if distrib == "uniform":
            pos = fake_decomp_sliced(npart, nfiles, ifile, DLE, DRE, **kws)
        else:
            raise ValueError(
                f"Unsupported value {distrib} for input parameter 'distrib'"
            )
    # Particles are assigned to files based on their location on a
    # Peano-Hilbert curve of order 6
    elif decomp.startswith("hilbert"):
        if decomp == "hilbert":
            kws["order"] = 6
        else:
            kws["order"] = int(decomp.split("hilbert")[-1])
        if distrib == "uniform":
            pos = fake_decomp_hilbert_uniform(npart, nfiles, ifile, DLE, DRE, **kws)
        elif distrib == "gaussian":
            makeall_decomp_hilbert_gaussian(
                npart, nfiles, DLE, DRE, fname_base=fname.split("file")[0], **kws
            )
            pos = fake_decomp(
                decomp,
                npart,
                nfiles,
                ifile,
                DLE,
                DRE,
                distrib=distrib,
                fname=fname,
                **kws,
            )
        else:
            raise ValueError(
                f"Unsupported value {distrib} for input parameter 'distrib'"
            )
    # Particles are assigned to files based on their location on a
    # Morton ordered Z-curve of order 6
    elif decomp.startswith("morton"):
        if decomp == "morton":
            kws["order"] = 6
        else:
            kws["order"] = int(decomp.split("morton")[-1])
        if distrib == "uniform":
            pos = fake_decomp_morton(npart, nfiles, ifile, DLE, DRE, **kws)
        else:
            raise ValueError(
                f"Unsupported value {distrib} for input parameter 'distrib'"
            )
    else:
        raise ValueError(f"Unsupported value {decomp} for input parameter 'decomp'")
    # Save
    if fname is not None:
        fd = open(fname, "wb")
        pickle.dump(pos, fd)
        fd.close()
    return pos
