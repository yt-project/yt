import os

import numpy as np

from yt.data_objects.static_output import ParticleDataset
from yt.funcs import mylog
from yt.geometry.particle_geometry_handler import ParticleIndex


class SPHDataset(ParticleDataset):
    default_kernel_name = "cubic"
    _sph_smoothing_styles = ["scatter", "gather"]
    _sph_smoothing_style = "scatter"
    _num_neighbors = 32
    _use_sph_normalization = True

    def __init__(
        self,
        filename,
        dataset_type=None,
        units_override=None,
        unit_system="cgs",
        index_order=None,
        index_filename=None,
        kdtree_filename=None,
        kernel_name=None,
        default_species_fields=None,
    ):
        if kernel_name is None:
            self.kernel_name = self.default_kernel_name
        else:
            self.kernel_name = kernel_name
        self.kdtree_filename = kdtree_filename
        super().__init__(
            filename,
            dataset_type=dataset_type,
            units_override=units_override,
            unit_system=unit_system,
            index_order=index_order,
            index_filename=index_filename,
            default_species_fields=default_species_fields,
        )

    @property
    def num_neighbors(self):
        return self._num_neighbors

    @num_neighbors.setter
    def num_neighbors(self, value):
        if value < 0:
            raise ValueError(f"Negative value not allowed: {value}")
        self._num_neighbors = value

    @property
    def sph_smoothing_style(self):
        return self._sph_smoothing_style

    @sph_smoothing_style.setter
    def sph_smoothing_style(self, value):
        if value not in self._sph_smoothing_styles:
            raise ValueError(
                "Smoothing style not implemented: %s, please "
                "select one of the following: " % value,
                self._sph_smoothing_styles,
            )

        self._sph_smoothing_style = value

    @property
    def use_sph_normalization(self):
        return self._use_sph_normalization

    @use_sph_normalization.setter
    def use_sph_normalization(self, value):
        if value is not True and value is not False:
            raise ValueError("SPH normalization needs to be True or False!")
        self._use_sph_normalization = value


class SPHParticleIndex(ParticleIndex):
    def _initialize_index(self):
        ds = self.dataset

        ds._file_hash = self._generate_hash()

        if hasattr(self.io, "_generate_smoothing_length"):
            self.io._generate_smoothing_length(self)

        super()._initialize_index()

    def _generate_kdtree(self, fname):
        from yt.utilities.lib.cykdtree import PyKDTree

        if fname is not None:
            if os.path.exists(fname):
                mylog.info("Loading KDTree from %s", os.path.basename(fname))
                kdtree = PyKDTree.from_file(fname)
                if kdtree.data_version != self.ds._file_hash:
                    mylog.info("Detected hash mismatch, regenerating KDTree")
                else:
                    self._kdtree = kdtree
                    return
        positions = []
        for data_file in self.data_files:
            for _, ppos in self.io._yield_coordinates(
                data_file, needed_ptype=self.ds._sph_ptypes[0]
            ):
                positions.append(ppos)
        if positions == []:
            self._kdtree = None
            return
        positions = np.concatenate(positions)
        mylog.info("Allocating KDTree for %s particles", positions.shape[0])
        num_neighbors = getattr(self.ds, "num_neighbors", 32)
        self._kdtree = PyKDTree(
            positions.astype("float64"),
            left_edge=self.ds.domain_left_edge,
            right_edge=self.ds.domain_right_edge,
            periodic=np.array(self.ds.periodicity),
            leafsize=2 * int(num_neighbors),
            data_version=self.ds._file_hash,
        )
        if fname is not None:
            self._kdtree.save(fname)

    @property
    def kdtree(self):
        if hasattr(self, "_kdtree"):
            return self._kdtree

        ds = self.ds

        if getattr(ds, "kdtree_filename", None) is None:
            if os.path.exists(ds.parameter_filename):
                fname = ds.parameter_filename + ".kdtree"
            else:
                # we don't want to write to disk for in-memory data
                fname = None
        else:
            fname = ds.kdtree_filename

        self._generate_kdtree(fname)

        return self._kdtree
