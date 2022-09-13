"""
Chimera data structures



"""

import os
import re

import numpy as np

from yt.data_objects.index_subobjects.unstructured_mesh import SemiStructuredMesh
from yt.data_objects.static_output import Dataset
from yt.geometry.geometry_handler import YTDataChunk
from yt.geometry.unstructured_mesh_handler import UnstructuredIndex
from yt.utilities.file_handler import HDF5FileHandler
from yt.utilities.io_handler import io_registry
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py

from .fields import ChimeraFieldInfo


class ChimeraMesh(SemiStructuredMesh):

    _index_offset = 0

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


def _find_files(filename_c):
    # Returns a list of all files that share a frame number with the input
    dirname, file = os.path.split(filename_c)
    match = re.match(r"chimera_\d+_grid", file)
    if match is None:
        raise RuntimeError(
            rf"Expected filename to be of form 'chimera_\d+_grid_*', got {file!r}"
        )
    prefix = match.group()
    frames = [f for f in os.listdir(dirname) if f.startswith(prefix)]
    index_filenames = [os.path.join(dirname, f) for f in sorted(frames)]
    return index_filenames


class ChimeraUNSIndex(UnstructuredIndex):
    def __init__(self, ds, dataset_type="chimera"):
        self._handle = ds._handle
        super().__init__(ds, dataset_type)
        self.directory = os.path.dirname(self.dataset.filename)
        self.dataset_type = dataset_type

    def _initialize_mesh(self):
        self.meshes = []
        index_filenames = _find_files(
            self.dataset.filename
        )  # Retrieves list of all datafiles with the same frame number
        # Detects Yin-Yang data format
        yy = any("grid_2" in file for file in index_filenames)
        for n, file in enumerate(index_filenames):
            with h5py.File(file, "r") as f:
                nmx, nmy, nmz = tuple(f["mesh"]["array_dimensions"][:])
                l = (
                    int(file[-5:-3]) - 1
                )  # Pulls the subgrid number from the data file name
                if nmz > 2:
                    k = f["fluid"]["entropy"].shape[0]
                    r = f["mesh"]["x_ef"][:-2]
                    theta = f["mesh"]["y_ef"][:]
                    phi = f["mesh"]["z_ef"][
                        k * l : k * (l + 1) + 1
                    ]  # Pulls only the individual subgrid's band of phi values
                elif f["mesh"]["z_ef"][-1] == f["mesh"]["z_ef"][0]:
                    r = f["mesh"]["x_ef"][
                        : f["mesh"]["radial_index_bound"][1]
                        - f["mesh"]["x_ef"].shape[0]
                    ]
                    theta = f["mesh"]["y_ef"][:]
                    phi = np.array([f["mesh"]["z_ef"][0], 2 * np.pi])
                else:
                    r = f["mesh"]["x_ef"][
                        : f["mesh"]["radial_index_bound"][1]
                        - f["mesh"]["x_ef"].shape[0]
                    ]
                    theta = f["mesh"]["y_ef"][:]
                    phi = f["mesh"]["z_ef"][:]

                # Creates variables to hold the size of dimensions
                nxd = r.size
                nyd = theta.size
                nzd = phi.size
                nyzd = nyd * nzd
                nyzd_ = (nyd - 1) * (nzd - 1)

                # Generates and fills coordinate array
                coords = np.zeros((nxd, nyd, nzd, 3), dtype="float64", order="C")
                coords[:, :, :, 0] = r[:, None, None]
                coords[:, :, :, 1] = theta[None, :, None]
                coords[:, :, :, 2] = phi[None, None, :]

                if yy:
                    mylog.warning(
                        "Yin-Yang File Detected; This data is not currently supported."
                    )
                coords.shape = (nxd * nyd * nzd, 3)
                # Connectivity is an array of rows, each of which corresponds to a grid cell.
                # The 8 elements of each row are integers representing the cell verticies.
                # These integers refrence the numerical index of the element of the
                # "coords" array which corresponds to the spatial coordinate.

                connectivity = np.zeros(
                    ((nyd - 1) * (nxd - 1) * (nzd - 1), 8), dtype="int64", order="C"
                )  # Creates scaffold array
                connectivity[0] = [
                    0,
                    1,
                    nzd,
                    (nzd + 1),
                    (nyzd),
                    (nyzd + 1),
                    (nyzd + nzd),
                    (nyzd + nzd + 1),
                ]  # Manually defines first coordinate set

                for p in range(
                    nzd - 1
                ):  # Increments first row around phi to define an arc of cells
                    if p > 0:
                        connectivity[p] = connectivity[p - 1] + 1
                for t in range(
                    nyd - 1
                ):  # Increments this arc around theta to define a shell
                    if t > 0:
                        connectivity[t * (nzd - 1) : (t + 1) * (nzd - 1)] = (
                            connectivity[(t - 1) * (nzd - 1) : t * (nzd - 1)] + nzd
                        )
                for r in range(
                    nxd - 1
                ):  # Increments this shell along r to define a sphere
                    if r > 0:
                        connectivity[r * (nyzd_) : (r + 1) * (nyzd_)] = (
                            connectivity[(r - 1) * (nyzd_) : r * (nyzd_)] + nyzd
                        )

                mesh = ChimeraMesh(
                    n, self.index_filename, connectivity, coords, self
                )  # Creates a mesh object

                if "grid_" in file:
                    mylog.info("Mesh %s generated", (n + 1) / len(index_filenames))
                    self.meshes.append(
                        mesh
                    )  # Adds new mesh to the list of generated meshes

    def _detect_output_fields(self):  # Reads in the available data fields
        with h5py.File(self.index_filename, "r") as f:
            fluids = [
                ("chimera", i)
                for i in f["fluid"]
                if np.shape(f["fluid"][i]) == np.shape(f["fluid"]["rho_c"])
            ]
            abundance = [
                ("chimera", i)
                for i in f["abundance"]
                if np.shape(f["abundance"][i]) == np.shape(f["fluid"]["rho_c"])
            ]
            e_rms = [("chimera", f"e_rms_{i+1}") for i in range(4)]
            lumin = [("chimera", f"lumin_{i+1}") for i in range(4)]
            num_lumin = [("chimera", f"num_lumin_{i+1}") for i in range(4)]
            a_name = [
                ("chimera", i.decode("utf-8").strip()) for i in f["abundance"]["a_name"]
            ]
            self.field_list = (
                fluids
                + abundance
                + e_rms
                + lumin
                + num_lumin
                + [("chimera", "abar")]
                + a_name
            )
            if np.shape(f["abundance"]["nse_c"]) != np.shape(f["fluid"]["rho_c"]):
                self.field_list += [("chimera", "nse_c")]

    def _chunk_io(self, dobj, cache=True, local_only=False):  # Creates Data chunk
        gobjs = getattr(dobj._current_chunk, "objs", dobj._chunk_info)
        for subset in gobjs:
            yield YTDataChunk(
                dobj, "io", [subset], self._count_selection(dobj, [subset]), cache=cache
            )

    def _setup_data_io(self):
        self.io = io_registry[self.dataset_type](self.dataset)


class ChimeraDataset(Dataset):
    _index_class = ChimeraUNSIndex  # ChimeraHierarchy
    _field_info_class = ChimeraFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="chimera",
        storage_filename=None,
        units_override=None,
    ):
        # refinement factor between a grid and its subgrid
        self.refine_by = 1  # Somewhat superfluous for Chimera, but left to avoid errors

        self.fluid_types += ("chimera",)
        super().__init__(filename, dataset_type, units_override=units_override)
        self.storage_filename = storage_filename
        self._handle = HDF5FileHandler(filename)

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        self.length_unit = self.quan(1.0, "cm")
        self.mass_unit = self.quan(1.0, "g")
        self.time_unit = self.quan(1.0, "s")
        self.time_unit = self.quan(1.0, "s")
        self.velocity_unit = self.quan(1.0, "cm/s")
        self.magnetic_unit = self.quan(1.0, "gauss")

    def _parse_parameter_file(self):
        with h5py.File(self.parameter_filename, "r") as f:
            # Reads in simulation time, number of dimensions and shape
            self.current_time = f["mesh"]["time"][()]
            self.dimensionality = 3
            self.domain_dimensions = f["mesh"]["array_dimensions"][()]

            self.geometry = "spherical"  # Uses default spherical geometry
            self._periodicity = (False, False, True)
            dle = [
                f["mesh"]["x_ef"][0],
                f["mesh"]["y_ef"][0],
                f["mesh"]["z_ef"][0],
            ]

            if (
                self.domain_dimensions[2] <= 2
                and f["mesh"]["z_ef"][-1] == f["mesh"]["z_ef"][0]
            ):
                dre = [
                    f["mesh"]["x_ef"][
                        f["mesh"]["radial_index_bound"][1] - f["mesh"]["x_ef"].shape[0]
                    ],
                    f["mesh"]["y_ef"][-1],
                    2 * np.pi,
                ]
            else:
                dre = [
                    f["mesh"]["x_ef"][
                        f["mesh"]["radial_index_bound"][1] - f["mesh"]["x_ef"].shape[0]
                    ],
                    f["mesh"]["y_ef"][-1],
                    f["mesh"]["z_ef"][-1],
                ]
            # Sets left and right bounds based on earlier definitions
            self.domain_right_edge = np.array(dre)
            self.domain_left_edge = np.array(dle)

        self.cosmological_simulation = 0  # Chimera is not a cosmological simulation

    @classmethod
    def _is_valid(self, *args, **kwargs):
        # This accepts a filename or a set of arguments and returns True or
        # False depending on if the file is of the type requested.
        try:
            fileh = HDF5FileHandler(args[0])
            if (
                "fluid" in fileh
                and "agr_c" in fileh["fluid"].keys()
                and "grav_x_c" in fileh["fluid"].keys()
            ):
                return True  # Numpy bless
        except (OSError, ImportError):
            pass
        return False
