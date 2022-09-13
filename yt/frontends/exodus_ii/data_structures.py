import numpy as np

from yt.data_objects.index_subobjects.unstructured_mesh import UnstructuredMesh
from yt.data_objects.static_output import Dataset
from yt.data_objects.unions import MeshUnion
from yt.funcs import setdefaultattr
from yt.geometry.unstructured_mesh_handler import UnstructuredIndex
from yt.utilities.file_handler import NetCDF4FileHandler, warn_netcdf
from yt.utilities.logger import ytLogger as mylog

from .fields import ExodusIIFieldInfo
from .util import get_num_pseudo_dims, load_info_records, sanitize_string


class ExodusIIUnstructuredMesh(UnstructuredMesh):
    _index_offset = 1

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class ExodusIIUnstructuredIndex(UnstructuredIndex):
    def __init__(self, ds, dataset_type="exodus_ii"):
        super().__init__(ds, dataset_type)

    def _initialize_mesh(self):
        coords = self.ds._read_coordinates()
        connectivity = self.ds._read_connectivity()
        self.meshes = []
        for mesh_id, conn_ind in enumerate(connectivity):
            displaced_coords = self.ds._apply_displacement(coords, mesh_id)
            mesh = ExodusIIUnstructuredMesh(
                mesh_id, self.index_filename, conn_ind, displaced_coords, self
            )
            self.meshes.append(mesh)
        self.mesh_union = MeshUnion("mesh_union", self.meshes)

    def _detect_output_fields(self):
        elem_names = self.dataset.parameters["elem_names"]
        node_names = self.dataset.parameters["nod_names"]
        fnames = elem_names + node_names
        self.field_list = []
        for i in range(1, len(self.meshes) + 1):
            self.field_list += [("connect%d" % i, fname) for fname in fnames]
        self.field_list += [("all", fname) for fname in fnames]


class ExodusIIDataset(Dataset):
    _index_class = ExodusIIUnstructuredIndex
    _field_info_class = ExodusIIFieldInfo

    def __init__(
        self,
        filename,
        step=0,
        displacements=None,
        dataset_type="exodus_ii",
        storage_filename=None,
        units_override=None,
    ):
        """

        A class used to represent an on-disk ExodusII dataset. The initializer takes
        two extra optional parameters, "step" and "displacements."

        Parameters
        ----------

        step : integer
            The step tells which time index to slice at. It throws an Error if
            the index is larger than the number of time outputs in the ExodusII
            file. Passing step=-1 picks out the last dataframe.
            Default is 0.

        displacements : dictionary of tuples
            This is a dictionary that controls whether or not displacement fields
            will be used with the meshes in this dataset. The keys of the
            displacements dictionary should the names of meshes in the file
            (e.g., "connect1", "connect2", etc... ), while the values should be
            tuples of the form (scale, offset), where "scale" is a floating point
            value and "offset" is an array-like with one component for each spatial
            dimension in the dataset. When the displacements for a given mesh are
            turned on, the coordinates of the vertices in that mesh get transformed
            as:

                  vertex_x = vertex_x + disp_x*scale + offset_x
                  vertex_y = vertex_y + disp_y*scale + offset_y
                  vertex_z = vertex_z + disp_z*scale + offset_z

            If no displacement
            fields (assumed to be named 'disp_x', 'disp_y', etc... ) are detected in
            the output file, then this dictionary is ignored.

        Examples
        --------

        This will load the Dataset at time index '0' with displacements turned off.

        >>> import yt
        >>> ds = yt.load("MOOSE_sample_data/mps_out.e")

        This will load the Dataset at the final index with displacements turned off.

        >>> import yt
        >>> ds = yt.load("MOOSE_sample_data/mps_out.e", step=-1)

        This will load the Dataset at index 10, turning on displacement fields for
        the 2nd mesh without applying any scale or offset:

        >>> import yt
        >>> ds = yt.load(
        ...     "MOOSE_sample_data/mps_out.e",
        ...     step=10,
        ...     displacements={"connect2": (1.0, [0.0, 0.0, 0.0])},
        ... )

        This will load the Dataset at index 10, scaling the displacements
        in the 2nd mesh by a factor of 5 while not applying an offset:

        >>> import yt
        >>> ds = yt.load(
        ...     "MOOSE_sample_data/mps_out.e",
        ...     step=10,
        ...     displacements={"connect2": (1.0, [0.0, 0.0, 0.0])},
        ... )

        This will load the Dataset at index 10, scaling the displacements for
        the 2nd mesh by a factor of 5.0 and shifting all the vertices in
        the first mesh by 1.0 unit in the z direction.

        >>> import yt
        >>> ds = yt.load(
        ...     "MOOSE_sample_data/mps_out.e",
        ...     step=10,
        ...     displacements={
        ...         "connect1": (0.0, [0.0, 0.0, 1.0]),
        ...         "connect2": (5.0, [0.0, 0.0, 0.0]),
        ...     },
        ... )

        """
        self.step = step
        if displacements is None:
            self.displacements = {}
        else:
            self.displacements = displacements
        self.storage_filename = storage_filename

        super().__init__(filename, dataset_type, units_override=units_override)

        self.fluid_types += self._get_fluid_types()
        self.default_field = [f for f in self.field_list if f[0] == "connect1"][-1]

    @property
    def index_filename(self):
        # historic alias
        return self.filename

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        setdefaultattr(self, "length_unit", self.quan(1.0, "cm"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "g"))
        setdefaultattr(self, "time_unit", self.quan(1.0, "s"))
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")

    def _parse_parameter_file(self):
        self._handle = NetCDF4FileHandler(self.parameter_filename)
        with self._handle.open_ds() as ds:
            self._read_glo_var()
            self.dimensionality = ds.variables["coor_names"].shape[0]
            self.parameters["info_records"] = self._load_info_records()
            self.num_steps = len(ds.variables["time_whole"])
            self.current_time = self._get_current_time()
            self.parameters["num_meshes"] = ds.variables["eb_status"].shape[0]
            self.parameters["elem_names"] = self._get_elem_names()
            self.parameters["nod_names"] = self._get_nod_names()
            self.domain_left_edge, self.domain_right_edge = self._load_domain_edge()
            self._periodicity = (False, False, False)

        # These attributes don't really make sense for unstructured
        # mesh data, but yt warns if they are not present, so we set
        # them to dummy values here.
        self.domain_dimensions = np.ones(3, "int32")
        self.cosmological_simulation = 0
        self.current_redshift = 0
        self.omega_lambda = 0
        self.omega_matter = 0
        self.hubble_constant = 0
        self.refine_by = 0

    def _get_fluid_types(self):
        with NetCDF4FileHandler(self.parameter_filename).open_ds() as ds:
            fluid_types = ()
            i = 1
            while True:
                ftype = "connect%d" % i
                if ftype in ds.variables:
                    fluid_types += (ftype,)
                    i += 1
                else:
                    break
            fluid_types += ("all",)
            return fluid_types

    def _read_glo_var(self):
        """
        Adds each global variable to the dict of parameters

        """
        names = self._get_glo_names()
        if not names:
            return
        with self._handle.open_ds() as ds:
            values = ds.variables["vals_glo_var"][:].transpose()
            for name, value in zip(names, values):
                self.parameters[name] = value

    def _load_info_records(self):
        """
        Returns parsed version of the info_records.
        """
        with self._handle.open_ds() as ds:
            try:
                return load_info_records(ds.variables["info_records"])
            except (KeyError, TypeError):
                mylog.warning("No info_records found")
                return []

    def _get_current_time(self):
        with self._handle.open_ds() as ds:
            try:
                return ds.variables["time_whole"][self.step]
            except IndexError as e:
                raise RuntimeError(
                    "Invalid step number, max is %d" % (self.num_steps - 1)
                ) from e
            except (KeyError, TypeError):
                return 0.0

    def _get_glo_names(self):
        """

        Returns the names of the global vars, if available.

        """

        with self._handle.open_ds() as ds:
            if "name_glo_var" not in ds.variables:
                mylog.warning("name_glo_var not found")
                return []
            else:
                return [
                    sanitize_string(v.tobytes()) for v in ds.variables["name_glo_var"]
                ]

    def _get_elem_names(self):
        """

        Returns the names of the element vars, if available.

        """

        with self._handle.open_ds() as ds:
            if "name_elem_var" not in ds.variables:
                mylog.warning("name_elem_var not found")
                return []
            else:
                return [
                    sanitize_string(v.tobytes()) for v in ds.variables["name_elem_var"]
                ]

    def _get_nod_names(self):
        """

        Returns the names of the node vars, if available

        """

        with self._handle.open_ds() as ds:
            if "name_nod_var" not in ds.variables:
                mylog.warning("name_nod_var not found")
                return []
            else:
                return [
                    sanitize_string(v.tobytes()) for v in ds.variables["name_nod_var"]
                ]

    def _read_coordinates(self):
        """

        Loads the coordinates for the mesh

        """

        coord_axes = "xyz"[: self.dimensionality]

        mylog.info("Loading coordinates")
        with self._handle.open_ds() as ds:
            if "coord" not in ds.variables:
                coords = (
                    np.array([ds.variables[f"coord{ax}"][:] for ax in coord_axes])
                    .transpose()
                    .astype("f8")
                )
            else:
                coords = (
                    np.array([coord for coord in ds.variables["coord"][:]])
                    .transpose()
                    .astype("f8")
                )
            return coords

    def _apply_displacement(self, coords, mesh_id):

        mesh_name = "connect%d" % (mesh_id + 1)
        if mesh_name not in self.displacements:
            new_coords = coords.copy()
            return new_coords

        new_coords = np.zeros_like(coords)
        fac = self.displacements[mesh_name][0]
        offset = self.displacements[mesh_name][1]

        coord_axes = "xyz"[: self.dimensionality]
        with self._handle.open_ds() as ds:
            for i, ax in enumerate(coord_axes):
                if f"disp_{ax}" in self.parameters["nod_names"]:
                    ind = self.parameters["nod_names"].index(f"disp_{ax}")
                    disp = ds.variables["vals_nod_var%d" % (ind + 1)][self.step]
                    new_coords[:, i] = coords[:, i] + fac * disp + offset[i]

            return new_coords

    def _read_connectivity(self):
        """
        Loads the connectivity data for the mesh
        """
        mylog.info("Loading connectivity")
        connectivity = []
        with self._handle.open_ds() as ds:
            for i in range(self.parameters["num_meshes"]):
                var = ds.variables["connect%d" % (i + 1)][:].astype("i8")
                try:
                    elem_type = var.elem_type.lower()
                    if elem_type == "nfaced":
                        raise NotImplementedError(
                            "3D arbitrary polyhedra are not implemented yet"
                        )
                    arbitrary_polyhedron = elem_type == "nsided"
                except AttributeError:
                    arbitrary_polyhedron = False

                conn = var[:]
                if arbitrary_polyhedron:
                    nodes_per_element = ds.variables[f"ebepecnt{i + 1}"]
                    npe = nodes_per_element[0]
                    if np.any(nodes_per_element != npe):
                        raise NotImplementedError("only equal-size polyhedra supported")
                    q, r = np.divmod(len(conn), npe)
                    assert r == 0
                    conn.shape = (q, npe)
                connectivity.append(conn)
            return connectivity

    def _load_domain_edge(self):
        """
        Loads the boundaries for the domain edge

        """

        coords = self._read_coordinates()
        connectivity = self._read_connectivity()

        mi = 1e300
        ma = -1e300
        for mesh_id, _ in enumerate(connectivity):
            displaced_coords = self._apply_displacement(coords, mesh_id)
            mi = np.minimum(displaced_coords.min(axis=0), mi)
            ma = np.maximum(displaced_coords.max(axis=0), ma)

        # pad domain boundaries
        width = ma - mi
        mi -= 0.1 * width
        ma += 0.1 * width

        # set up pseudo-3D for lodim datasets here
        for _ in range(self.dimensionality, 3):
            mi = np.append(mi, 0.0)
            ma = np.append(ma, 1.0)

        num_pseudo_dims = get_num_pseudo_dims(coords)
        self.dimensionality -= num_pseudo_dims
        for i in range(self.dimensionality, 3):
            mi[i] = 0.0
            ma[i] = 1.0

        return mi, ma

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs):
        warn_netcdf(filename)
        try:
            from netCDF4 import Dataset

            # We use keepweakref here to avoid holding onto the file handle
            # which can interfere with other is_valid calls.
            with Dataset(filename, keepweakref=True) as f:
                f.variables["connect1"]
            return True
        except Exception:
            pass
