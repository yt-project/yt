import json
import os
import stat

import meshio  # should be on demand import: need to use a fork for now...
import numpy as np
import xmltodict  # should be on demand import

from yt.data_objects.index_subobjects.unstructured_mesh import UnstructuredMesh
from yt.data_objects.static_output import Dataset
from yt.data_objects.unions import MeshUnion
from yt.funcs import setdefaultattr
from yt.geometry.unstructured_mesh_handler import UnstructuredIndex
from yt.utilities.logger import ytLogger as mylog

from .fields import ASPECTFieldInfo


class ASPECTUnstructuredMesh(UnstructuredMesh):
    _index_offset = 1

    def __init__(self, *args, **kwargs):
        super(ASPECTUnstructuredMesh, self).__init__(*args, **kwargs)


class ASPECTUnstructuredIndex(UnstructuredIndex):
    """
    for on-disk vtu files. treats each vtu file as a separate mesh.
    """

    def __init__(self, ds, dataset_type="aspect"):
        super(ASPECTUnstructuredIndex, self).__init__(ds, dataset_type)

    def _initialize_mesh(self):
        connectivity, coords = self.ds._read_pieces()
        self.meshes = []
        for mesh_id, conn_ind in enumerate(connectivity):
            # aspect does not displace between vtu
            # displaced_coords = self.ds._apply_displacement(coords, mesh_id)
            mesh_coords = coords[mesh_id]
            mesh_file = self.ds.parameters["vtu_files"][mesh_id]
            mesh = ASPECTUnstructuredMesh(
                mesh_id, mesh_file, conn_ind, mesh_coords, self
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


class ASPECTDataset(Dataset):
    _index_class = ASPECTUnstructuredIndex
    _field_info_class = ASPECTFieldInfo

    def __init__(
        self,
        filename,
        step=0,
        displacements=None,
        dataset_type="aspect",
        storage_filename=None,
        units_override=None,
    ):
        """

        A class used to represent an on-disk ASPECT dataset. The initializer takes
        two extra optional parameters, "step" and "displacements."

        Parameters
        ----------

        step : integer
            The step tells which time index to slice at. It throws an Error if
            the index is larger than the number of time outputs in the ExodusII
            file. Passing step=-1 picks out the last dataframe.
            Default is 0.



        """
        self.parameter_filename = filename
        self.data_dir = os.path.dirname(filename)
        self.fluid_types += self._get_fluid_types()
        self.step = step
        if displacements is None:
            self.displacements = {}
        else:
            self.displacements = displacements
        super(ASPECTDataset, self).__init__(
            filename, dataset_type, units_override=units_override
        )
        self.index_filename = filename
        self.storage_filename = storage_filename
        self.default_field = [f for f in self.field_list if f[0] == "connect1"][-1]
        self.parameter_info = self._read_sidecar()

    def _get_sidecar_file(self):
        param_file = os.path.basename(self.parameter_filename)
        return os.path.join(
            self.data_dir, "." + os.path.splitext(param_file)[0] + ".json"
        )

    def _init_sidecar(self):
        with open(self._get_sidecar_file(), "w") as shandle:
            json.dump([{"pvtu": self.parameter_filename}], shandle)

    def _read_sidecar(self):
        # initialize the hidden sidecar file for this .pvtu if it is not there
        sidecar = self._get_sidecar_file()
        if not os.path.isfile(sidecar):
            self._init_sidecar()

        # read in parameter info
        with open(sidecar, "r") as shandle:
            param_info = json.load(shandle)

        return param_info[0]

    def _write_sidecar(self):
        # writes current parameter_info to json file
        with open(self._get_sidecar_file(), "w") as shandle:
            json.dump([self.parameter_info], shandle)

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        setdefaultattr(self, "length_unit", self.quan(1.0, "m"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "kg"))
        setdefaultattr(self, "time_unit", self.quan(1.0, "s"))
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "m/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")

    def _parse_parameter_file(self):

        # store the top level info
        with open(self.parameter_filename) as pvtu_fi:
            self.parameters["pXML"] = xmltodict.parse(pvtu_fi.read())

        pieces = self.parameters["pXML"]["VTKFile"]["PUnstructuredGrid"]["Piece"]
        self.parameters["vtu_files"] = [piece["@Source"] for piece in pieces]
        self.periodicity = (False, False, False)  # might not be true...
        self.unique_identifier = self._get_unique_identifier()
        self.dimensionality = int(
            self.parameters["pXML"]["VTKFile"]["PUnstructuredGrid"]["PPoints"][
                "PDataArray"
            ]["@NumberOfComponents"]
        )
        self.parameters["nod_names"] = self._get_nod_names()
        self.parameters["elem_names"] = self._get_elem_names()
        self.parameters["num_meshes"] = len(pieces)
        self.current_time = self._get_current_time()
        self.domain_left_edge, self.domain_right_edge = self._load_domain_edge()
        # self.parameters["info_records"] = self._load_info_records()
        # self.num_steps = len(ds.variables["time_whole"])
        self._set_dummy_parameters()

    def _set_dummy_parameters(self):
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
        with open(self.parameter_filename) as pvtu_fi:
            pXML = xmltodict.parse(pvtu_fi.read())
        n_pieces = len(pXML["VTKFile"]["PUnstructuredGrid"]["Piece"])
        fluid_types = tuple([f"connect{pcid}" for pcid in range(n_pieces)])
        fluid_types += ("all",)
        return fluid_types

    def _get_unique_identifier(self):
        return int(os.stat(self.parameter_filename)[stat.ST_CTIME])

    def _get_current_time(self):
        # might do something here...
        return 0.0

    def _get_elem_names(self):
        """

        Returns the names of the element vars, if available.

        """
        # everything right now is nodal data...
        return []

    def _get_nod_names(self):
        """

        Returns the names of the node vars

        """
        fields = self.parameters["pXML"]["VTKFile"]["PUnstructuredGrid"]["PPointData"][
            "PDataArray"
        ]
        return [field["@Name"] for field in fields]

    def _read_piece(self, srcFi, connectivity_offset=0):
        meshPiece = meshio.read(srcFi)  # read it in with meshio
        coords = meshPiece.points.astype(np.float64)  # coords are already global
        # need to generalize the cell type
        cell_types = list(meshPiece.cells_dict.keys())
        if len(cell_types) > 1:
            mylog.warn(
                (
                    "meshio detected more than one cell type, but the "
                    "current aspect vtu loader\n    can only handle a single "
                    f"cell type. Detected cell types: {cell_types}. \n"
                    "    Will continue with the first and hope for the best..."
                )
            )
        self.parameters["cell_type"] = cell_types[0]
        mylog.info(f"detected cell type is {cell_types[0]}.")

        # do an element type check here, reduce it if necessary.

        connectivity = meshPiece.cells_dict[cell_types[0]]  # 2D connectivity array
        # offset the connectivity matrix to global value
        connectivity = np.array(connectivity) + connectivity_offset
        return [connectivity, coords]

    def _read_pieces(self):
        # reads coordinates and connectivity from pieces
        with open(self.parameter_filename) as pvtu_fi:
            pXML = xmltodict.parse(pvtu_fi.read())
        pieces = pXML["VTKFile"]["PUnstructuredGrid"]["Piece"]

        conlist = []  # list of 2D connectivity arrays
        coordlist = []  # global, concatenated coordinate array
        con_offset = -1
        for mesh_id, src in enumerate(pieces):
            # mesh_name = f"connect{mesh_id + 1}"  # connect1, connect2, etc.
            srcfi = os.path.join(
                self.data_dir, src["@Source"]
            )  # full path to .vtu file
            [con, coord] = self._read_piece(srcfi, con_offset + 1)
            # con_offset = con.max()  # reading by file, don't need to concat?
            conlist.append(con.astype("i8"))
            coordlist.append(coord.astype("f8"))

        return conlist, coordlist

    def _load_domain_edge(self):
        """
        Loads the global boundaries for the domain edge

        """

        # check our sidecar file first:
        self._init_sidecar()  # temporary til we know it works...
        self.parameter_info = self._read_sidecar()
        left_edge = self.parameter_info.get("domain_left_edge", None)
        right_edge = self.parameter_info.get("domain_right_edge", None)

        if not left_edge:
            _, coord = self._read_pieces()
            coord = np.vstack(coord)
            left_edge = [coord[:, i].min() for i in range(self.dimensionality)]
            right_edge = [coord[:, i].max() for i in range(self.dimensionality)]
            self.parameter_info["domain_left_edge"] = left_edge
            self.parameter_info["domain_right_edge"] = right_edge
            self._write_sidecar()

        return np.array(left_edge), np.array(right_edge)

    @classmethod
    def _is_valid(self, datafile, *args, **kwargs):

        if datafile.split(".")[-1] == "pvtu":
            try:
                with open(datafile) as data:
                    xml = xmltodict.parse(data.read())
                    if "VTKFile" in xml.keys():
                        return True

            except Exception:
                pass

        return False
