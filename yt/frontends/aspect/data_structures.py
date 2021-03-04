import json
import os
import stat

import numpy as np

from yt.data_objects.index_subobjects.unstructured_mesh import UnstructuredMesh
from yt.data_objects.static_output import Dataset
from yt.data_objects.unions import MeshUnion
from yt.funcs import setdefaultattr
from yt.geometry.unstructured_mesh_handler import UnstructuredIndex
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _xmltodict as xmltodict

from .fields import ASPECTFieldInfo
from .util import decode_piece


class ASPECTUnstructuredMesh(UnstructuredMesh):
    _type_name = "aspect_unstructured_mesh"
    _con_args = (
        "mesh_id",
        "filename",
        "filenames",
        "vtu_offsets",
        "connectivity_indices",
        "connectivity_coords",
        "nodes_per_cell",
    )

    def __init__(
        self,
        mesh_id,
        filename,
        filenames,
        node_offsets,
        element_count,
        connectivity_indices,
        connectivity_coords,
        nodes_per_cell,
        index,
    ):
        super().__init__(
            mesh_id, filename, connectivity_indices, connectivity_coords, index
        )
        self.filenames = filenames
        self.node_offsets = node_offsets
        self.element_count = element_count
        self.nodes_per_cell = nodes_per_cell

        # some error checking


class ASPECTUnstructuredIndex(UnstructuredIndex):
    """
    for on-disk vtu files. assumes we have a single mesh split across vtu files.
    """

    def __init__(self, ds, dataset_type="aspect"):
        super().__init__(ds, dataset_type)

        # for when the mesh pieces are treated as independent meshes
        self.global_connectivity = None
        self.global_coords = None

    def _initialize_mesh(self):
        (
            connectivity,
            coords,
            node_offsets,
            el_count,
            nodes_per_cell,
        ) = self.ds._read_pieces()
        self.meshes = []
        mesh_files = self.ds.parameters["vtu_files"]
        pvtu_file = self.ds.parameter_filename
        mesh_id = 0

        # single global mesh
        mesh = ASPECTUnstructuredMesh(
            mesh_id,
            pvtu_file,
            mesh_files,
            node_offsets,
            el_count,
            connectivity,
            coords,
            nodes_per_cell,
            self,
        )
        self.meshes.append(mesh)
        self.mesh_union = MeshUnion("mesh_union", self.meshes)

    def _detect_output_fields(self):
        elem_names = self.dataset.parameters["elem_names"]
        node_names = self.dataset.parameters["nod_names"]
        fnames = elem_names + node_names

        if "velocity" in fnames:
            fnames.remove("velocity")
            for dim in ["x", "y", "z"]:
                fnames.append("velocity_" + dim)

        self.field_list = []
        for i in range(0, len(self.meshes)):
            self.field_list += [("connect%d" % i, fname) for fname in fnames]
        self.field_list += [("all", fname) for fname in fnames]

    def _get_mesh_union(self):
        """
        returns global connectivity and coordinate arrays across the meshes,
        when the mesh is split up across chunks and treated as independent
        meshes (no longer using this for now...).
        """

        if self.global_connectivity is None:
            indices = np.concatenate(
                [
                    mesh.connectivity_indices + mesh._index_offset
                    for mesh in self.mesh_union
                ]
            )
            coords = np.concatenate(
                [mesh.connectivity_coords for mesh in self.mesh_union]
            )
            self.global_connectivity = indices
            self.global_coords = coords
        offset = 0
        return self.global_connectivity, self.global_coords, offset


class ASPECTDataset(Dataset):
    _index_class = ASPECTUnstructuredIndex
    _field_info_class = ASPECTFieldInfo

    def __init__(
        self,
        filename,
        dataset_type="aspect",
        storage_filename=None,
        units_override=None,
    ):
        """

        A class used to represent a single timestep of a on-disk ASPECT dataset

        """
        self.parameter_filename = filename
        self.data_dir = os.path.dirname(filename)
        self.fluid_types += self._get_fluid_types()
        super().__init__(filename, dataset_type, units_override=units_override)
        self.index_filename = filename
        self.storage_filename = storage_filename
        self.default_field = ("all", "T")
        self.sidecar_info = self._read_sidecar()

    def _get_sidecar_file(self):
        # returns the sidecar filename
        param_file = os.path.basename(self.parameter_filename)
        return os.path.join(
            self.data_dir, "." + os.path.splitext(param_file)[0] + ".json"
        )

    def _read_sidecar(self):
        # reads in sidecar file, returns empty dict if no file.

        # initialize for this .pvtu if it is not there
        sidecar = self._get_sidecar_file()
        if os.path.isfile(sidecar):
            with open(sidecar) as shandle:
                return json.load(shandle)[0]
        else:
            return {}

    def _write_sidecar(self):
        # writes current sidecar_info to json file (overwrite everything)
        with open(self._get_sidecar_file(), "w") as shandle:
            json.dump([self.sidecar_info], shandle)

    def _set_code_unit_attributes(self):
        setdefaultattr(self, "length_unit", self.quan(1.0, "m"))
        setdefaultattr(self, "mass_unit", self.quan(1.0, "kg"))
        setdefaultattr(self, "time_unit", self.quan(1.0, "s"))

    def _setup_coordinate_handler(self):
        # ensure correct ordering of axes so plots aren't rotated (z should always be
        # on the vertical axis).
        super()._setup_coordinate_handler()
        self.coordinates._x_pairs = (("x", "y"), ("y", "x"), ("z", "x"))
        self.coordinates._y_pairs = (("x", "z"), ("y", "z"), ("z", "y"))

    def _parse_parameter_file(self):

        # store the top level pvtu info
        with open(self.parameter_filename) as pvtu_fi:
            self.parameters["pXML"] = xmltodict.parse(pvtu_fi.read())

        pieces = self.parameters["pXML"]["VTKFile"]["PUnstructuredGrid"]["Piece"]
        if not isinstance(pieces, list):
            pieces = [pieces]
        self.parameters["vtu_files"] = [piece["@Source"] for piece in pieces]

        # periodicity: setting to False. Info is not in the pvtu files. Might
        # be able to check original.prm file if it exists.
        self._periodicity = (False, False, False)

        self.unique_identifier = self._get_unique_identifier()
        self.dimensionality = int(
            self.parameters["pXML"]["VTKFile"]["PUnstructuredGrid"]["PPoints"][
                "PDataArray"
            ]["@NumberOfComponents"]
        )
        self.domain_dimensions = np.ones(3, "int32")
        if self.dimensionality != 3:
            raise NotImplementedError(
                "The ASPECT frontend currently only "
                "supports 3D datasets, but this dataset "
                f"has {self.dimensionality} dimensions."
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
        # These attributes don't really make sense for ASPECT data, but yt warns
        # if they are not present, so we set them to dummy values here.
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
        fluid_types = tuple([f"connect{pcid}" for pcid in range(0, n_pieces)])
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

    def _add_piece_to_field_map(self, src_file, xmlPieces):
        """
        This constructs a look-up dictionary, field_to_piece_index:

        field_to_piece_index[vtufilename][piece_id][fieldname] --> piece index

        where piece index is the index of xmlPiece["PointData"]["DataArray"]
        for the desired field for a given vtu file and piece of the vtu file.
        Allows quick lookup of the field rather than looping through
        xmlPiece["PointData"]["DataArray"] until the desired field is found.
        Used in aspect.io._read_fluid_selection

        Parameters:
            src_file : str
                vtu filename
            xmlPieces : list
                list containing the pieces of the mesh for this vtu file

        stores field_to_piece_index in the self.parameters dictionary.
        """

        if "field_to_piece_index" not in self.parameters.keys():
            self.parameters["field_to_piece_index"] = {}

        if src_file not in self.parameters["field_to_piece_index"].keys():
            src_map = {}
            for piece_id in range(0, len(xmlPieces)):
                src_map[piece_id] = {}
                xmlPiece = xmlPieces[piece_id]
                for field_id, data_array in enumerate(
                    xmlPiece["PointData"]["DataArray"]
                ):
                    fname = data_array["@Name"]
                    src_map[piece_id][fname] = field_id

            self.parameters["field_to_piece_index"][src_file] = src_map

    def _read_pieces_from_single_vtu(self, src_file):
        with open(src_file) as data:
            xml = xmltodict.parse(data.read())
            xmlPieces = xml["VTKFile"]["UnstructuredGrid"]["Piece"]

        # initialize the containers for each piece of this vtu file. These get
        # combined into single arrays after reading each piece.
        conns = []  # connectivity (node --> index of coordinate arrays)
        alloffs = []  # element offsets
        x, y, z = [], [], []  # node coordinates
        n_p_cells = []  # nodes per cell of each piece
        pieceoff = 0  # connectivity index offset between pieces

        # do some piece sanitization and documentation
        if type(xmlPieces) != list:
            # handles the case where we have a single piece in this vtu file
            xmlPieces = [xmlPieces]
        self._add_piece_to_field_map(src_file, xmlPieces)

        # loop over pieces of this vtu, decode and add to containers
        for piece_id in range(0, len(xmlPieces)):
            coords, conn, offsets, cell_types = decode_piece(xmlPieces[piece_id])

            if len(np.unique(cell_types)) > 1:
                # check for multiple cell types within a single piece
                raise NotImplementedError(
                    f"multiple cell types in piece {piece_id} of vtu {src_file}. "
                    f"The ASPECT frontend can only handle single cell types at present."
                )

            n_nodes = conn.size
            n_cells = int(xmlPieces[piece_id]["@NumberOfCells"])
            n_p_cells.append(n_nodes // n_cells)

            x.extend(coords[:, 0])
            y.extend(coords[:, 1])
            z.extend(coords[:, 2])
            alloffs.extend(offsets)
            conns.extend(conn + pieceoff)  # connectivity across pieces

            # the connectivity array can repeat index references to coordinates
            # so update the offset number by the number of coords points
            # rather than the number of elements in conn:
            pieceoff += coords.shape[0]

        # concatenate some of the containers
        coords = np.array(np.column_stack([x, y, z]))
        alloffs = np.array(alloffs)
        conns = np.array(conns)

        # check that cell types are the same across the pieces. Not handling
        # multiple cell types right now as ASPECT generally does not mix
        # elements. May want to generalize this for wider vtu applicability...
        n_p_cell_vals = np.unique(n_p_cells)
        if len(n_p_cell_vals) > 1:
            raise NotImplementedError(
                f"{src_file} contains multiple cell types. The ASPECT frontend "
                f"can only handle a single cell type at present"
            )
        npc = n_p_cell_vals[0]

        conns = conns.reshape((conns.size // npc, npc))

        return coords, conns, alloffs, npc

    def _read_pieces(self):
        # reads coordinates and connectivity from pieces, returns the mesh
        # concatenated across all pieces.
        with open(self.parameter_filename) as pvtu_fi:
            pXML = xmltodict.parse(pvtu_fi.read())

        # the separate vtu files
        vtu_pieces = pXML["VTKFile"]["PUnstructuredGrid"]["Piece"]
        if not isinstance(vtu_pieces, list):
            vtu_pieces = [vtu_pieces]

        # containers that will be concatenated across vtu files
        conlist = []  # list of 2D connectivity arrays for each vtu
        coordlist = []  # global coordinates in each vtu
        node_offsets = []  # the offset within each vtu
        element_count = []  # total elements in each vtu
        nodes_per_cell = []  # nodes per cell for each vtu

        current_offset = 0  # the current NODE offset to global coordinate index
        for src in vtu_pieces:
            src_file = os.path.join(self.data_dir, src["@Source"])
            coord, con, all_offs, npc = self._read_pieces_from_single_vtu(src_file)
            nodes_per_cell.append(npc)

            node_offsets.append(current_offset)
            element_count.append(con.shape[0])
            con = con + current_offset  # offset to global
            conlist.append(con.astype("i8"))
            coordlist.append(coord.astype("f8"))

            # the connectivity array can repeat index references to coordinates
            # so update the offset number by the number of coords points
            # rather than the number of elements in conn:
            current_offset += coord.shape[0]

        # check that we have the same cell types across vtu files
        nodes_per_cell = np.unique(nodes_per_cell)
        if len(nodes_per_cell) > 1:
            raise NotImplementedError(
                "Found different cell types across vtu files, can only "
                "handle single cell types at present"
            )
        nodes_per_cell = nodes_per_cell[-1]
        if nodes_per_cell > 8:
            mylog.info(
                "Found higher order elements: most operations will only use first order components."
            )

        # concatenate across into a single mesh.
        coordlist = np.vstack(coordlist)
        conlist = np.vstack(conlist).astype("i8")
        return conlist, coordlist, node_offsets, element_count, nodes_per_cell

    def _load_domain_edge(self):
        """
        Loads the global boundaries for the domain edge

        """

        # check our sidecar file first:
        self.sidecar_info = self._read_sidecar()
        left_edge = self.sidecar_info.get("domain_left_edge", None)
        right_edge = self.sidecar_info.get("domain_right_edge", None)

        if not left_edge or not right_edge:
            _, coord, _, _, _ = self._read_pieces()
            left_edge = [coord[:, i].min() for i in range(self.dimensionality)]
            right_edge = [coord[:, i].max() for i in range(self.dimensionality)]
            self.sidecar_info["domain_left_edge"] = left_edge
            self.sidecar_info["domain_right_edge"] = right_edge
            self._write_sidecar()

        return np.array(left_edge), np.array(right_edge)

    @classmethod
    def _is_valid(self, datafile, *args, **kwargs):

        is_aspect = False
        if datafile.split(".")[-1] == "pvtu":
            try:
                # first we check if this was generated by deal.ii
                # pvtu files are small (~4K), so just parse the whole file and
                # check for the deal.II string
                with open(datafile) as data:
                    is_aspect = "generated by the deal.II library" in data.read()
            except OSError:
                return False

        # now we check if xmltodict is available
        if is_aspect:
            try:
                _ = hasattr(xmltodict, "parse")
                return True
            except ImportError:
                mylog.warning(
                    "It looks like you are trying to load an ASPECT pvtu file, "
                    "but xmltodict is not installed. Install with:"
                    "'pip install xmltodict' and then try loading your file again"
                )

        return False
