import abc
import sys
import warnings
from functools import wraps
from types import ModuleType
from typing import Optional, Union

import numpy as np

from yt.config import ytcfg
from yt.data_objects.image_array import ImageArray
from yt.funcs import ensure_numpy_array, is_sequence, mylog
from yt.geometry.grid_geometry_handler import GridIndex
from yt.geometry.oct_geometry_handler import OctreeIndex
from yt.utilities.amr_kdtree.api import AMRKDTree
from yt.utilities.configure import YTConfig, configuration_callbacks
from yt.utilities.lib.bounding_volume_hierarchy import BVH
from yt.utilities.lib.misc_utilities import zlines, zpoints
from yt.utilities.lib.octree_raytracing import OctreeRayTracing
from yt.utilities.lib.partitioned_grid import PartitionedGrid
from yt.utilities.on_demand_imports import NotAModule
from yt.utilities.parallel_tools.parallel_analysis_interface import (
    ParallelAnalysisInterface,
)
from yt.visualization.image_writer import apply_colormap

from .transfer_function_helper import TransferFunctionHelper
from .transfer_functions import (
    ColorTransferFunction,
    ProjectionTransferFunction,
    TransferFunction,
)
from .utils import (
    data_source_or_all,
    get_corners,
    new_interpolated_projection_sampler,
    new_mesh_sampler,
    new_projection_sampler,
    new_volume_render_sampler,
)
from .zbuffer_array import ZBuffer

if sys.version_info >= (3, 8):
    from typing import Literal
else:
    from typing_extensions import Literal

OptionalModule = Union[ModuleType, NotAModule]
mesh_traversal: OptionalModule = NotAModule("pyembree")
mesh_construction: OptionalModule = NotAModule("pyembree")


def set_raytracing_engine(
    engine: Literal["yt", "embree"],
) -> None:
    """
    Safely switch raytracing engines at runtime.

    Parameters
    ----------

    engine: 'yt' or 'embree'
      - 'yt' selects the default engine.
      - 'embree' requires extra installation steps, see
        https://yt-project.org/doc/visualizing/unstructured_mesh_rendering.html?highlight=pyembree#optional-embree-installation

    Raises
    ------

    UserWarning
      Raised if the required engine is not available.
      In this case, the default engine is restored.

    """
    from yt.config import ytcfg

    global mesh_traversal, mesh_construction

    if engine == "embree":
        try:
            from yt.utilities.lib.embree_mesh import mesh_construction  # type: ignore
            from yt.utilities.lib.embree_mesh import mesh_traversal  # type: ignore
        except (ImportError, ValueError) as exc:
            # Catch ValueError in case size of objects in Cython change
            warnings.warn(
                "Failed to switch to embree raytracing engine. "
                f"The following error was raised:\n{exc}",
                stacklevel=2,
            )
            mesh_traversal = NotAModule("pyembree")
            mesh_construction = NotAModule("pyembree")
            ytcfg["yt", "ray_tracing_engine"] = "yt"
        else:
            ytcfg["yt", "ray_tracing_engine"] = "embree"
    else:
        mesh_traversal = NotAModule("pyembree")
        mesh_construction = NotAModule("pyembree")
        ytcfg["yt", "ray_tracing_engine"] = "yt"


def _init_raytracing_engine(ytcfg: YTConfig) -> None:
    # validate option from configuration file or fall back to default engine
    set_raytracing_engine(engine=ytcfg["yt", "ray_tracing_engine"])


configuration_callbacks.append(_init_raytracing_engine)


def invalidate_volume(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        ret = f(*args, **kwargs)
        obj = args[0]
        if isinstance(obj._transfer_function, ProjectionTransferFunction):
            obj.sampler_type = "projection"
            obj._log_field = False
            obj._use_ghost_zones = False
        del obj.volume
        obj._volume_valid = False
        return ret

    return wrapper


def validate_volume(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        obj = args[0]
        fields = [obj.field]
        log_fields = [obj.log_field]
        if obj.weight_field is not None:
            fields.append(obj.weight_field)
            log_fields.append(obj.log_field)
        if not obj._volume_valid:
            obj.volume.set_fields(
                fields, log_fields, no_ghost=(not obj.use_ghost_zones)
            )
        obj._volume_valid = True
        return f(*args, **kwargs)

    return wrapper


class RenderSource(ParallelAnalysisInterface, abc.ABC):
    """Base Class for Render Sources.

    Will be inherited for volumes, streamlines, etc.

    """

    volume_method: Optional[str] = None

    def __init__(self):
        super().__init__()
        self.opaque = False
        self.zbuffer = None

    @abc.abstractmethod
    def render(self, camera, zbuffer=None):
        pass

    @abc.abstractmethod
    def _validate(self):
        pass


class OpaqueSource(RenderSource):
    """A base class for opaque render sources.

    Will be inherited from for LineSources, BoxSources, etc.

    """

    def __init__(self):
        super().__init__()
        self.opaque = True

    def set_zbuffer(self, zbuffer):
        self.zbuffer = zbuffer


def create_volume_source(data_source, field):
    data_source = data_source_or_all(data_source)
    ds = data_source.ds
    index_class = ds.index.__class__
    if issubclass(index_class, GridIndex):
        return KDTreeVolumeSource(data_source, field)
    elif issubclass(index_class, OctreeIndex):
        return OctreeVolumeSource(data_source, field)
    else:
        raise NotImplementedError


class VolumeSource(RenderSource, abc.ABC):
    """A class for rendering data from a volumetric data source

    Examples of such sources include a sphere, cylinder, or the
    entire computational domain.

    A :class:`VolumeSource` provides the framework to decompose an arbitrary
    yt data source into bricks that can be traversed and volume rendered.

    Parameters
    ----------
    data_source: :class:`AMR3DData` or :class:`Dataset`, optional
        This is the source to be rendered, which can be any arbitrary yt
        data object or dataset.
    field : string
        The name of the field to be rendered.

    Examples
    --------

    The easiest way to make a VolumeSource is to use the volume_render
    function, so that the VolumeSource gets created automatically. This
    example shows how to do this and then access the resulting source:

    >>> import yt
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> im, sc = yt.volume_render(ds)
    >>> volume_source = sc.get_source(0)

    You can also create VolumeSource instances by hand and add them to Scenes.
    This example manually creates a VolumeSource, adds it to a scene, sets the
    camera, and renders an image.

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import (
    ...     Camera, Scene, create_volume_source)
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
    >>> sc = Scene()
    >>> source = create_volume_source(ds.all_data(), "density")
    >>> sc.add_source(source)
    >>> sc.add_camera()
    >>> im = sc.render()

    """

    _image = None
    data_source = None

    def __init__(self, data_source, field):
        r"""Initialize a new volumetric source for rendering."""
        super().__init__()
        self.data_source = data_source_or_all(data_source)
        field = self.data_source._determine_fields(field)[0]
        self.current_image = None
        self.check_nans = False
        self.num_threads = 0
        self.num_samples = 10
        self.sampler_type = "volume-render"

        self._volume_valid = False

        # these are caches for properties, defined below
        self._volume = None
        self._transfer_function = None
        self._field = field
        self._log_field = self.data_source.ds.field_info[field].take_log
        self._use_ghost_zones = False
        self._weight_field = None

        self.tfh = TransferFunctionHelper(self.data_source.pf)
        self.tfh.set_field(self.field)

    @property
    def transfer_function(self):
        """The transfer function associated with this VolumeSource"""
        if self._transfer_function is not None:
            return self._transfer_function

        if self.tfh.tf is not None:
            self._transfer_function = self.tfh.tf
            return self._transfer_function

        mylog.info("Creating transfer function")
        self.tfh.set_field(self.field)
        self.tfh.set_log(self.log_field)
        self.tfh.build_transfer_function()
        self.tfh.setup_default()
        self._transfer_function = self.tfh.tf

        return self._transfer_function

    @transfer_function.setter
    def transfer_function(self, value):
        self.tfh.tf = None
        valid_types = (
            TransferFunction,
            ColorTransferFunction,
            ProjectionTransferFunction,
            type(None),
        )
        if not isinstance(value, valid_types):
            raise RuntimeError(
                "transfer_function not a valid type, "
                "received object of type %s" % type(value)
            )
        if isinstance(value, ProjectionTransferFunction):
            self.sampler_type = "projection"
            if self._volume is not None:
                fields = [self.field]
                if self.weight_field is not None:
                    fields.append(self.weight_field)
                self._volume_valid = False
        self._transfer_function = value

    @property
    def volume(self):
        """The abstract volume associated with this VolumeSource

        This object does the heavy lifting to access data in an efficient manner
        using a KDTree
        """
        return self._get_volume()

    @volume.setter
    def volume(self, value):
        assert isinstance(value, AMRKDTree)
        del self._volume
        self._field = value.fields
        self._log_field = value.log_fields
        self._volume = value
        assert self._volume_valid

    @volume.deleter
    def volume(self):
        del self._volume
        self._volume = None

    @property
    def field(self):
        """The field to be rendered"""
        return self._field

    @field.setter  # type: ignore
    @invalidate_volume
    def field(self, value):
        field = self.data_source._determine_fields(value)
        if len(field) > 1:
            raise RuntimeError(
                "VolumeSource.field can only be a single field but received "
                "multiple fields: %s"
            ) % field
        field = field[0]
        if self._field != field:
            log_field = self.data_source.ds.field_info[field].take_log
            self.tfh.bounds = None
        else:
            log_field = self._log_field
        self._log_field = log_field
        self._field = value
        self.transfer_function = None
        self.tfh.set_field(value)
        self.tfh.set_log(log_field)

    @property
    def log_field(self):
        """Whether or not the field rendering is computed in log space"""
        return self._log_field

    @log_field.setter  # type: ignore
    @invalidate_volume
    def log_field(self, value):
        self.transfer_function = None
        self.tfh.set_log(value)
        self._log_field = value

    @property
    def use_ghost_zones(self):
        """Whether or not ghost zones are used to estimate vertex-centered data
        values at grid boundaries"""
        return self._use_ghost_zones

    @use_ghost_zones.setter  # type: ignore
    @invalidate_volume
    def use_ghost_zones(self, value):
        self._use_ghost_zones = value

    @property
    def weight_field(self):
        """The weight field for the rendering

        Currently this is only used for off-axis projections.
        """
        return self._weight_field

    @weight_field.setter  # type: ignore
    @invalidate_volume
    def weight_field(self, value):
        self._weight_field = value

    def set_transfer_function(self, transfer_function):
        """Set transfer function for this source"""
        self.transfer_function = transfer_function
        return self

    def _validate(self):
        """Make sure that all dependencies have been met"""
        if self.data_source is None:
            raise RuntimeError("Data source not initialized")

    def set_volume(self, volume):
        """Associates an AMRKDTree with the VolumeSource"""
        self.volume = volume
        return self

    def set_field(self, field):
        """Set the source's field to render

        Parameters
        ----------

        field: field name
            The field to render
        """
        self.field = field
        return self

    def set_log(self, log_field):
        """Set whether the rendering of the source's field is done in log space

        Generally volume renderings of data whose values span a large dynamic
        range should be done on log space and volume renderings of data with
        small dynamic range should be done in linear space.

        Parameters
        ----------

        log_field: boolean
            If True, the volume rendering will be done in log space, and if False
            will be done in linear space.
        """
        self.log_field = log_field
        return self

    def set_weight_field(self, weight_field):
        """Set the source's weight field

        .. note::

          This is currently only used for renderings using the
          ProjectionTransferFunction

        Parameters
        ----------

        weight_field: field name
            The weight field to use in the rendering
        """
        self.weight_field = weight_field
        return self

    def set_use_ghost_zones(self, use_ghost_zones):
        """Set whether or not interpolation at grid edges uses ghost zones

        Parameters
        ----------

        use_ghost_zones: boolean
            If True, the AMRKDTree estimates vertex centered data using ghost
            zones, which can eliminate seams in the resulting volume rendering.
            Defaults to False for performance reasons.

        """
        self.use_ghost_zones = use_ghost_zones
        return self

    def set_sampler(self, camera, interpolated=True):
        """Sets a volume render sampler

        The type of sampler is determined based on the ``sampler_type`` attribute
        of the VolumeSource. Currently the ``volume_render`` and ``projection``
        sampler types are supported.

        The 'interpolated' argument is only meaningful for projections. If True,
        the data is first interpolated to the cell vertices, and then
        tri-linearly interpolated to the ray sampling positions. If False, then
        the cell-centered data is simply accumulated along the
        ray. Interpolation is always performed for volume renderings.

        """
        if self.sampler_type == "volume-render":
            sampler = new_volume_render_sampler(camera, self)
        elif self.sampler_type == "projection" and interpolated:
            sampler = new_interpolated_projection_sampler(camera, self)
        elif self.sampler_type == "projection":
            sampler = new_projection_sampler(camera, self)
        else:
            NotImplementedError(f"{self.sampler_type} not implemented yet")
        self.sampler = sampler
        assert self.sampler is not None

    @abc.abstractmethod
    def _get_volume(self):
        """The abstract volume associated with this VolumeSource

        This object does the heavy lifting to access data in an efficient manner
        using a KDTree
        """
        pass

    @abc.abstractmethod
    @validate_volume
    def render(self, camera, zbuffer=None):
        """Renders an image using the provided camera

        Parameters
        ----------
        camera: :class:`yt.visualization.volume_rendering.camera.Camera` instance
            A volume rendering camera. Can be any type of camera.
        zbuffer: :class:`yt.visualization.volume_rendering.zbuffer_array.Zbuffer` instance  # noqa: E501
            A zbuffer array. This is used for opaque sources to determine the
            z position of the source relative to other sources. Only useful if
            you are manually calling render on multiple sources. Scene.render
            uses this internally.

        Returns
        -------
        A :class:`yt.data_objects.image_array.ImageArray` instance containing
        the rendered image.

        """
        pass

    def finalize_image(self, camera, image):
        """Parallel reduce the image.

        Parameters
        ----------
        camera: :class:`yt.visualization.volume_rendering.camera.Camera` instance
            The camera used to produce the volume rendering image.
        image: :class:`yt.data_objects.image_array.ImageArray` instance
            A reference to an image to fill
        """
        image.shape = camera.resolution[0], camera.resolution[1], 4
        # If the call is from VR, the image is rotated by 180 to get correct
        # up direction
        if not self.transfer_function.grey_opacity:
            image[:, :, 3] = 1
        return image

    def __repr__(self):
        disp = f"<Volume Source>:{str(self.data_source)} "
        disp += f"transfer_function:{str(self._transfer_function)}"
        return disp


class KDTreeVolumeSource(VolumeSource):
    volume_method = "KDTree"

    def _get_volume(self):
        """The abstract volume associated with this VolumeSource

        This object does the heavy lifting to access data in an efficient manner
        using a KDTree
        """

        if self._volume is None:
            mylog.info("Creating volume")
            volume = AMRKDTree(self.data_source.ds, data_source=self.data_source)
            self._volume = volume

        return self._volume

    @validate_volume
    def render(self, camera, zbuffer=None):
        """Renders an image using the provided camera

        Parameters
        ----------
        camera: :class:`yt.visualization.volume_rendering.camera.Camera`
            A volume rendering camera. Can be any type of camera.
        zbuffer: :class:`yt.visualization.volume_rendering.zbuffer_array.Zbuffer`
            A zbuffer array. This is used for opaque sources to determine the
            z position of the source relative to other sources. Only useful if
            you are manually calling render on multiple sources. Scene.render
            uses this internally.

        Returns
        -------
        A :class:`yt.data_objects.image_array.ImageArray` containing
        the rendered image.

        """
        self.zbuffer = zbuffer
        self.set_sampler(camera)
        assert self.sampler is not None

        mylog.debug("Casting rays")
        total_cells = 0
        if self.check_nans:
            for brick in self.volume.bricks:
                for data in brick.my_data:
                    if np.any(np.isnan(data)):
                        raise RuntimeError

        for brick in self.volume.traverse(camera.lens.viewpoint):
            mylog.debug("Using sampler %s", self.sampler)
            self.sampler(brick, num_threads=self.num_threads)
            total_cells += np.prod(brick.my_data[0].shape)
        mylog.debug("Done casting rays")
        self.current_image = self.finalize_image(camera, self.sampler.aimage)

        if zbuffer is None:
            self.zbuffer = ZBuffer(
                self.current_image, np.full(self.current_image.shape[:2], np.inf)
            )

        return self.current_image

    def finalize_image(self, camera, image):
        if self._volume is not None:
            image = self.volume.reduce_tree_images(image, camera.lens.viewpoint)

        return super().finalize_image(camera, image)


class OctreeVolumeSource(VolumeSource):
    volume_method = "Octree"

    def __init__(self, *args, **kwa):
        super().__init__(*args, **kwa)
        self.set_use_ghost_zones(True)

    def _get_volume(self):
        """The abstract volume associated with this VolumeSource

        This object does the heavy lifting to access data in an efficient manner
        using an octree.
        """

        if self._volume is None:
            mylog.info("Creating volume")
            volume = OctreeRayTracing(self.data_source)
            self._volume = volume

        return self._volume

    @validate_volume
    def render(self, camera, zbuffer=None):
        """Renders an image using the provided camera

        Parameters
        ----------
        camera: :class:`yt.visualization.volume_rendering.camera.Camera` instance
            A volume rendering camera. Can be any type of camera.
        zbuffer: :class:`yt.visualization.volume_rendering.zbuffer_array.Zbuffer` instance  # noqa: E501
            A zbuffer array. This is used for opaque sources to determine the
            z position of the source relative to other sources. Only useful if
            you are manually calling render on multiple sources. Scene.render
            uses this internally.

        Returns
        -------
        A :class:`yt.data_objects.image_array.ImageArray` instance containing
        the rendered image.

        """
        self.zbuffer = zbuffer
        self.set_sampler(camera)
        if self.sampler is None:
            raise RuntimeError(
                "No sampler set. This is likely a bug as it should never happen."
            )

        data = self.data_source

        dx = data["dx"].to("unitary").value[:, None]
        xyz = np.stack([data[_].to("unitary").value for _ in "x y z".split()], axis=-1)
        LE = xyz - dx / 2
        RE = xyz + dx / 2

        mylog.debug("Gathering data")
        dt = np.stack(list(self.volume.data) + [*LE.T, *RE.T], axis=-1).reshape(
            1, len(dx), 14, 1
        )
        mask = np.full(dt.shape[1:], 1, dtype=np.uint8)
        dims = np.array([1, 1, 1], dtype="int64")
        pg = PartitionedGrid(0, dt, mask, LE.flatten(), RE.flatten(), dims, n_fields=1)

        mylog.debug("Casting rays")
        self.sampler(pg, oct=self.volume.octree)
        mylog.debug("Done casting rays")

        self.current_image = self.finalize_image(camera, self.sampler.aimage)

        if zbuffer is None:
            self.zbuffer = ZBuffer(
                self.current_image, np.full(self.current_image.shape[:2], np.inf)
            )

        return self.current_image


class MeshSource(OpaqueSource):
    """A source for unstructured mesh data.

    This functionality requires the embree ray-tracing engine and the
    associated pyembree python bindings to be installed in order to
    function.

    A :class:`MeshSource` provides the framework to volume render
    unstructured mesh data.

    Parameters
    ----------
    data_source: :class:`AMR3DData` or :class:`Dataset`, optional
        This is the source to be rendered, which can be any arbitrary yt
        data object or dataset.
    field : string
        The name of the field to be rendered.

    Examples
    --------
    >>> source = MeshSource(ds, ("connect1", "convected"))
    """

    _image = None
    data_source = None

    def __init__(self, data_source, field):
        r"""Initialize a new unstructured mesh source for rendering."""
        super().__init__()
        self.data_source = data_source_or_all(data_source)
        field = self.data_source._determine_fields(field)[0]
        self.field = field
        self.volume = None
        self.current_image = None
        self.engine = ytcfg.get("yt", "ray_tracing_engine")

        # default color map
        self._cmap = ytcfg.get("yt", "default_colormap")
        self._color_bounds = None

        # default mesh annotation options
        self._annotate_mesh = False
        self._mesh_line_color = None
        self._mesh_line_alpha = 1.0

        # Error checking
        assert self.field is not None
        assert self.data_source is not None
        if self.field[0] == "all":
            raise NotImplementedError(
                "Mesh unions are not implemented for 3D rendering"
            )

        if self.engine == "embree":
            self.volume = mesh_traversal.YTEmbreeScene()
            self.build_volume_embree()
        elif self.engine == "yt":
            self.build_volume_bvh()
        else:
            raise NotImplementedError(
                "Invalid ray-tracing engine selected. Choices are 'embree' and 'yt'."
            )

    @property
    def cmap(self):
        """
        This is the name of the colormap that will be used when rendering
        this MeshSource object. Should be a string, like 'cmyt.arbre', or 'cmyt.dusk'.

        """
        return self._cmap

    @cmap.setter
    def cmap(self, cmap_name):
        self._cmap = cmap_name
        if hasattr(self, "data"):
            self.current_image = self.apply_colormap()

    @property
    def color_bounds(self):
        """
        These are the bounds that will be used with the colormap to the display
        the rendered image. Should be a (vmin, vmax) tuple, like (0.0, 2.0). If
        None, the bounds will be automatically inferred from the max and min of
        the rendered data.

        """
        return self._color_bounds

    @color_bounds.setter
    def color_bounds(self, bounds):
        self._color_bounds = bounds
        if hasattr(self, "data"):
            self.current_image = self.apply_colormap()

    def _validate(self):
        """Make sure that all dependencies have been met"""
        if self.data_source is None:
            raise RuntimeError("Data source not initialized.")

        if self.volume is None:
            raise RuntimeError("Volume not initialized.")

    def build_volume_embree(self):
        """

        This constructs the mesh that will be ray-traced by pyembree.

        """
        ftype, fname = self.field
        mesh_id = int(ftype[-1]) - 1
        index = self.data_source.ds.index
        offset = index.meshes[mesh_id]._index_offset
        field_data = self.data_source[self.field].d  # strip units

        vertices = index.meshes[mesh_id].connectivity_coords
        indices = index.meshes[mesh_id].connectivity_indices - offset

        # if this is an element field, promote to 2D here
        if len(field_data.shape) == 1:
            field_data = np.expand_dims(field_data, 1)

        # Here, we decide whether to render based on high-order or
        # low-order geometry. Right now, high-order geometry is only
        # implemented for 20-point hexes.
        if indices.shape[1] == 20 or indices.shape[1] == 10:
            self.mesh = mesh_construction.QuadraticElementMesh(
                self.volume, vertices, indices, field_data
            )
        else:
            # if this is another type of higher-order element, we demote
            # to 1st order here, for now.
            if indices.shape[1] == 27:
                # hexahedral
                mylog.warning("27-node hexes not yet supported, dropping to 1st order.")
                field_data = field_data[:, 0:8]
                indices = indices[:, 0:8]

            self.mesh = mesh_construction.LinearElementMesh(
                self.volume, vertices, indices, field_data
            )

    def build_volume_bvh(self):
        """

        This constructs the mesh that will be ray-traced.

        """
        ftype, fname = self.field
        mesh_id = int(ftype[-1]) - 1
        index = self.data_source.ds.index
        offset = index.meshes[mesh_id]._index_offset
        field_data = self.data_source[self.field].d  # strip units

        vertices = index.meshes[mesh_id].connectivity_coords
        indices = index.meshes[mesh_id].connectivity_indices - offset

        # if this is an element field, promote to 2D here
        if len(field_data.shape) == 1:
            field_data = np.expand_dims(field_data, 1)

        # Here, we decide whether to render based on high-order or
        # low-order geometry.
        if indices.shape[1] == 27:
            # hexahedral
            mylog.warning("27-node hexes not yet supported, dropping to 1st order.")
            field_data = field_data[:, 0:8]
            indices = indices[:, 0:8]

        self.volume = BVH(vertices, indices, field_data)

    def render(self, camera, zbuffer=None):
        """Renders an image using the provided camera

        Parameters
        ----------
        camera: :class:`yt.visualization.volume_rendering.camera.Camera`
            A volume rendering camera. Can be any type of camera.
        zbuffer: :class:`yt.visualization.volume_rendering.zbuffer_array.Zbuffer`
            A zbuffer array. This is used for opaque sources to determine the
            z position of the source relative to other sources. Only useful if
            you are manually calling render on multiple sources. Scene.render
            uses this internally.

        Returns
        -------
        A :class:`yt.data_objects.image_array.ImageArray` containing
        the rendered image.

        """

        shape = (camera.resolution[0], camera.resolution[1], 4)
        if zbuffer is None:
            empty = np.empty(shape, dtype="float64")
            z = np.empty(empty.shape[:2], dtype="float64")
            empty[:] = 0.0
            z[:] = np.inf
            zbuffer = ZBuffer(empty, z)
        elif zbuffer.rgba.shape != shape:
            zbuffer = ZBuffer(zbuffer.rgba.reshape(shape), zbuffer.z.reshape(shape[:2]))
        self.zbuffer = zbuffer

        self.sampler = new_mesh_sampler(camera, self, engine=self.engine)

        mylog.debug("Casting rays")
        self.sampler(self.volume)
        mylog.debug("Done casting rays")

        self.finalize_image(camera)
        self.current_image = self.apply_colormap()

        zbuffer += ZBuffer(self.current_image.astype("float64"), self.sampler.azbuffer)
        zbuffer.rgba = ImageArray(zbuffer.rgba)
        self.zbuffer = zbuffer
        self.current_image = self.zbuffer.rgba

        if self._annotate_mesh:
            self.current_image = self.annotate_mesh_lines(
                self._mesh_line_color, self._mesh_line_alpha
            )

        return self.current_image

    def finalize_image(self, camera):
        sam = self.sampler

        # reshape data
        Nx = camera.resolution[0]
        Ny = camera.resolution[1]
        self.data = sam.aimage[:, :, 0].reshape(Nx, Ny)

    def annotate_mesh_lines(self, color=None, alpha=1.0):
        r"""

        Modifies this MeshSource by drawing the mesh lines.
        This modifies the current image by drawing the element
        boundaries and returns the modified image.

        Parameters
        ----------
        color: array_like of shape (4,), optional
            The RGBA value to use to draw the mesh lines.
            Default is black.
        alpha : float, optional
            The opacity of the mesh lines. Default is 255 (solid).

        """

        self.annotate_mesh = True
        self._mesh_line_color = color
        self._mesh_line_alpha = alpha

        if color is None:
            color = np.array([0, 0, 0, alpha])

        locs = (self.sampler.amesh_lines == 1,)

        self.current_image[:, :, 0][locs] = color[0]
        self.current_image[:, :, 1][locs] = color[1]
        self.current_image[:, :, 2][locs] = color[2]
        self.current_image[:, :, 3][locs] = color[3]

        return self.current_image

    def apply_colormap(self):
        """

        Applies a colormap to the current image without re-rendering.

        Returns
        -------
        current_image : A new image with the specified color scale applied to
            the underlying data.


        """

        image = (
            apply_colormap(
                self.data, color_bounds=self._color_bounds, cmap_name=self._cmap
            )
            / 255.0
        )
        alpha = image[:, :, 3]
        alpha[self.sampler.aimage_used == -1] = 0.0
        image[:, :, 3] = alpha
        return image

    def __repr__(self):
        disp = f"<Mesh Source>:{str(self.data_source)} "
        return disp


class PointSource(OpaqueSource):
    r"""A rendering source of opaque points in the scene.

    This class provides a mechanism for adding points to a scene; these
    points will be opaque, and can also be colored.

    Parameters
    ----------
    positions: array_like of shape (N, 3)
        The positions of points to be added to the scene. If specified with no
        units, the positions will be assumed to be in code units.
    colors : array_like of shape (N, 4), optional
        The colors of the points, including an alpha channel, in floating
        point running from 0..1.
    color_stride : int, optional
        The stride with which to access the colors when putting them on the
        scene.
    radii : array_like of shape (N), optional
        The radii of the points in the final image, in pixels (int)

    Examples
    --------

    This example creates a volume rendering and adds 1000 random points to
    the image:

    >>> import yt
    >>> import numpy as np
    >>> from yt.visualization.volume_rendering.api import PointSource
    >>> from yt.units import kpc
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    >>> im, sc = yt.volume_render(ds)

    >>> npoints = 1000
    >>> vertices = np.random.random([npoints, 3]) * 1000 * kpc
    >>> colors = np.random.random([npoints, 4])
    >>> colors[:, 3] = 1.0

    >>> points = PointSource(vertices, colors=colors)
    >>> sc.add_source(points)

    >>> im = sc.render()

    """

    _image = None
    data_source = None

    def __init__(self, positions, colors=None, color_stride=1, radii=None):
        assert positions.ndim == 2 and positions.shape[1] == 3
        if colors is not None:
            assert colors.ndim == 2 and colors.shape[1] == 4
            assert colors.shape[0] == positions.shape[0]
        if not is_sequence(radii):
            if radii is not None:  # broadcast the value
                radii = radii * np.ones(positions.shape[0], dtype="int64")
            else:  # default radii to 0 pixels (i.e. point is 1 pixel wide)
                radii = np.zeros(positions.shape[0], dtype="int64")
        else:
            assert radii.ndim == 1
            assert radii.shape[0] == positions.shape[0]
        self.positions = positions
        # If colors aren't individually set, make black with full opacity
        if colors is None:
            colors = np.ones((len(positions), 4))
        self.colors = colors
        self.color_stride = color_stride
        self.radii = radii

    def _validate(self):
        pass

    def render(self, camera, zbuffer=None):
        """Renders an image using the provided camera

        Parameters
        ----------
        camera: :class:`yt.visualization.volume_rendering.camera.Camera`
            A volume rendering camera. Can be any type of camera.
        zbuffer: :class:`yt.visualization.volume_rendering.zbuffer_array.Zbuffer`
            A zbuffer array. This is used for opaque sources to determine the
            z position of the source relative to other sources. Only useful if
            you are manually calling render on multiple sources. Scene.render
            uses this internally.

        Returns
        -------
        A :class:`yt.data_objects.image_array.ImageArray` containing
        the rendered image.

        """
        vertices = self.positions
        if zbuffer is None:
            empty = camera.lens.new_image(camera)
            z = np.empty(empty.shape[:2], dtype="float64")
            empty[:] = 0.0
            z[:] = np.inf
            zbuffer = ZBuffer(empty, z)
        else:
            empty = zbuffer.rgba
            z = zbuffer.z

        # DRAW SOME POINTS
        camera.lens.setup_box_properties(camera)
        px, py, dz = camera.lens.project_to_plane(camera, vertices)

        zpoints(empty, z, px, py, dz, self.colors, self.radii, self.color_stride)

        self.zbuffer = zbuffer
        return zbuffer

    def __repr__(self):
        disp = "<Point Source>"
        return disp


class LineSource(OpaqueSource):
    r"""A render source for a sequence of opaque line segments.

    This class provides a mechanism for adding lines to a scene; these
    points will be opaque, and can also be colored.

    .. note::

        If adding a LineSource to your rendering causes the image to appear
        blank or fades a VolumeSource, try lowering the values specified in
        the alpha channel of the ``colors`` array.

    Parameters
    ----------
    positions: array_like of shape (N, 2, 3)
        The positions of the starting and stopping points for each line.
        For example,positions[0][0] and positions[0][1] would give the (x, y, z)
        coordinates of the beginning and end points of the first line,
        respectively. If specified with no units, assumed to be in code units.
    colors : array_like of shape (N, 4), optional
        The colors of the points, including an alpha channel, in floating
        point running from 0..1.  The four channels correspond to r, g, b, and
        alpha values. Note that they correspond to the line segment succeeding
        each point; this means that strictly speaking they need only be (N-1)
        in length.
    color_stride : int, optional
        The stride with which to access the colors when putting them on the
        scene.

    Examples
    --------

    This example creates a volume rendering and then adds some random lines
    to the image:

    >>> import yt
    >>> import numpy as np
    >>> from yt.visualization.volume_rendering.api import LineSource
    >>> from yt.units import kpc
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    >>> im, sc = yt.volume_render(ds)

    >>> nlines = 4
    >>> vertices = np.random.random([nlines, 2, 3]) * 600 * kpc
    >>> colors = np.random.random([nlines, 4])
    >>> colors[:, 3] = 1.0

    >>> lines = LineSource(vertices, colors)
    >>> sc.add_source(lines)

    >>> im = sc.render()

    """

    _image = None
    data_source = None

    def __init__(self, positions, colors=None, color_stride=1):
        super().__init__()

        assert positions.ndim == 3
        assert positions.shape[1] == 2
        assert positions.shape[2] == 3
        if colors is not None:
            assert colors.ndim == 2
            assert colors.shape[1] == 4

        # convert the positions to the shape expected by zlines, below
        N = positions.shape[0]
        self.positions = positions.reshape((2 * N, 3))

        # If colors aren't individually set, make black with full opacity
        if colors is None:
            colors = np.ones((len(positions), 4))
        self.colors = colors
        self.color_stride = color_stride

    def _validate(self):
        pass

    def render(self, camera, zbuffer=None):
        """Renders an image using the provided camera

        Parameters
        ----------
        camera: :class:`yt.visualization.volume_rendering.camera.Camera`
            A volume rendering camera. Can be any type of camera.
        zbuffer: :class:`yt.visualization.volume_rendering.zbuffer_array.Zbuffer`
            z position of the source relative to other sources. Only useful if
            you are manually calling render on multiple sources. Scene.render
            uses this internally.

        Returns
        -------
        A :class:`yt.data_objects.image_array.ImageArray` containing
        the rendered image.

        """
        vertices = self.positions
        if zbuffer is None:
            empty = camera.lens.new_image(camera)
            z = np.empty(empty.shape[:2], dtype="float64")
            empty[:] = 0.0
            z[:] = np.inf
            zbuffer = ZBuffer(empty, z)
        else:
            empty = zbuffer.rgba
            z = zbuffer.z

        # DRAW SOME LINES
        camera.lens.setup_box_properties(camera)
        px, py, dz = camera.lens.project_to_plane(camera, vertices)

        px = px.astype("int64")
        py = py.astype("int64")

        if len(px.shape) == 1:
            zlines(
                empty, z, px, py, dz, self.colors.astype("float64"), self.color_stride
            )
        else:
            # For stereo-lens, two sets of pos for each eye are contained
            # in px...pz
            zlines(
                empty,
                z,
                px[0, :],
                py[0, :],
                dz[0, :],
                self.colors.astype("float64"),
                self.color_stride,
            )
            zlines(
                empty,
                z,
                px[1, :],
                py[1, :],
                dz[1, :],
                self.colors.astype("float64"),
                self.color_stride,
            )

        self.zbuffer = zbuffer
        return zbuffer

    def __repr__(self):
        disp = "<Line Source>"
        return disp


class BoxSource(LineSource):
    r"""A render source for a box drawn with line segments.
    This render source will draw a box, with transparent faces, in data
    space coordinates.  This is useful for annotations.

    Parameters
    ----------
    left_edge: array-like of shape (3,), float
        The left edge coordinates of the box.
    right_edge : array-like of shape (3,), float
        The right edge coordinates of the box.
    color : array-like of shape (4,), float, optional
        The colors (including alpha) to use for the lines.
        Default is black with an alpha of 1.0.

    Examples
    --------

    This example shows how to use BoxSource to add an outline of the
    domain boundaries to a volume rendering.

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import BoxSource
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    >>> im, sc = yt.volume_render(ds)

    >>> box_source = BoxSource(
    ...     ds.domain_left_edge, ds.domain_right_edge, [1.0, 1.0, 1.0, 1.0]
    ... )
    >>> sc.add_source(box_source)

    >>> im = sc.render()

    """

    def __init__(self, left_edge, right_edge, color=None):

        assert left_edge.shape == (3,)
        assert right_edge.shape == (3,)

        if color is None:
            color = np.array([1.0, 1.0, 1.0, 1.0])

        color = ensure_numpy_array(color)
        color.shape = (1, 4)
        corners = get_corners(left_edge.copy(), right_edge.copy())
        order = [0, 1, 1, 2, 2, 3, 3, 0]
        order += [4, 5, 5, 6, 6, 7, 7, 4]
        order += [0, 4, 1, 5, 2, 6, 3, 7]
        vertices = np.empty([24, 3])
        for i in range(3):
            vertices[:, i] = corners[order, i, ...].ravel(order="F")
        vertices = vertices.reshape((12, 2, 3))

        super().__init__(vertices, color, color_stride=24)

    def _validate(self):
        pass


class GridSource(LineSource):
    r"""A render source for drawing grids in a scene.

    This render source will draw blocks that are within a given data
    source, by default coloring them by their level of resolution.

    Parameters
    ----------
    data_source: :class:`~yt.data_objects.api.DataContainer`
        The data container that will be used to identify grids to draw.
    alpha : float
        The opacity of the grids to draw.
    cmap : color map name
        The color map to use to map resolution levels to color.
    min_level : int, optional
        Minimum level to draw
    max_level : int, optional
        Maximum level to draw

    Examples
    --------

    This example makes a volume rendering and adds outlines of all the
    AMR grids in the simulation:

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import GridSource
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    >>> im, sc = yt.volume_render(ds)

    >>> grid_source = GridSource(ds.all_data(), alpha=1.0)

    >>> sc.add_source(grid_source)

    >>> im = sc.render()

    This example does the same thing, except it only draws the grids
    that are inside a sphere of radius (0.1, "unitary") located at the
    domain center:

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import GridSource
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    >>> im, sc = yt.volume_render(ds)

    >>> dd = ds.sphere("c", (0.1, "unitary"))
    >>> grid_source = GridSource(dd, alpha=1.0)

    >>> sc.add_source(grid_source)

    >>> im = sc.render()

    """

    def __init__(
        self, data_source, alpha=0.3, cmap=None, min_level=None, max_level=None
    ):
        self.data_source = data_source_or_all(data_source)
        corners = []
        levels = []
        for block, _mask in self.data_source.blocks:
            block_corners = np.array(
                [
                    [block.LeftEdge[0], block.LeftEdge[1], block.LeftEdge[2]],
                    [block.RightEdge[0], block.LeftEdge[1], block.LeftEdge[2]],
                    [block.RightEdge[0], block.RightEdge[1], block.LeftEdge[2]],
                    [block.LeftEdge[0], block.RightEdge[1], block.LeftEdge[2]],
                    [block.LeftEdge[0], block.LeftEdge[1], block.RightEdge[2]],
                    [block.RightEdge[0], block.LeftEdge[1], block.RightEdge[2]],
                    [block.RightEdge[0], block.RightEdge[1], block.RightEdge[2]],
                    [block.LeftEdge[0], block.RightEdge[1], block.RightEdge[2]],
                ],
                dtype="float64",
            )
            corners.append(block_corners)
            levels.append(block.Level)
        corners = np.dstack(corners)
        levels = np.array(levels)
        if cmap is None:
            cmap = ytcfg.get("yt", "default_colormap")

        if max_level is not None:
            subset = levels <= max_level
            levels = levels[subset]
            corners = corners[:, :, subset]
        if min_level is not None:
            subset = levels >= min_level
            levels = levels[subset]
            corners = corners[:, :, subset]

        colors = (
            apply_colormap(
                levels * 1.0,
                color_bounds=[0, self.data_source.ds.index.max_level],
                cmap_name=cmap,
            )[0, :, :]
            / 255.0
        )
        colors[:, 3] = alpha

        order = [0, 1, 1, 2, 2, 3, 3, 0]
        order += [4, 5, 5, 6, 6, 7, 7, 4]
        order += [0, 4, 1, 5, 2, 6, 3, 7]

        vertices = np.empty([corners.shape[2] * 2 * 12, 3])
        for i in range(3):
            vertices[:, i] = corners[order, i, ...].ravel(order="F")
        vertices = vertices.reshape((corners.shape[2] * 12, 2, 3))

        super().__init__(vertices, colors, color_stride=24)


class CoordinateVectorSource(OpaqueSource):
    r"""Draw coordinate vectors on the scene.

    This will draw a set of coordinate vectors on the camera image.  They
    will appear in the lower right of the image.

    Parameters
    ----------
    colors: array-like of shape (3,4), optional
        The RGBA values to use to draw the x, y, and z vectors. The default is
        [[1, 0, 0, alpha], [0, 1, 0, alpha], [0, 0, 1, alpha]]  where ``alpha``
        is set by the parameter below. If ``colors`` is set then ``alpha`` is
        ignored.
    alpha : float, optional
        The opacity of the vectors.

    Examples
    --------

    >>> import yt
    >>> from yt.visualization.volume_rendering.api import \
    ...     CoordinateVectorSource
    >>> ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")

    >>> im, sc = yt.volume_render(ds)

    >>> coord_source = CoordinateVectorSource()

    >>> sc.add_source(coord_source)

    >>> im = sc.render()

    """

    def __init__(self, colors=None, alpha=1.0):
        super().__init__()
        # If colors aren't individually set, make black with full opacity
        if colors is None:
            colors = np.zeros((3, 4))
            colors[0, 0] = 1.0  # x is red
            colors[1, 1] = 1.0  # y is green
            colors[2, 2] = 1.0  # z is blue
            colors[:, 3] = alpha
        self.colors = colors

    def _validate(self):
        pass

    def render(self, camera, zbuffer=None):
        """Renders an image using the provided camera

        Parameters
        ----------
        camera: :class:`yt.visualization.volume_rendering.camera.Camera`
            A volume rendering camera. Can be any type of camera.
        zbuffer: :class:`yt.visualization.volume_rendering.zbuffer_array.Zbuffer`
            A zbuffer array. This is used for opaque sources to determine the
            z position of the source relative to other sources. Only useful if
            you are manually calling render on multiple sources. Scene.render
            uses this internally.

        Returns
        -------
        A :class:`yt.data_objects.image_array.ImageArray` containing
        the rendered image.

        """
        camera.lens.setup_box_properties(camera)
        center = camera.focus
        # Get positions at the focus
        positions = np.zeros([6, 3])
        positions[:] = center

        # Create vectors in the x,y,z directions
        for i in range(3):
            positions[2 * i + 1, i] += camera.width.in_units("code_length").d[i] / 16.0

        # Project to the image plane
        px, py, dz = camera.lens.project_to_plane(camera, positions)

        if len(px.shape) == 1:
            dpx = px[1::2] - px[::2]
            dpy = py[1::2] - py[::2]

            # Set the center of the coordinates to be in the lower left of the image
            lpx = camera.resolution[0] / 8
            lpy = camera.resolution[1] - camera.resolution[1] / 8  # Upside-downsies

            # Offset the pixels according to the projections above
            px[::2] = lpx
            px[1::2] = lpx + dpx
            py[::2] = lpy
            py[1::2] = lpy + dpy
            dz[:] = 0.0
        else:
            # For stereo-lens, two sets of pos for each eye are contained in px...pz
            dpx = px[:, 1::2] - px[:, ::2]
            dpy = py[:, 1::2] - py[:, ::2]

            lpx = camera.resolution[0] / 16
            lpy = camera.resolution[1] - camera.resolution[1] / 8  # Upside-downsies

            # Offset the pixels according to the projections above
            px[:, ::2] = lpx
            px[:, 1::2] = lpx + dpx
            px[1, :] += camera.resolution[0] / 2
            py[:, ::2] = lpy
            py[:, 1::2] = lpy + dpy
            dz[:, :] = 0.0

        # Create a zbuffer if needed
        if zbuffer is None:
            empty = camera.lens.new_image(camera)
            z = np.empty(empty.shape[:2], dtype="float64")
            empty[:] = 0.0
            z[:] = np.inf
            zbuffer = ZBuffer(empty, z)
        else:
            empty = zbuffer.rgba
            z = zbuffer.z

        # Draw the vectors

        px = px.astype("int64")
        py = py.astype("int64")

        if len(px.shape) == 1:
            zlines(empty, z, px, py, dz, self.colors.astype("float64"))
        else:
            # For stereo-lens, two sets of pos for each eye are contained
            # in px...pz
            zlines(
                empty, z, px[0, :], py[0, :], dz[0, :], self.colors.astype("float64")
            )
            zlines(
                empty, z, px[1, :], py[1, :], dz[1, :], self.colors.astype("float64")
            )

        # Set the new zbuffer
        self.zbuffer = zbuffer
        return zbuffer

    def __repr__(self):
        disp = "<Coordinates Source>"
        return disp
