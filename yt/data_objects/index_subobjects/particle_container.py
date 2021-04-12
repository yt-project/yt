import contextlib

from more_itertools import always_iterable

from yt.data_objects.data_containers import YTFieldData
from yt.data_objects.selection_objects.data_selection_objects import (
    YTSelectionContainer,
)
from yt.utilities.exceptions import (
    YTDataSelectorNotImplemented,
    YTNonIndexedDataContainer,
)


def _non_indexed(name):
    def _func_non_indexed(self, *args, **kwargs):
        raise YTNonIndexedDataContainer(self)

    return _func_non_indexed


class ParticleContainer(YTSelectionContainer):
    _spatial = False
    _type_name = "particle_container"
    _skip_add = True
    _con_args = ("base_region", "base_selector", "data_files", "overlap_files")

    def __init__(
        self, base_region, base_selector, data_files, overlap_files=None, domain_id=-1
    ):
        if overlap_files is None:
            overlap_files = []
        self.field_data = YTFieldData()
        self.field_parameters = {}
        self.data_files = list(always_iterable(data_files))
        self.overlap_files = list(always_iterable(overlap_files))
        self.ds = self.data_files[0].ds
        self._last_mask = None
        self._last_selector_id = None
        self._current_particle_type = "all"
        # self._current_fluid_type = self.ds.default_fluid_type

        self.base_region = base_region
        self.base_selector = base_selector
        self._octree = None
        self._temp_spatial = False
        if isinstance(base_region, ParticleContainer):
            self._temp_spatial = base_region._temp_spatial
            self._octree = base_region._octree
        # To ensure there are not domains if global octree not used
        self.domain_id = -1

    def __reduce__(self):
        # we need to override the __reduce__ from data_containers as this method is not
        # a registered dataset method (i.e., ds.particle_container does not exist)
        arg_tuple = tuple(getattr(self, attr) for attr in self._con_args)
        return (self.__class__, arg_tuple)

    @property
    def selector(self):
        raise YTDataSelectorNotImplemented(self.oc_type_name)

    def select_particles(self, selector, x, y, z):
        mask = selector.select_points(x, y, z)
        return mask

    @contextlib.contextmanager
    def _expand_data_files(self):
        old_data_files = self.data_files
        old_overlap_files = self.overlap_files
        self.data_files = list(set(self.data_files + self.overlap_files))
        self.data_files.sort()
        self.overlap_files = []
        yield self
        self.data_files = old_data_files
        self.overlap_files = old_overlap_files

    def retrieve_ghost_zones(self, ngz, coarse_ghosts=False):
        gz_oct = self.octree.retrieve_ghost_zones(ngz, coarse_ghosts=coarse_ghosts)
        gz = ParticleContainer(
            gz_oct.base_region,
            gz_oct.base_selector,
            gz_oct.data_files,
            overlap_files=gz_oct.overlap_files,
            selector_mask=gz_oct.selector_mask,
            domain_id=gz_oct.domain_id,
        )
        gz._octree = gz_oct
        return gz

    select_blocks = _non_indexed("select_blocks")
    deposit = _non_indexed("deposit")
    smooth = _non_indexed("smooth")
    select_icoords = _non_indexed("select_icoords")
    select_fcoords = _non_indexed("select_fcoords")
    select_fwidth = _non_indexed("select_fwidth")
    select_ires = _non_indexed("select_ires")
    select = _non_indexed("select")
    count = _non_indexed("count")
    count_particles = _non_indexed("count_particles")
