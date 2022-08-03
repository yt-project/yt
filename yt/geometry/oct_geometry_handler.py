from typing import Optional, Tuple

import numpy as np

from yt.fields.field_detector import FieldDetector
from yt.geometry.geometry_handler import Index
from yt.utilities.logger import ytLogger as mylog


class OctreeIndex(Index):
    """The Index subclass for oct AMR datasets"""

    def _setup_geometry(self):
        mylog.debug("Initializing Octree Geometry Handler.")
        self._initialize_oct_handler()

    def get_smallest_dx(self):
        """
        Returns (in code units) the smallest cell size in the simulation.
        """
        return (
            self.dataset.domain_width
            / (self.dataset.domain_dimensions * 2 ** (self.max_level))
        ).min()

    def convert(self, unit):
        return self.dataset.conversion_factors[unit]

    def _add_mesh_sampling_particle_field(self, deposit_field, ftype, ptype):
        units = self.ds.field_info[ftype, deposit_field].units
        take_log = self.ds.field_info[ftype, deposit_field].take_log
        field_name = f"cell_{ftype}_{deposit_field}"

        def _cell_index(field, data):
            # Get the position of the particles
            pos = data[ptype, "particle_position"]
            Npart = pos.shape[0]
            ret = np.zeros(Npart)
            tmp = np.zeros(Npart)

            if isinstance(data, FieldDetector):
                return ret

            remaining = np.ones(Npart, dtype=bool)
            Nremaining = Npart

            Nobjs = len(data._current_chunk.objs)
            Nbits = int(np.ceil(np.log2(Nobjs)))

            # Sort objs by decreasing number of octs
            enumerated_objs = sorted(
                enumerate(data._current_chunk.objs),
                key=lambda arg: arg[1].oct_handler.nocts,
                reverse=True,
            )
            for i, obj in enumerated_objs:
                if Nremaining == 0:
                    break
                icell = (
                    obj[("index", "ones")].T.reshape(-1).astype(np.int64).cumsum().value
                    - 1
                )
                mesh_data = ((icell << Nbits) + i).astype(np.float64)
                # Access the mesh data and attach them to their particles
                tmp[:Nremaining] = obj.mesh_sampling_particle_field(
                    pos[remaining], mesh_data
                )

                ret[remaining] = tmp[:Nremaining]

                remaining[remaining] = np.isnan(tmp[:Nremaining])
                Nremaining = remaining.sum()

            return data.ds.arr(ret.astype(np.float64), units="1")

        def _mesh_sampling_particle_field(field, data):
            """
            Create a grid field for particle quantities using given method.
            """
            ones = data[ptype, "particle_ones"]

            # Access "cell_index" field
            Npart = ones.shape[0]
            ret = np.zeros(Npart)
            cell_index = np.array(data[ptype, "cell_index"], np.int64)

            if isinstance(data, FieldDetector):
                return ret

            # The index of the obj is stored on the first bits
            Nobjs = len(data._current_chunk.objs)
            Nbits = int(np.ceil(np.log2(Nobjs)))
            icell = cell_index >> Nbits
            iobj = cell_index - (icell << Nbits)
            for i, subset in enumerate(data._current_chunk.objs):
                mask = iobj == i

                subset.field_parameters = data.field_parameters

                cell_data = subset[ftype, deposit_field].T.reshape(-1)

                ret[mask] = cell_data[icell[mask]]

            return data.ds.arr(ret, units=cell_data.units)

        if (ptype, "cell_index") not in self.ds.derived_field_list:
            self.ds.add_field(
                (ptype, "cell_index"),
                function=_cell_index,
                sampling_type="particle",
                units="1",
            )

        self.ds.add_field(
            (ptype, field_name),
            function=_mesh_sampling_particle_field,
            sampling_type="particle",
            units=units,
            take_log=take_log,
        )

    def _icoords_to_fcoords(
        self,
        icoords: np.ndarray,
        ires: np.ndarray,
        axes: Optional[Tuple[int, ...]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Accepts icoords and ires and returns appropriate fcoords and fwidth.
        Mostly useful for cases where we have irregularly spaced or structured
        grids.
        """
        dds = self.ds.domain_width[(axes,)] / (
            self.ds.domain_dimensions[
                axes,
            ]
            * self.ds.refine_by ** ires[:, None]
        )
        pos = (0.5 + icoords) * dds + self.ds.domain_left_edge[
            axes,
        ]
        return pos, dds
