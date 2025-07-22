import abc

import numpy as np

import yt.utilities.lib.interpolators as lib


class _LinearInterpolator(abc.ABC):
    _ndim: int

    def __init__(self, table, field_names, truncate=False, *, store_table=True):
        if store_table:
            self.table = table.astype("float64", copy=False)
        else:
            self.table = None
        self.table_shape = table.shape
        self.truncate = truncate
        self._field_names = field_names

    def _validate_table(self, table_override):
        if table_override is None:
            if self.table is None:
                msg = (
                    f"You must either store the table used when initializing a new "
                    f"{type(self).__name__} (set `store_table=True`) or you provide a `table` when "
                    f"calling {type(self).__name__}"
                )
                raise ValueError(msg)
            return self.table

        if table_override.shape != self.table_shape:
            msg = f"The table_override shape, {table_override.shape}, must match the base table shape, {self.table_shape}"
            raise ValueError(msg)

        return table_override.astype("float64", copy=False)

    def _get_digitized_arrays(self, data_object):
        return_arrays = []
        for dim in "xyzw"[: self._ndim]:
            dim_name = getattr(self, f"{dim}_name")
            dim_bins = getattr(self, f"{dim}_bins")

            dim_vals = data_object[dim_name].astype("float64").ravel()
            dim_i = (np.digitize(dim_vals, dim_bins) - 1).astype("int64")
            if np.any((dim_i == -1) | (dim_i == len(dim_bins) - 1)):
                if not self.truncate:
                    msg = (
                        f"The dimension values for {dim_name} and data object {data_object} are outside the bounds "
                        f"of the table! You can avoid this error by providing truncate=True to the interpolator or "
                        f"by adjusting the data object to remain inside the table bounds. But for now, dunno what "
                        f"to do, so dying."
                    )
                    raise ValueError(msg)
                else:
                    dim_i = np.minimum(np.maximum(dim_i, 0), len(dim_bins) - 2)
            return_arrays.append(dim_vals)
            return_arrays.append(dim_i)

        return return_arrays

    def _validate_bin_boundaries(self, boundaries):
        # boundaries: tuple of ndarrays
        for idim in range(self._ndim):
            if boundaries[idim].size != self.table_shape[idim]:
                msg = f"{self._field_names[idim]} bins array not the same length as the data."
                raise ValueError(msg)


class UnilinearFieldInterpolator(_LinearInterpolator):
    _ndim = 1

    def __init__(
        self, table, boundaries, field_names, truncate=False, *, store_table=True
    ):
        r"""Initialize a 1D interpolator for field data.

        table : array
            The data table over which interpolation is performed.
        boundaries: tuple or array
            If a tuple, this should specify the upper and lower bounds
            for the bins of the data table.  This assumes the bins are
            evenly spaced.  If an array, this specifies the bins
            explicitly.
        field_names: str
            Name of the field to be used as input data for interpolation.
        truncate: bool
            If False (default), an exception is raised if the input values are
            outside the bounds of the table.  If True, extrapolation is
            performed.
        store_table: bool
            If True (default), the full table is stored in the interpolator.
            If False, only the shape of the input table is stored and a full
            table must be provided when calling the interpolator.

        Examples
        --------

        >>> ad = ds.all_data()
        >>> table_data = np.random.random(64)
        >>> interp = UnilinearFieldInterpolator(table_data, (0.0, 1.0), "x",
                                                truncate=True)
        >>> field_data = interp(ad)

        If you want to re-use the interpolator with table_data of the same shape
        but different values, you can also supply the `table` keyword argument when
        calling the interpolator:

        >>> new_table_data = np.random.random(64)
        >>> field_data = interp(ad, table=new_table_data)

        """
        super().__init__(table, field_names, truncate=truncate, store_table=store_table)
        self.x_name = field_names
        if isinstance(boundaries, np.ndarray):
            self._validate_bin_boundaries((boundaries,))
            self.x_bins = boundaries
        else:
            x0, x1 = boundaries
            self.x_bins = np.linspace(x0, x1, table.shape[0], dtype="float64")

    def __call__(self, data_object, *, table=None):
        table = self._validate_table(table)
        orig_shape = data_object[self.x_name].shape
        x_vals, x_i = self._get_digitized_arrays(data_object)
        my_vals = np.zeros(x_vals.shape, dtype="float64")
        lib.UnilinearlyInterpolate(table, x_vals, self.x_bins, x_i, my_vals)
        my_vals.shape = orig_shape
        return my_vals


class BilinearFieldInterpolator(_LinearInterpolator):
    _ndim = 2

    def __init__(
        self, table, boundaries, field_names, truncate=False, *, store_table=True
    ):
        r"""Initialize a 2D interpolator for field data.

        table : array
            The data table over which interpolation is performed.
        boundaries: tuple
            Either a tuple of lower and upper bounds for the x and y bins
            given as (x0, x1, y0, y1) or a tuple of two arrays containing the
            x and y bins.
        field_names: list
            Names of the fields to be used as input data for interpolation.
        truncate: bool
            If False, an exception is raised if the input values are
            outside the bounds of the table.  If True, extrapolation is
            performed.
        store_table: bool
            If True (default), the full table is stored in the interpolator.
            If False, only the shape of the input table is stored and a full
            table must be provided when calling the interpolator.

        Examples
        --------

        >>> ad = ds.all_data()
        >>> table_data = np.random.random((64, 64))
        >>> interp = BilinearFieldInterpolator(table_data, (0.0, 1.0, 0.0, 1.0),
                                              ["x", "y"],
                                              truncate=True)
        >>> field_data = interp(ad)

        If you want to re-use the interpolator with table_data of the same shape
        but different values, you can also supply the `table` keyword argument when
        calling the interpolator:

        >>> new_table_data = np.random.random((64, 64))
        >>> field_data = interp(ad, table=new_table_data)

        """
        super().__init__(table, field_names, truncate=truncate, store_table=store_table)
        self.x_name, self.y_name = field_names
        if len(boundaries) == 4:
            x0, x1, y0, y1 = boundaries
            self.x_bins = np.linspace(x0, x1, table.shape[0], dtype="float64")
            self.y_bins = np.linspace(y0, y1, table.shape[1], dtype="float64")
        elif len(boundaries) == 2:
            self._validate_bin_boundaries(boundaries)
            self.x_bins, self.y_bins = boundaries
        else:
            msg = "Boundaries must be given as (x0, x1, y0, y1) or as (x_bins, y_bins)"
            raise ValueError(msg)

    def __call__(self, data_object, *, table=None):
        table = self._validate_table(table)

        orig_shape = data_object[self.x_name].shape
        x_vals, x_i, y_vals, y_i = self._get_digitized_arrays(data_object)
        my_vals = np.zeros(x_vals.shape, dtype="float64")
        lib.BilinearlyInterpolate(
            table, x_vals, y_vals, self.x_bins, self.y_bins, x_i, y_i, my_vals
        )
        my_vals.shape = orig_shape
        return my_vals


class TrilinearFieldInterpolator(_LinearInterpolator):
    _ndim = 3

    def __init__(
        self, table, boundaries, field_names, truncate=False, *, store_table=True
    ):
        r"""Initialize a 3D interpolator for field data.

        table : array
            The data table over which interpolation is performed.
        boundaries: tuple
            Either a tuple of lower and upper bounds for the x, y, and z bins
            given as (x0, x1, y0, y1, z0, z1) or a tuple of three arrays
            containing the x, y, and z bins.
        field_names: list
            Names of the fields to be used as input data for interpolation.
        truncate: bool
            If False (default), an exception is raised if the input values are
            outside the bounds of the table.  If True, extrapolation is
            performed.
        store_table: bool
            If True (default), the full table is stored in the interpolator.
            If False, only the shape of the input table is stored and a full
            table must be provided when calling the interpolator.

        Examples
        --------

        >>> ad = ds.all_data()
        >>> table_data = np.random.random((64, 64, 64))
        >>> interp = TrilinearFieldInterpolator(table_data,
                                               (0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
                                               ["x", "y", "z"],
                                               truncate=True)
        >>> field_data = interp(ad)

        If you want to re-use the interpolator with table_data of the same shape
        but different values, you can also supply the `table` keyword argument when
        calling the interpolator:

        >>> new_table_data = np.random.random((64, 64, 64))
        >>> field_data = interp(ad, table=new_table_data)

        """
        super().__init__(table, field_names, truncate=truncate, store_table=store_table)
        self.x_name, self.y_name, self.z_name = field_names
        if len(boundaries) == 6:
            x0, x1, y0, y1, z0, z1 = boundaries
            self.x_bins = np.linspace(x0, x1, table.shape[0], dtype="float64")
            self.y_bins = np.linspace(y0, y1, table.shape[1], dtype="float64")
            self.z_bins = np.linspace(z0, z1, table.shape[2], dtype="float64")
        elif len(boundaries) == 3:
            self._validate_bin_boundaries(boundaries)
            self.x_bins, self.y_bins, self.z_bins = boundaries
        else:
            msg = "Boundaries must be given as (x0, x1, y0, y1, z0, z1) or as (x_bins, y_bins, z_bins)"
            raise ValueError(msg)

    def __call__(self, data_object, *, table=None):
        table = self._validate_table(table)

        orig_shape = data_object[self.x_name].shape
        x_vals, x_i, y_vals, y_i, z_vals, z_i = self._get_digitized_arrays(data_object)

        my_vals = np.zeros(x_vals.shape, dtype="float64")
        lib.TrilinearlyInterpolate(
            table,
            x_vals,
            y_vals,
            z_vals,
            self.x_bins,
            self.y_bins,
            self.z_bins,
            x_i,
            y_i,
            z_i,
            my_vals,
        )
        my_vals.shape = orig_shape
        return my_vals


class QuadrilinearFieldInterpolator(_LinearInterpolator):
    _ndim = 4

    def __init__(
        self, table, boundaries, field_names, truncate=False, *, store_table=True
    ):
        r"""Initialize a 4D interpolator for field data.

        table : array
            The data table over which interpolation is performed.
        boundaries: tuple
            Either a tuple of lower and upper bounds for the x, y, z, and w bins
            given as (x0, x1, y0, y1, z0, z1, w0, w1) or a tuple of four arrays
            containing the x, y, z, and w bins.
        field_names: list
            Names of the fields to be used as input data for interpolation.
        truncate: bool
            If False (default), an exception is raised if the input values are
            outside the bounds of the table.  If True, extrapolation is
            performed.
        store_table: bool
            If True (default), the full table is stored in the interpolator.
            If False, only the shape of the input table is stored and a full
            table must be provided when calling the interpolator.

        Examples
        --------
        >>> ad = ds.all_data()
        >>> table_data = np.random.random((64, 64, 64, 64))
        >>> interp = QuadrilinearFieldInterpolator(table_data,
                                                  (0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0),
                                                  ["x", "y", "z", "w"],
                                                  truncate=True)
        >>> field_data = interp(ad)

        If you want to re-use the interpolator with table_data of the same shape
        but different values, you can also supply the `table` keyword argument when
        calling the interpolator:

        >>> new_table_data = np.random.random((64, 64, 64, 64))
        >>> field_data = interp(ad, table=new_table_data)

        """
        super().__init__(table, field_names, truncate=truncate, store_table=store_table)
        self.x_name, self.y_name, self.z_name, self.w_name = field_names
        if len(boundaries) == 8:
            x0, x1, y0, y1, z0, z1, w0, w1 = boundaries
            self.x_bins = np.linspace(x0, x1, table.shape[0]).astype("float64")
            self.y_bins = np.linspace(y0, y1, table.shape[1]).astype("float64")
            self.z_bins = np.linspace(z0, z1, table.shape[2]).astype("float64")
            self.w_bins = np.linspace(w0, w1, table.shape[3]).astype("float64")
        elif len(boundaries) == 4:
            self._validate_bin_boundaries(boundaries)
            self.x_bins, self.y_bins, self.z_bins, self.w_bins = boundaries
        else:
            msg = "Boundaries must be given as (x0, x1, y0, y1, z0, z1, w0, w1) or as (x_bins, y_bins, z_bins, w_bins)"
            raise ValueError(msg)

    def __call__(self, data_object, *, table=None):
        table = self._validate_table(table)

        orig_shape = data_object[self.x_name].shape
        x_vals, x_i, y_vals, y_i, z_vals, z_i, w_vals, w_i = self._get_digitized_arrays(
            data_object
        )

        my_vals = np.zeros(x_vals.shape, dtype="float64")
        lib.QuadrilinearlyInterpolate(
            table,
            x_vals,
            y_vals,
            z_vals,
            w_vals,
            self.x_bins,
            self.y_bins,
            self.z_bins,
            self.w_bins,
            x_i,
            y_i,
            z_i,
            w_i,
            my_vals,
        )
        my_vals.shape = orig_shape
        return my_vals


def get_centers(ds, filename, center_cols, radius_col, unit="1"):
    """
    Return an iterator over EnzoSphere objects generated from the appropriate
    columns in *filename*.  Optionally specify the *unit* radius is in.
    """
    for line in open(filename):
        if line.startswith("#"):
            continue
        vals = line.split()
        x, y, z = (float(vals[i]) for i in center_cols)
        r = float(vals[radius_col])
        yield ds.sphere([x, y, z], r / ds[unit])
