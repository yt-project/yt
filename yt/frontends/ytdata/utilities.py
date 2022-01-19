from yt.units.yt_array import YTArray
from yt.utilities.logger import ytLogger as mylog
from yt.utilities.on_demand_imports import _h5py as h5py


def save_as_dataset(ds, filename, data, field_types=None, extra_attrs=None):
    r"""Export a set of field arrays to a reloadable yt dataset.

    This function can be used to create a yt loadable dataset from a
    set of arrays.  The field arrays can either be associated with a
    loaded dataset or, if not, a dictionary of dataset attributes can
    be provided that will be used as metadata for the new dataset.  The
    resulting dataset can be reloaded as a yt dataset.

    Parameters
    ----------
    ds : dataset or dict
        The dataset associated with the fields or a dictionary of
        parameters.
    filename : str
        The name of the file to be written.
    data : dict
        A dictionary of field arrays to be saved.
    field_types: dict, optional
        A dictionary denoting the group name to which each field is to
        be saved. When the resulting dataset is reloaded, this will be
        the field type for this field. If not given, "data" will be
        used.
    extra_attrs: dict, optional
        A dictionary of additional attributes to be saved.

    Returns
    -------
    filename : str
        The name of the file that has been created.

    Examples
    --------

    >>> import numpy as np
    >>> import yt
    >>> ds = yt.load("enzo_tiny_cosmology/DD0046/DD0046")
    >>> sphere = ds.sphere([0.5] * 3, (10, "Mpc"))
    >>> sphere_density = sphere[("gas", "density")]
    >>> region = ds.box([0.0] * 3, [0.25] * 3)
    >>> region_density = region[("gas", "density")]
    >>> data = {}
    >>> data["sphere_density"] = sphere_density
    >>> data["region_density"] = region_density
    >>> yt.save_as_dataset(ds, "density_data.h5", data)
    >>> new_ds = yt.load("density_data.h5")
    >>> print(new_ds.data["region_density"])
    [  7.47650434e-32   7.70370740e-32   9.74692941e-32 ...,   1.22384547e-27
       5.13889063e-28   2.91811974e-28] g/cm**3
    >>> print(new_ds.data["sphere_density"])
    [  4.46237613e-32   4.86830178e-32   4.46335118e-32 ...,   6.43956165e-30
       3.57339907e-30   2.83150720e-30] g/cm**3
    >>> data = {
    ...     "density": yt.YTArray(1e-24 * np.ones(10), "g/cm**3"),
    ...     "temperature": yt.YTArray(1000.0 * np.ones(10), "K"),
    ... }
    >>> ds_data = {"current_time": yt.YTQuantity(10, "Myr")}
    >>> yt.save_as_dataset(ds_data, "random_data.h5", data)
    >>> new_ds = yt.load("random_data.h5")
    >>> print(new_ds.data[("gas", "temperature")])
    [ 1000.  1000.  1000.  1000.  1000.  1000.  1000.  1000.  1000.  1000.] K

    """

    mylog.info("Saving field data to yt dataset: %s.", filename)

    if extra_attrs is None:
        extra_attrs = {}
    base_attrs = [
        "dimensionality",
        "domain_left_edge",
        "domain_right_edge",
        "current_redshift",
        "current_time",
        "domain_dimensions",
        "geometry",
        "periodicity",
        "cosmological_simulation",
        "omega_lambda",
        "omega_matter",
        "hubble_constant",
        "length_unit",
        "mass_unit",
        "time_unit",
        "velocity_unit",
        "magnetic_unit",
    ]

    fh = h5py.File(filename, mode="w")
    if ds is None:
        ds = {}

    if hasattr(ds, "parameters") and isinstance(ds.parameters, dict):
        for attr, val in ds.parameters.items():
            _yt_array_hdf5_attr(fh, attr, val)

    if hasattr(ds, "unit_registry"):
        _yt_array_hdf5_attr(fh, "unit_registry_json", ds.unit_registry.to_json())

    if hasattr(ds, "unit_system"):
        _yt_array_hdf5_attr(fh, "unit_system_name", ds.unit_system.name.split("_")[0])

    for attr in base_attrs:
        if isinstance(ds, dict):
            my_val = ds.get(attr, None)
        else:
            my_val = getattr(ds, attr, None)
        if my_val is None:
            continue
        _yt_array_hdf5_attr(fh, attr, my_val)

    for attr in extra_attrs:
        my_val = extra_attrs[attr]
        _yt_array_hdf5_attr(fh, attr, my_val)
    if "data_type" not in extra_attrs:
        fh.attrs["data_type"] = "yt_array_data"

    for field in data:
        if field_types is None:
            field_type = "data"
        else:
            field_type = field_types[field]
        if field_type not in fh:
            fh.create_group(field_type)

        if isinstance(field, tuple):
            field_name = field[1]
        else:
            field_name = field

        # for python3
        if data[field].dtype.kind == "U":
            data[field] = data[field].astype("|S")

        _yt_array_hdf5(fh[field_type], field_name, data[field])
        if "num_elements" not in fh[field_type].attrs:
            fh[field_type].attrs["num_elements"] = data[field].size
    fh.close()
    return filename


def _hdf5_yt_array(fh, field, ds=None):
    r"""Load an hdf5 dataset as a YTArray.

    Reads in a dataset from an open hdf5 file or group and uses the
    "units" attribute, if it exists, to apply units.

    Parameters
    ----------
    fh : an open hdf5 file or hdf5 group
        The hdf5 file or group in which the dataset exists.
    field : str
        The name of the field to be loaded.
    ds : yt Dataset
        If not None, the unit_registry of the dataset
        is used to apply units.

    Returns
    -------
    A YTArray of the requested field.

    """

    if ds is None:
        new_arr = YTArray
    else:
        new_arr = ds.arr
    units = ""
    if "units" in fh[field].attrs:
        units = fh[field].attrs["units"]
    if units == "dimensionless":
        units = ""
    return new_arr(fh[field][()], units)


def _yt_array_hdf5(fh, field, data):
    r"""Save a YTArray to an open hdf5 file or group.

    Save a YTArray to an open hdf5 file or group, and save the
    units to a "units" attribute.

    Parameters
    ----------
    fh : an open hdf5 file or hdf5 group
        The hdf5 file or group to which the data will be written.
    field : str
        The name of the field to be saved.
    data : YTArray
        The data array to be saved.

    Returns
    -------
    dataset : hdf5 dataset
        The created hdf5 dataset.

    """

    dataset = fh.create_dataset(str(field), data=data)
    units = ""
    if isinstance(data, YTArray):
        units = str(data.units)
    dataset.attrs["units"] = units
    return dataset


def _yt_array_hdf5_attr(fh, attr, val):
    r"""Save a YTArray or YTQuantity as an hdf5 attribute.

    Save an hdf5 attribute.  If it has units, save an
    additional attribute with the units.

    Parameters
    ----------
    fh : an open hdf5 file, group, or dataset
        The hdf5 file, group, or dataset to which the
        attribute will be written.
    attr : str
        The name of the attribute to be saved.
    val : anything
        The value to be saved.

    """

    if val is None:
        val = "None"
    if hasattr(val, "units"):
        fh.attrs[f"{attr}_units"] = str(val.units)
    try:
        fh.attrs[str(attr)] = val
    # This is raised if no HDF5 equivalent exists.
    # In that case, save its string representation.
    except TypeError:
        fh.attrs[str(attr)] = repr(val)
