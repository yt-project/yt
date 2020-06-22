import tempfile
import os
import shutil
from yt.testing import fake_random_ds, requires_module
from yt.convenience import load
from numpy.testing import \
    assert_equal, \
    assert_allclose
from yt.utilities.on_demand_imports import _astropy
from yt.visualization.fits_image import \
    FITSImageData, FITSProjection, \
    FITSSlice, FITSOffAxisSlice, \
    FITSOffAxisProjection, \
    assert_same_wcs
from yt.visualization.volume_rendering.off_axis_projection import \
    off_axis_projection


@requires_module("astropy")
def test_fits_image():
    tmpdir = tempfile.mkdtemp()
    curdir = os.getcwd()
    os.chdir(tmpdir)

    fields = ("density", "temperature")
    units = ('g/cm**3', 'K',)
    ds = fake_random_ds(64, fields=fields, units=units, nprocs=16,
                        length_unit=100.0)

    prj = ds.proj("density", 2)
    prj_frb = prj.to_frb((0.5, "unitary"), 128)

    fid1 = prj_frb.to_fits_data(fields=[("gas", "density"), ("gas", "temperature")],
                                length_unit="cm")
    fits_prj = FITSProjection(ds, "z", [ds.fields.gas.density,"temperature"],
                              image_res=128, width=(0.5, "unitary"))

    assert_equal(fid1["density"].data, fits_prj["density"].data)
    assert_equal(fid1["temperature"].data, fits_prj["temperature"].data)

    fid1.writeto("fid1.fits", overwrite=True)
    new_fid1 = FITSImageData.from_file("fid1.fits")

    assert_equal(fid1["density"].data, new_fid1["density"].data)
    assert_equal(fid1["temperature"].data, new_fid1["temperature"].data)
    assert_equal(fid1.length_unit, new_fid1.length_unit)
    assert_equal(fid1.time_unit, new_fid1.time_unit)
    assert_equal(fid1.mass_unit, new_fid1.mass_unit)
    assert_equal(fid1.velocity_unit, new_fid1.velocity_unit)
    assert_equal(fid1.magnetic_unit, new_fid1.magnetic_unit)
    assert_equal(fid1.current_time, new_fid1.current_time)

    ds2 = load("fid1.fits")
    ds2.index

    assert ("fits", "density") in ds2.field_list
    assert ("fits", "temperature") in ds2.field_list

    dw_cm = ds2.domain_width.in_units("cm")

    assert dw_cm[0].v == 50.
    assert dw_cm[1].v == 50.

    slc = ds.slice(2, 0.5)
    slc_frb = slc.to_frb((0.5, "unitary"), 128)

    fid2 = slc_frb.to_fits_data(fields=[("gas", "density"), ("gas", "temperature")],
                                length_unit="cm")
    fits_slc = FITSSlice(ds, "z", [("gas", "density"), ("gas", "temperature")],
                         image_res=128, width=(0.5,"unitary"))

    assert_equal(fid2["density"].data, fits_slc["density"].data)
    assert_equal(fid2["temperature"].data, fits_slc["temperature"].data)

    dens_img = fid2.pop("density")
    temp_img = fid2.pop("temperature")

    combined_fid = FITSImageData.from_images([dens_img, temp_img])
    assert_equal(combined_fid.length_unit, dens_img.length_unit)
    assert_equal(combined_fid.time_unit, dens_img.time_unit)
    assert_equal(combined_fid.mass_unit, dens_img.mass_unit)
    assert_equal(combined_fid.velocity_unit, dens_img.velocity_unit)
    assert_equal(combined_fid.magnetic_unit, dens_img.magnetic_unit)
    assert_equal(combined_fid.current_time, dens_img.current_time)

    cut = ds.cutting([0.1, 0.2, -0.9], [0.5, 0.42, 0.6])
    cut_frb = cut.to_frb((0.5, "unitary"), 128)

    fid3 = cut_frb.to_fits_data(fields=[("gas", "density"),
                                        ds.fields.gas.temperature],
                                length_unit="cm")
    fits_cut = FITSOffAxisSlice(ds, [0.1, 0.2, -0.9], ["density", "temperature"],
                                image_res=128, center=[0.5, 0.42, 0.6],
                                width=(0.5, "unitary"))

    assert_equal(fid3["density"].data, fits_cut["density"].data)
    assert_equal(fid3["temperature"].data, fits_cut["temperature"].data)

    fid3.create_sky_wcs([30.,45.], (1.0,"arcsec/kpc"))
    fid3.writeto("fid3.fits", overwrite=True)
    new_fid3 = FITSImageData.from_file("fid3.fits")
    assert_same_wcs(fid3.wcs, new_fid3.wcs)
    assert new_fid3.wcs.wcs.cunit[0] == "deg"
    assert new_fid3.wcs.wcs.cunit[1] == "deg"
    assert new_fid3.wcs.wcs.ctype[0] == "RA---TAN"
    assert new_fid3.wcs.wcs.ctype[1] == "DEC--TAN"

    buf = off_axis_projection(ds, ds.domain_center, [0.1, 0.2, -0.9],
                              0.5, 128, "density").swapaxes(0, 1)
    fid4 = FITSImageData(buf, fields="density", width=100.0)
    fits_oap = FITSOffAxisProjection(ds, [0.1, 0.2, -0.9], "density",
                                     width=(0.5, "unitary"), image_res=128,
                                     depth=(0.5, "unitary"))

    assert_equal(fid4["density"].data, fits_oap["density"].data)

    fid4.create_sky_wcs([30., 45.], (1.0, "arcsec/kpc"), replace_old_wcs=False)
    assert fid4.wcs.wcs.cunit[0] == "cm"
    assert fid4.wcs.wcs.cunit[1] == "cm"
    assert fid4.wcs.wcs.ctype[0] == "linear"
    assert fid4.wcs.wcs.ctype[1] == "linear"
    assert fid4.wcsa.wcs.cunit[0] == "deg"
    assert fid4.wcsa.wcs.cunit[1] == "deg"
    assert fid4.wcsa.wcs.ctype[0] == "RA---TAN"
    assert fid4.wcsa.wcs.ctype[1] == "DEC--TAN"

    cvg = ds.covering_grid(ds.index.max_level, [0.25, 0.25, 0.25],
                           [32, 32, 32], fields=["density", "temperature"])
    fid5 = cvg.to_fits_data(fields=["density", "temperature"])
    assert fid5.dimensionality == 3

    fid5.update_header("density", "time", 0.1)
    fid5.update_header("all", "units", "cgs")

    assert fid5["density"].header["time"] == 0.1
    assert fid5["temperature"].header["units"] == "cgs"
    assert fid5["density"].header["units"] == "cgs"

    fid6 = FITSImageData.from_images(fid5)

    fid5.change_image_name("density", "mass_per_volume")
    assert fid5["mass_per_volume"].name == "mass_per_volume"
    assert fid5["mass_per_volume"].header["BTYPE"] == "mass_per_volume"
    assert "mass_per_volume" in fid5.fields
    assert "mass_per_volume" in fid5.field_units
    assert "density" not in fid5.fields
    assert "density" not in fid5.field_units

    assert "density" in fid6.fields
    assert_equal(fid6["density"].data, fid5["mass_per_volume"].data)

    fid7 = FITSImageData.from_images(fid4)
    fid7.convolve("density", (3.0, "cm"))

    sigma = 3.0/fid7.wcs.wcs.cdelt[0]
    kernel = _astropy.conv.Gaussian2DKernel(x_stddev=sigma)
    data_conv = _astropy.conv.convolve(fid4["density"].data.d, kernel)
    assert_allclose(data_conv, fid7["density"].data.d)

    os.chdir(curdir)
    shutil.rmtree(tmpdir)
