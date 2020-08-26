from yt.testing import assert_array_equal, fake_amr_ds, fake_particle_ds, fake_random_ds


def test_center_squeeze():
    # tests that selected values match when supplying center arrays of different shapes
    # to the data container.

    # list of fields to populate fake datasets with
    fldz = ("Density",)

    # create and test amr, random and particle data
    check_single_ds(fake_amr_ds(fields=fldz), True)
    check_single_ds(fake_random_ds(16, fields=fldz), True)
    check_single_ds(fake_particle_ds(npart=100), False)


def check_single_ds(ds, check_morton):
    # compares values for range of data containers using different center array shapes

    center = ds.domain_center  # reference center array

    # build some data containers
    sp0 = ds.sphere(center, 0.25)
    sl0 = ds.slice(0, 0.25, center=center)
    reg0 = ds.region(center, [-0.25, -0.25, -0.25], [0.25, 0.25, 0.25])

    # store morton indices of each
    i_sp0 = None
    i_sl0 = None
    i_reg0 = None
    if check_morton:
        i_sp0 = sp0["index", "morton_index"]
        i_sp0.sort()
        i_sl0 = sl0["index", "morton_index"]
        i_sl0.sort()
        i_reg0 = reg0["index", "morton_index"]
        i_reg0.sort()

    # create new containers for different shapes of the center array
    for test_shape in [(1, 3), (1, 1, 3)]:
        new_center = center.reshape(test_shape)
        sp = ds.sphere(new_center, 0.25)
        sl = ds.slice(0, 0.25, center=new_center)
        reg = ds.region(new_center, [-0.25, -0.25, -0.25], [0.25, 0.25, 0.25])

        # compare each to the reference containers
        for ob, ob0, i_ob0 in [(sp, sp0, i_sp0), (sl, sl0, i_sl0), (reg, reg0, i_reg0)]:

            # check that selection field values match the reference
            for fld in ds.field_list:
                assert_array_equal(ob[fld], ob0[fld])

            if check_morton:
                # check that morton indices match the reference
                i_ob = ob["index", "morton_index"]
                i_ob.sort()
                assert_array_equal(i_ob, i_ob0)
