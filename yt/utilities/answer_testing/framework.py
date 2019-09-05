"""
Title: framework.py
Purpose: Contains answer tests that are used by yt's various frontends
Notes:
    Copyright (c) 2013, yt Development Team.
    Distributed under the terms of the Modified BSD License.
    The full license is in the file COPYING.txt, distributed with this
    software.

    * The downside of hashing the dictionaries the way I have is that
    the keys and values get mashed together, meaning that if a test
    fails, there's no way to know if it's because the actual data is
    different or because just the key name has changed. I might want
    to hash them separately.
"""
from collections import OrderedDict
import os
import sys
import tempfile
import zlib

import matplotlib.image as mpimg
import numpy as np
import pytest

from yt.analysis_modules.halo_analysis.api import HaloCatalog
from yt.analysis_modules.halo_mass_function.api import HaloMassFcn
from . import utils


#============================================
#                 AnswerTest
#============================================
@pytest.mark.usefixtures("cli_testing_opts")
class AnswerTest():
    """
    Contains the various answer tests.

    Attributes:
    -----------
        pass

    Methods:
    --------
        pass
    """

    #-----
    # grid_hierarchy_test
    #-----
    def grid_hierarchy_test(self, ds):
        """
        Tests various aspects of the data set's grids.

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        result = OrderedDict()
        result["grid_dimensions"] = ds.index.grid_dimensions
        result["grid_left_edges"] = ds.index.grid_left_edge
        result["grid_right_edges"] = ds.index.grid_right_edge
        result["grid_levels"] = ds.index.grid_levels
        result["grid_particle_count"] = ds.index.grid_particle_count
        # Put result into a hashable form
        s = b''
        for k, v in result.items():
            s += bytes(k.encode('utf-8')) + bytes(v)
        return s

    #-----
    # parentage_relationships_test
    #-----
    def parentage_relationships_test(self, ds):
        """
        Makes sure the nested grids are properly related (I think).

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        result = OrderedDict()
        result["parents"] = []
        result["children"] = []
        for g in ds.index.grids:
            p = g.Parent
            if p is None:
                result["parents"].append(-1)
            elif hasattr(p, "id"):
                result["parents"].append(p.id)
            else:
                result["parents"].append([pg.id for pg in p])
            result["children"].append([c.id for c in g.Children])
        # Check result for compatibility with hashing
        try:
            result = utils.check_result_hashability(result)
        except ValueError:
            print("Could not put result in a hashable form!")
            sys.exit()
        # Put result in a hashable form
        s = b''
        for k, v in result.items():
            # It's possible for v to contain None (if p is None above),
            # and in that case None cannot be converted to bytes, so
            # I'll use -1 as a default value, though it might make more
            # sense to define a global DEFAULT_NONE_AS_INT constant and
            # use that in case using 0 ever needs to be changed for
            # some reason
            v = np.array(v)
            s += bytes(k.encode('utf-8')) + bytes(v)
        return s

    #-----
    # grid_values_test
    #-----
    def grid_values_test(self, ds, field):
        """
        Tests the actual data stored in each grid.

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        result = OrderedDict()
        for g in ds.index.grids:
            result[g.id] = g[field].tostring()
            g.clear_data()
        # Put result into a hashable form
        s = b''
        for k, v in result.items():
            # v is already a bytes array from above. Casting to a bytes
            # array is the same as doing array.tostring(). Also, g.id
            # is an integer. I did not know you could key a dictionary
            # by an int. Anyways, it means no encoding is needed on the
            # key
            s += bytes(k) + v
        return s

    #-----
    # projection_values_test
    #-----
    def projection_values_test(self, ds, axis, field, weight_field, dobj_type):
        """
        Ensures that projections of various fields using various
        weights are unchanged.

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        if dobj_type is not None:
            dobj = utils.create_obj(ds, dobj_type)
        else:
            dobj = None
        if ds.domain_dimensions[axis] == 1:
            # This original returned None, but None can't be converted
            # to a bytes array, so use -1 as a string, since ints can't
            # be converted to bytes either
            return bytes(str(-1).encode('utf-8'))
        proj = ds.proj(field,
                    axis,
                    weight_field=weight_field,
                    data_source=dobj
                )
        # proj.field_data is an instance of the YTFieldData class,
        # which is basically just an alias for dict, as it inherits
        # from dict and does nothing else. As such, it can be converted
        # to a bytestring the same way as the other dicts in these
        # tests. However, in order to ensure that the keys are written
        # to the bytearray in the same order every time, I need an
        # OrderedDict. As such, I'm going to store the entries 
        # alphabetically by key
        #field_data = utils.convert_to_ordered_dict(proj.field_data)
        s = b''
        for k, v in proj.field_data.items():
            # The key is a tuple, hence the join
            s += bytes(''.join(k).encode('utf-8')) + bytes(v)
        return s

    def field_values_test(self, ds, field, obj_type=None, particle_type=False):
        """
        Tests that the average, minimum, and maximum values of a field
        remain unchanged.

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass 
        """
        # If needed build an instance of the dataset type
        obj = utils.create_obj(ds, obj_type)
        determined_field = obj._determine_fields(field)[0]
        # Get the proper weight field depending on if we're looking at
        # particles or not
        if particle_type:
            weight_field = (determined_field[0], "particle_ones")
        else:
            weight_field = ("index", "ones")
        # Get the average, min, and max
        avg = obj.quantities.weighted_average_quantity(
            determined_field,
            weight=weight_field)
        minimum, maximum = obj.quantities.extrema(field)
        # Return as a hashable bytestring
        return np.array([avg, minimum, maximum]).tostring()
    #-----
    # pixelized_projection_values_test
    #-----
    def pixelized_projection_values_test(self, ds, axis, field,
        weight_field=None, dobj_type=None):
        """
        Tests aspects of a pixelized projection.

        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        if dobj_type is not None:
            obj = utils.create_obj(ds, dobj_type)
        else:
            obj = None
        proj = ds.proj(field, axis, weight_field=weight_field, data_source=obj)
        frb = proj.to_frb((1.0, 'unitary'), 256)
        frb[field]
        if weight_field is not None:
            frb[weight_field]
        d = frb.data
        for f in proj.field_data:
            # Sometimes f will be a tuple.
            d["%s_sum" % (f,)] = proj.field_data[f].sum(dtype="float64")
        # Put the dictionary into a hashable form
        s = b''
        for k, v in d.items():
            if isinstance(k, str):
                s += bytes(k.encode('utf-8')) + bytes(v)
            # Tuples can't be encoded, so the string entries are joined
            # together as one string and then encoded
            elif isinstance(k, tuple):
                s += bytes(''.join(k[:]).encode('utf-8')) + bytes(v)
            else:
                err_msg = "Key isn't string or tuple in pixelized " + \
                    "projection values test."
                raise KeyError(err_msg)
        return s

    #-----
    # check_color
    #-----
    def color_conservation_test(self, ds):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        species_names = ds.field_info.species_names
        dd = ds.all_data()
        dens_yt = dd["density"].copy()
        # Enumerate our species here
        for s in sorted(species_names):
            if s == "El": continue
            dens_yt -= dd["%s_density" % s]
        dens_yt -= dd["metal_density"]
        delta_yt = np.abs(dens_yt / dd["density"])
        # Now we compare color conservation to Enzo's color conservation
        dd = ds.all_data()
        dens_enzo = dd["Density"].copy()
        for f in sorted(ds.field_list):
            ff = f[1]
            if not ff.endswith("_Density"):
                continue
            start_strings = ["Electron_", "SFR_", "Forming_Stellar_",
                             "Dark_Matter", "Star_Particle_"]
            if any([ff.startswith(ss) for ss in start_strings]):
                continue
            dens_enzo -= dd[f]
        delta_enzo = np.abs(dens_enzo / dd["Density"])
        np.testing.assert_almost_equal(delta_yt, delta_enzo)

    #-----
    # simulated_halo_mass_function_test
    #-----
    def simulated_halo_mass_function_test(self, ds, finder):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        hc = HaloCatalog(data_ds=ds, finder_method=finder)
        hc.create()
        hmf = HaloMassFcn(halos_ds=hc.halos_ds)
        result = np.empty((2, hmf.masses_sim.size))
        result[0] = hmf.masses_sim.d
        result[1] = hmf.n_cumulative_sim.d
        # Put in hashable form
        s = result[0].tostring() + result[1].tostring()
        return s

    #-----
    # analytic_halo_mass_function_test
    def analytic_halo_mass_function_test(self, ds, fit):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        hmf = HaloMassFcn(simulation_ds=ds, fitting_function=fit)
        result = np.empty((2, hmf.masses_analytic.size))
        result[0] = hmf.masses_analytic.d
        result[1] = hmf.n_cumulative_analytic.d
        # Put in hashable form
        s = result[0].tostring() + result[1].tostring()
        return s

    #-----
    # small_patch_amr
    #-----
    def small_patch_amr(self, ds, fields, weights, axes, ds_objs):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        # Set up keys of test names
        test_keys = ['grid_hierarchy',
            'parentage_relationships',
            'grid_values',
            'projection_values',
            'field_values'
        ]
        # Grid hierarchy test
        gh_hd = self.grid_hierarchy_test(ds)
        # Parentage relationships test
        pr_hd = self.parentage_relationships_test(ds)
        # Grid values, projection values, and field values tests
        gv_hd = b''
        pv_hd = b''
        fv_hd = b''
        for field in fields:
            gv_hd += self.grid_values_test(ds, field)
            for axis in axes:
                for dobj_name in ds_objs:
                    for weight_field in weights:
                        pv_hd += self.projection_values_test(ds,
                            axis,
                            field,
                            weight_field,
                            dobj_name
                        )
                    fv_hd += self.field_values_test(ds,
                        field,
                        dobj_name
                    )
        # Hash the final byte arrays
        hex_digests = [gh_hd, pr_hd, gv_hd, pv_hd, fv_hd]
        hashes = [utils.generate_hash(hd) for hd in hex_digests]
        hash_dict = {}
        for key, value in zip(test_keys, hashes):
            hash_dict[key] = value
        return hash_dict

    #-----
    # big_patch_amr
    #-----
    def big_patch_amr(self, ds, fields, weights, axes, ds_objs):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        # Set up keys of test names
        test_keys = ['grid_hierarchy',
            'parentage_relationships',
            'grid_values',
            'pixelized_projection_values',
        ]
        # Grid hierarchy test
        gh_hd = self.grid_hierarchy_test(ds)
        # Parentage relationships test
        pr_hd = self.parentage_relationships_test(ds)
        # Grid values and pixelized projection values tests
        # For these tests it might be possible for a frontend to not
        # run them for every field, axis, or data source, so the tests
        # are written to work with one instance of each of those. These
        # loops therefore combine the results of each desired combo.
        # This mirrors the way that it was done in the original test
        gv_hd = b''
        ppv_hd = b''
        for field in fields:
            gv_hd += self.grid_values_test(ds, field)
            for axis in axes:
                for dobj_name in ds_objs:
                    for weight_field in weights:
                        ppv_hd += self.pixelized_projection_values_test(ds,
                            axis,
                            field,
                            weight_field,
                            dobj_name
                        )
        # Hash the final byte arrays
        hex_digests = [gh_hd, pr_hd, gv_hd, ppv_hd]
        hashes = [utils.generate_hash(hd) for hd in hex_digests]
        hash_dict = {}
        for key, value in zip(test_keys, hashes):
            hash_dict[key] = value
        return hash_dict

    #-----
    # generic_array_test
    #-----
    def generic_array_test(self, ds_fn, func, args=None, kwargs=None):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        if args is None:
            args = []
        if kwargs is None:
            kwargs = {}
        return func(args, kwargs).tostring()

    #-----
    # sph_answer
    #-----
    def sph_answer(self, ds, ds_str_repr, ds_nparticles, fields):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        # Make sure we're dealing with the right dataset
        assert str(ds) == ds_str_repr
        # Set up keys of test names
        test_keys = ['pixelized_projection_values', 'field_values']
        # Set up by bytestrings for tests
        ppv_hd = b''
        fv_hd = b''
        # Set up data objects to test with
        dso = [None, ("sphere", ("c", (0.1, 'unitary')))]
        dd = ds.all_data()
        assert dd["particle_position"].shape == (ds_nparticles, 3)
        tot = sum(dd[ptype, "particle_position"].shape[0]
                  for ptype in ds.particle_types if ptype != "all")
        # Check
        assert tot == ds_nparticles
        for dobj_name in dso:
            dobj = utils.create_obj(ds, dobj_name)
            s1 = dobj["ones"].sum()
            s2 = sum(mask.sum() for block, mask in dobj.blocks)
            assert s1 == s2
            for field, weight_field in fields.items():
                if field[0] in ds.particle_types:
                    particle_type = True
                else:
                    particle_type = False
                for axis in [0, 1, 2]:
                    if particle_type is False:
                        ppv_hd += self.pixelized_projection_values_test(
                            ds, axis, field, weight_field, dobj_name)
                fv_hd += self.field_values_test(ds, field, dobj_name,
                                      particle_type=particle_type)
        # Hash the final byte arrays
        hex_digests = [ppv_hd, fv_hd]
        hashes = [utils.generate_hash(hd) for hd in hex_digests]
        hash_dict = {}
        for key, value in zip(test_keys, hashes):
            hash_dict[key] = value
        return hash_dict

    #-----
    # yt_field_test
    #-----
    def yt_field_test(self, ds, field, geometric):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        if geometric:
            obj = ds.all_data()
        else:
            obj = ds.data
        return np.array([obj[field].size, obj[field].mean()]).tostring()

    #-----
    # plot_window_attribute_test
    #-----
    def plot_window_attribute_test(self, ds, plot_field, plot_axis, attr_name,
        attr_args, decimals, plot_type='SlicePlot', callback_id='',
        callback_runners=[]):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        plot = utils.create_plot(ds, plot_type, plot_field, plot_axis, {})
        for r in self.callback_runners:
            r(plot_field, plot)
        attr = getattr(plot, attr_name)
        attr(*attr_args[0], **attr_args[1])
        tmpfd, tmpname = tempfile.mkstemp(suffix='.png')
        os.close(tmpfd)
        plot.save(name=tmpname)
        image = mpimg.imread(tmpname)
        os.remove(tmpname)
        return zlib.compress(image.dumps())

    #-----
    # phase_plot_attribute_test
    #-----
    def phase_plot_attribute_test(self, ds_fn, x_field, y_field, z_field,
                 attr_name, attr_args, decimals, plot_type='PhasePlot',
                 plot_kwargs={}):
        """
        Parameters:
        -----------
            pass

        Raises:
        -------
            pass

        Returns:
        --------
            pass
        """
        plot = utils.create_plot2(ds_fn, x_field, y_field,
                                z_field, plot_type, plot_kwargs)
        attr = getattr(plot, attr_name)
        attr(*attr_args[0], **attr_args[1])
        tmpfd, tmpname = tempfile.mkstemp(suffix='.png')
        os.close(tmpfd)
        plot.save(name=tmpname)
        image = mpimg.imread(tmpname)
        os.remove(tmpname)
        return zlib.compress(image.dumps())

    #-----
    # generic_image_test
    #-----
    def generic_image_test(self, img_fname):
        img_data = mpimg.imread(img_fname)
        return zlib.compress(img_data.dumps())
