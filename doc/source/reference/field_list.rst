
.. _field-list:

Field List
==========

This is a list of many of the fields available in yt.  We have attempted to
include most of the fields that are accessible through the plugin system, as 
well as the fields that are known by the frontends, however it is possible to 
generate many more permutations, particularly through vector operations. For 
more information about the fields framework, see :ref:`fields`.

Some fields are recognized by specific frontends only. These are typically 
fields like density and temperature that have their own names and units in 
the different frontend datasets. Often, these fields are aliased to their 
yt-named counterpart fields (typically 'gas' fieldtypes). For example, in 
the ``FLASH`` frontend, the ``dens`` field (i.e. ``(flash, dens)``) is aliased 
to the gas field density (i.e. ``(gas, density)``), similarly ``(flash, velx)`` 
is aliased to ``(gas, velocity_x)``, and so on. In what follows, if a field 
is aliased it will be noted.

Try using the ``ds.field_list`` and ``ds.derived_field_list`` to view the
native and derived fields available for your dataset respectively. For example
to display the native fields in alphabetical order:

.. notebook-cell::

  import yt
  ds = yt.load("Enzo_64/DD0043/data0043")
  for i in sorted(ds.field_list):
    print i

To figure out out what all of the field types here mean, see
:ref:`known-field-types`.

.. contents:: Table of Contents
   :depth: 1
   :local:
   :backlinks: none

.. _yt-fields:

Universal Fields
----------------

('all', 'mesh_id')
^^^^^^^^^^^^^^^^^^

   * Particle Type: True

**Field Source**

.. code-block:: python

      def particle_mesh_ids(field, data):
          pos = data[ptype, coord_name]
          ids = np.zeros(pos.shape[0], dtype="float64") - 1
          # This is float64 in name only.  It will be properly cast inside the
          # deposit operation.
          #_ids = ids.view("float64")
          data.deposit(pos, [ids], method = "mesh_id")
          return data.apply_units(ids, "")
  

('all', 'particle_angular_momentum_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('all', 'particle_angular_momentum_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_angular_momentum_x(field, data):
          return data[ptype, "particle_mass"] * \
                 data[ptype, "particle_specific_angular_momentum_x"]
  

('all', 'particle_angular_momentum_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_angular_momentum_y(field, data):
          return data[ptype, "particle_mass"] * \
                 data[ptype, "particle_specific_angular_momentum_y"]
  

('all', 'particle_angular_momentum_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_angular_momentum_z(field, data):
          return data[ptype, "particle_mass"] * \
                 data[ptype, "particle_specific_angular_momentum_z"]
  

('all', 'particle_ones')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

**Field Source**

.. code-block:: python

      def particle_ones(field, data):
          v = np.ones(data[ptype, mass_name].shape, dtype="float64")
          return data.apply_units(v, field.units)
  

('all', 'particle_position')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{code}~\rm{length}`
   * Particle Type: True

**Field Source**

.. code-block:: python

          def particle_vectors(field, data):
              v = [data[_ptype, name].in_units(field.units)
                    for name in names]
              c = np.column_stack(v)
              return data.apply_units(c, field.units)
  

('all', 'particle_radial_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_velocity_radius(field, data):
          """
          Radial component of the particles' velocity vectors in spherical coords
          based on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          vel = svel
          vel = YTArray([data[ptype, vel % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          vel = vel - np.reshape(bv, (3, 1))
          sphr = get_sph_r_component(vel, theta, phi, normal)
          return sphr
  

('all', 'particle_radius')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_radius(field, data):
          return get_radius(data, "particle_position_")
  

('all', 'particle_specific_angular_momentum')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_specific_angular_momentum(field, data):
          """
          Calculate the angular of a particle velocity.  Returns a vector for each
          particle.
          """
          if data.has_field_parameter("bulk_velocity"):
              bv = data.get_field_parameter("bulk_velocity")
          else: bv = np.zeros(3, dtype=np.float64)
          xv = data[ptype, svel % 'x'] - bv[0]
          yv = data[ptype, svel % 'y'] - bv[1]
          zv = data[ptype, svel % 'z'] - bv[2]
          center = data.get_field_parameter('center')
          coords = YTArray([data[ptype, spos % 'x'],
                             data[ptype, spos % 'y'],
                             data[ptype, spos % 'z']], dtype=np.float64)
          new_shape = tuple([3] + [1]*(len(coords.shape)-1))
          r_vec = coords - np.reshape(center,new_shape)
          v_vec = YTArray([xv,yv,zv], dtype=np.float64)
          return np.cross(r_vec, v_vec, axis=0)
  

('all', 'particle_specific_angular_momentum_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('all', 'particle_specific_angular_momentum_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_specific_angular_momentum_x(field, data):
          if data.has_field_parameter("bulk_velocity"):
              bv = data.get_field_parameter("bulk_velocity")
          else: bv = np.zeros(3, dtype=np.float64)
          center = data.get_field_parameter('center')
          y = data[ptype, spos % "y"] - center[1]
          z = data[ptype, spos % "z"] - center[2]
          yv = data[ptype, svel % "y"] - bv[1]
          zv = data[ptype, svel % "z"] - bv[2]
          return yv*z - zv*y
  

('all', 'particle_specific_angular_momentum_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_specific_angular_momentum_y(field, data):
          if data.has_field_parameter("bulk_velocity"):
              bv = data.get_field_parameter("bulk_velocity")
          else: bv = np.zeros(3, dtype=np.float64)
          center = data.get_field_parameter('center')
          x = data[ptype, spos % "x"] - center[0]
          z = data[ptype, spos % "z"] - center[2]
          xv = data[ptype, svel % "x"] - bv[0]
          zv = data[ptype, svel % "z"] - bv[2]
          return -(xv*z - zv*x)
  

('all', 'particle_specific_angular_momentum_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_specific_angular_momentum_z(field, data):
          if data.has_field_parameter("bulk_velocity"):
              bv = data.get_field_parameter("bulk_velocity")
          else: bv = np.zeros(3, dtype=np.float64)
          center = data.get_field_parameter('center')
          x = data[ptype, spos % "x"] - center[0]
          y = data[ptype, spos % "y"] - center[1]
          xv = data[ptype, svel % "x"] - bv[0]
          yv = data[ptype, svel % "y"] - bv[1]
          return xv*y - yv*x
  

('all', 'particle_spherical_position_phi')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_position_phi(field, data):
          """
          Phi component of the particles' position vectors in spherical coords
          on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          sphp = get_sph_phi_component(pos, phi, normal)
          return sphp
  

('all', 'particle_spherical_position_radius')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_position_radius(field, data):
          """
          Radial component of the particles' position vectors in spherical coords
          on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          sphr = get_sph_r_component(pos, theta, phi, normal)
          return sphr
  

('all', 'particle_spherical_position_theta')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_position_theta(field, data):
          """
          Theta component of the particles' position vectors in spherical coords
          on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          spht = get_sph_theta_component(pos, theta, phi, normal)
          return spht
  

('all', 'particle_spherical_velocity_phi')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_velocity_phi(field, data):
          """
          Phi component of the particles' velocity vectors in spherical coords
          based on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = YTArray([data[ptype, spos % ax] for ax in "xyz"])
          vel = YTArray([data[ptype, svel % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          vel = vel - np.reshape(bv, (3, 1))
          sphp = get_sph_phi_component(vel, phi, normal)
          return sphp
  

('all', 'particle_spherical_velocity_radius')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_velocity_radius(field, data):
          """
          Radial component of the particles' velocity vectors in spherical coords
          based on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          vel = svel
          vel = YTArray([data[ptype, vel % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          vel = vel - np.reshape(bv, (3, 1))
          sphr = get_sph_r_component(vel, theta, phi, normal)
          return sphr
  

('all', 'particle_spherical_velocity_theta')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_velocity_theta(field, data):
          """
          Theta component of the particles' velocity vectors in spherical coords
          based on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          vel = svel
          vel = YTArray([data[ptype, vel % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          vel = vel - np.reshape(bv, (3, 1))
          spht = get_sph_theta_component(vel, theta, phi, normal)
          return spht
  

('all', 'particle_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

          def particle_vectors(field, data):
              v = [data[_ptype, name].in_units(field.units)
                    for name in names]
              c = np.column_stack(v)
              return data.apply_units(c, field.units)
  

('all', 'particle_velocity_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_velocity_magnitude(field, data):
          """ M{|v|} """
          bulk_velocity = data.get_field_parameter("bulk_velocity")
          if bulk_velocity is None:
              bulk_velocity = np.zeros(3)
          return np.sqrt((data[ptype, svel % 'x'] - bulk_velocity[0])**2
                       + (data[ptype, svel % 'y'] - bulk_velocity[1])**2
                       + (data[ptype, svel % 'z'] - bulk_velocity[2])**2 )
  

('deposit', 'all_cic')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def particle_cic(field, data):
          pos = data[ptype, coord_name]
          d = data.deposit(pos, [data[ptype, mass_name]], method = "cic")
          d = data.apply_units(d, data[ptype, mass_name].units)
          d /= data["index", "cell_volume"]
          return d
  

('deposit', 'all_count')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def particle_count(field, data):
          pos = data[ptype, coord_name]
          d = data.deposit(pos, method = "count")
          d = data.ds.arr(d, input_units = "cm**-3")
          return data.apply_units(d, field.units)
  

('deposit', 'all_density')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def particle_density(field, data):
          pos = data[ptype, coord_name]
          mass = data[ptype, mass_name]
          pos.convert_to_units("code_length")
          mass.convert_to_units("code_mass")
          d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
          d = data.ds.arr(d, "code_mass")
          d /= data["index", "cell_volume"]
          return d
  

('deposit', 'all_mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def particle_mass(field, data):
          pos = data[ptype, coord_name]
          pmass = data[ptype, mass_name].in_units(field.units)
          d = data.deposit(pos, [pmass], method = "sum")
          return data.apply_units(d, field.units)
  

('deposit', 'io_cic')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def particle_cic(field, data):
          pos = data[ptype, coord_name]
          d = data.deposit(pos, [data[ptype, mass_name]], method = "cic")
          d = data.apply_units(d, data[ptype, mass_name].units)
          d /= data["index", "cell_volume"]
          return d
  

('deposit', 'io_count')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def particle_count(field, data):
          pos = data[ptype, coord_name]
          d = data.deposit(pos, method = "count")
          d = data.ds.arr(d, input_units = "cm**-3")
          return data.apply_units(d, field.units)
  

('deposit', 'io_density')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def particle_density(field, data):
          pos = data[ptype, coord_name]
          mass = data[ptype, mass_name]
          pos.convert_to_units("code_length")
          mass.convert_to_units("code_mass")
          d = data.deposit(pos, [data[ptype, mass_name]], method = "sum")
          d = data.ds.arr(d, "code_mass")
          d /= data["index", "cell_volume"]
          return d
  

('deposit', 'io_mass')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def particle_mass(field, data):
          pos = data[ptype, coord_name]
          pmass = data[ptype, mass_name].in_units(field.units)
          d = data.deposit(pos, [pmass], method = "sum")
          return data.apply_units(d, field.units)
  

('gas', 'angular_momentum_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('gas', 'angular_momentum_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _angular_momentum_x(field, data):
          return data[ftype, "cell_mass"] \
               * data[ftype, "specific_angular_momentum_x"]
  

('gas', 'angular_momentum_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _angular_momentum_y(field, data):
          return data[ftype, "cell_mass"] \
               * data[ftype, "specific_angular_momentum_y"]
  

('gas', 'angular_momentum_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _angular_momentum_z(field, data):
          return data[ftype, "cell_mass"] \
               * data[ftype, "specific_angular_momentum_z"]
  

('gas', 'averaged_density')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _averaged_field(field, data):
          nx, ny, nz = data[(ftype, basename)].shape
          new_field = data.ds.arr(np.zeros((nx-2, ny-2, nz-2), dtype=np.float64),
                                  (just_one(data[(ftype, basename)]) *
                                   just_one(data[(ftype, weight)])).units)
          weight_field = data.ds.arr(np.zeros((nx-2, ny-2, nz-2), dtype=np.float64),
                                     data[(ftype, weight)].units)
          i_i, j_i, k_i = np.mgrid[0:3, 0:3, 0:3]
  
          for i, j, k in zip(i_i.ravel(), j_i.ravel(), k_i.ravel()):
              sl = [slice(i, nx-(2-i)), slice(j, ny-(2-j)), slice(k, nz-(2-k))]
              new_field += data[(ftype, basename)][sl] * \
                data[(ftype, weight)][sl]
              weight_field += data[(ftype, weight)][sl]
  
          # Now some fancy footwork
          new_field2 = data.ds.arr(np.zeros((nx, ny, nz)), 
                                   data[(ftype, basename)].units)
          new_field2[1:-1, 1:-1, 1:-1] = new_field / weight_field
          return new_field2
  

('gas', 'baroclinic_vorticity_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('gas', 'baroclinic_vorticity_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _baroclinic_vorticity_x(field, data):
          rho2 = data[ftype, "density"].astype(np.float64)**2
          return (data[ftype, "pressure_gradient_y"] *
                  data[ftype, "density_gradient_z"] -
                  data[ftype, "pressure_gradient_z"] *
                  data[ftype, "density_gradient_z"]) / rho2
  

('gas', 'baroclinic_vorticity_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _baroclinic_vorticity_y(field, data):
          rho2 = data[ftype, "density"].astype(np.float64)**2
          return (data[ftype, "pressure_gradient_z"] *
                  data[ftype, "density_gradient_x"] -
                  data[ftype, "pressure_gradient_x"] *
                  data[ftype, "density_gradient_z"]) / rho2
  

('gas', 'baroclinic_vorticity_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _baroclinic_vorticity_z(field, data):
          rho2 = data[ftype, "density"].astype(np.float64)**2
          return (data[ftype, "pressure_gradient_x"] *
                  data[ftype, "density_gradient_y"] -
                  data[ftype, "pressure_gradient_y"] *
                  data[ftype, "density_gradient_x"]) / rho2
  

('gas', 'baryon_overdensity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _baryon_overdensity(field, data):
          if not hasattr(data.ds, "cosmological_simulation") or \
            not data.ds.cosmological_simulation:
              raise NeedsConfiguration("cosmological_simulation", 1)
          omega_baryon = data.get_field_parameter("omega_baryon")
          if omega_baryon is None:
              raise NeedsParameter("omega_baryon")
          co = data.ds.cosmology
          # critical_density(z) ~ omega_lambda + omega_matter * (1 + z)^3
          # mean matter density(z) ~ omega_matter * (1 + z)^3
          return data[ftype, "density"] / omega_baryon / co.critical_density(0.0) / \
            (1.0 + data.ds.current_redshift)**3
  

('gas', 'cell_mass')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _cell_mass(field, data):
          return data[ftype, "density"] * data["index", "cell_volume"]
  

('gas', 'chandra_emissivity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _chandra_emissivity(field, data):
          logT0 = np.log10(data[ftype, "temperature"].to_ndarray().astype(np.float64)) - 7
          # we get rid of the units here since this is a fit and not an 
          # analytical expression
          return data.ds.arr(data[ftype, "number_density"].to_ndarray().astype(np.float64)**2
                             * (10**(- 0.0103 * logT0**8 + 0.0417 * logT0**7
                                     - 0.0636 * logT0**6 + 0.1149 * logT0**5
                                     - 0.3151 * logT0**4 + 0.6655 * logT0**3
                                     - 1.1256 * logT0**2 + 1.0026 * logT0**1
                                     - 0.6984 * logT0)
                               + data[ftype, "metallicity"].to_ndarray() *
                               10**(  0.0305 * logT0**11 - 0.0045 * logT0**10
                                      - 0.3620 * logT0**9  + 0.0513 * logT0**8
                                      + 1.6669 * logT0**7  - 0.3854 * logT0**6
                                      - 3.3604 * logT0**5  + 0.4728 * logT0**4
                                      + 4.5774 * logT0**3  - 2.3661 * logT0**2
                                      - 1.6667 * logT0**1  - 0.2193 * logT0)),
                             "") # add correct units here
  

('gas', 'courant_time_step')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _courant_time_step(field, data):
          t1 = data["index", "dx"] / (data[ftype, "sound_speed"]
                          + np.abs(data[ftype, "velocity_x"]))
          t2 = data["index", "dy"] / (data[ftype, "sound_speed"]
                          + np.abs(data[ftype, "velocity_y"]))
          t3 = data["index", "dz"] / (data[ftype, "sound_speed"]
                          + np.abs(data[ftype, "velocity_z"]))
          tr = np.minimum(np.minimum(t1, t2), t3)
          return tr
  

('gas', 'cutting_plane_velocity_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def _cp_val(field, data):
              vec = data.get_field_parameter("cp_%s_vec" % (ax))
              bv = data.get_field_parameter("bulk_%s" % basename)
              if bv == None: bv = np.zeros(3)
              tr  = (data[xn] - bv[0]) * vec[0]
              tr += (data[yn] - bv[1]) * vec[1]
              tr += (data[zn] - bv[2]) * vec[2]
              return tr
  

('gas', 'cutting_plane_velocity_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def _cp_val(field, data):
              vec = data.get_field_parameter("cp_%s_vec" % (ax))
              bv = data.get_field_parameter("bulk_%s" % basename)
              if bv == None: bv = np.zeros(3)
              tr  = (data[xn] - bv[0]) * vec[0]
              tr += (data[yn] - bv[1]) * vec[1]
              tr += (data[zn] - bv[2]) * vec[2]
              return tr
  

('gas', 'cutting_plane_velocity_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def _cp_val(field, data):
              vec = data.get_field_parameter("cp_%s_vec" % (ax))
              bv = data.get_field_parameter("bulk_%s" % basename)
              if bv == None: bv = np.zeros(3)
              tr  = (data[xn] - bv[0]) * vec[0]
              tr += (data[yn] - bv[1]) * vec[1]
              tr += (data[zn] - bv[2]) * vec[2]
              return tr
  

('gas', 'cylindrical_radial_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _cylindrical_radial(field, data):
          normal = data.get_field_parameter("normal")
          vectors = obtain_rv_vec(data, (xn, yn, zn),
                                  "bulk_%s" % basename)
          theta = resize_vector(data["index", 'cylindrical_theta'], vectors)
          return get_cyl_r_component(vectors, theta, normal)
  

('gas', 'cylindrical_radial_velocity_absolute')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _cylindrical_radial_absolute(field, data):
          return np.abs(_cylindrical_radial(field, data))
  

('gas', 'cylindrical_tangential_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _cylindrical_tangential(field, data):
          normal = data.get_field_parameter("normal")
          vectors = obtain_rv_vec(data, (xn, yn, zn),
                                  "bulk_%s" % basename)
          theta = data["index", 'cylindrical_theta'].copy()
          theta = np.tile(theta, (3,) + (1,)*len(theta.shape))
          return get_cyl_theta_component(vectors, theta, normal)
  

('gas', 'cylindrical_tangential_velocity_absolute')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _cylindrical_tangential_absolute(field, data):
          return np.abs(_cylindrical_tangential(field, data))
  

('gas', 'dark_matter_density')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'density')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'density_gradient_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{4}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('gas', 'density_gradient_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{4}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def func(field, data):
              ds = div_fac * data["index", "dx"]
              f  = data[grad_field][slice_3dr]/ds[slice_3d]
              f -= data[grad_field][slice_3dl]/ds[slice_3d]
              new_field = data.ds.arr(np.zeros_like(data[grad_field], dtype=np.float64),
                                      f.units)
              new_field[slice_3d] = f
              return new_field
  

('gas', 'density_gradient_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{4}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def func(field, data):
              ds = div_fac * data["index", "dx"]
              f  = data[grad_field][slice_3dr]/ds[slice_3d]
              f -= data[grad_field][slice_3dl]/ds[slice_3d]
              new_field = data.ds.arr(np.zeros_like(data[grad_field], dtype=np.float64),
                                      f.units)
              new_field[slice_3d] = f
              return new_field
  

('gas', 'density_gradient_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{4}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def func(field, data):
              ds = div_fac * data["index", "dx"]
              f  = data[grad_field][slice_3dr]/ds[slice_3d]
              f -= data[grad_field][slice_3dl]/ds[slice_3d]
              new_field = data.ds.arr(np.zeros_like(data[grad_field], dtype=np.float64),
                                      f.units)
              new_field[slice_3d] = f
              return new_field
  

('gas', 'di_density')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'dii_density')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'dynamical_time')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dynamical_time(field, data):
          """
          sqrt(3 pi / (16 G rho))
          """
          return np.sqrt(3.0 * np.pi / (16.0 * G * data[ftype, "density"]))
  

('gas', 'entropy')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{keV}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _entropy(field, data):
          mw = data.get_field_parameter("mu")
          if mw is None:
              mw = 1.0
          mw *= mh
          gammam1 = 2./3.
          tr = data[ftype,"kT"] / ((data[ftype, "density"]/mw)**gammam1)
          return data.apply_units(tr, field.units)
  

('gas', 'h2i_density')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'h2ii_density')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'h2m_density')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'hdi_density')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'hei_density')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'heii_density')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'heiii_density')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'hi_density')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'hii_density')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'jeans_mass')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _jeans_mass(field, data):
          MJ_constant = (((5.0 * kboltz) / (G * mh)) ** (1.5)) * \
            (3.0 / (4.0 * np.pi)) ** (0.5)
          u = (MJ_constant * \
               ((data[ftype, "temperature"] /
                 data[ftype, "mean_molecular_weight"])**(1.5)) * \
               (data[ftype, "density"]**(-0.5)))
          return u
  

('gas', 'kT')
^^^^^^^^^^^^^

   * Units: :math:`\rm{keV}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _kT(field, data):
          return (kboltz*data[ftype, "temperature"]).in_units("keV")
  

('gas', 'kinetic_energy')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm} \cdot \rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _kin_energy(field, data):
          return 0.5*data[ftype, "density"] * ( data[ftype, "velocity_x"]**2.0
                                                + data[ftype, "velocity_y"]**2.0
                                                + data[ftype, "velocity_z"]**2.0 )
  

('gas', 'mach_number')
^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _mach_number(field, data):
          """ M{|v|/c_sound} """
          return data[ftype, "velocity_magnitude"] / data[ftype, "sound_speed"]
  

('gas', 'magnetic_energy')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm} \cdot \rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnetic_energy(field,data):
          """This assumes that your front end has provided Bx, By, Bz in
          units of Gauss. If you use MKS, make sure to write your own
          magnetic_energy field to deal with non-unitary \mu_0.
          """
          return (data[ftype,"magnetic_field_x"]**2 +
                  data[ftype,"magnetic_field_y"]**2 +
                  data[ftype,"magnetic_field_z"]**2)/(8*np.pi)
  

('gas', 'magnetic_field_poloidal')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\sqrt{\rm{g}}}{\sqrt{\rm{cm}} \cdot \rm{s}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnetic_field_poloidal(field,data):
          normal = data.get_field_parameter("normal")
          d = data[ftype,'magnetic_field_x']
          Bfields = data.ds.arr(
                      [data[ftype,'magnetic_field_x'],
                       data[ftype,'magnetic_field_y'],
                       data[ftype,'magnetic_field_z']],
                       d.units)
          
          theta = data["index", 'spherical_theta']
          phi   = data["index", 'spherical_phi']
          
          return get_sph_theta_component(Bfields, theta, phi, normal)
  

('gas', 'magnetic_field_strength')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\sqrt{\rm{g}}}{\sqrt{\rm{cm}} \cdot \rm{s}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnetic_field_strength(field,data):
          return np.sqrt(8.*np.pi*data[ftype,"magnetic_energy"])
  

('gas', 'magnetic_field_toroidal')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\sqrt{\rm{g}}}{\sqrt{\rm{cm}} \cdot \rm{s}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnetic_field_toroidal(field,data):
          normal = data.get_field_parameter("normal")
          d = data[ftype,'magnetic_field_x']
          Bfields = data.ds.arr(
                      [data[ftype,'magnetic_field_x'],
                       data[ftype,'magnetic_field_y'],
                       data[ftype,'magnetic_field_z']],
                       d.units)
          
          phi = data["index", 'spherical_phi']
          return get_sph_phi_component(Bfields, phi, normal)
  

('gas', 'magnetic_field_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\sqrt{\rm{g}}}{\sqrt{\rm{cm}} \cdot \rm{s}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'magnetic_field_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\sqrt{\rm{g}}}{\sqrt{\rm{cm}} \cdot \rm{s}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'magnetic_field_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\sqrt{\rm{g}}}{\sqrt{\rm{cm}} \cdot \rm{s}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'magnetic_pressure')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm} \cdot \rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnetic_pressure(field,data):
          return data[ftype,'magnetic_energy']
  

('gas', 'matter_density')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _matter_density(field, data):
          return data[ftype, "density"] + \
            data[ftype, "dark_matter_density"]
  

('gas', 'matter_mass')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _matter_mass(field, data):
          return data[ftype, "matter_density"] * data["index", "cell_volume"]
  

('gas', 'matter_overdensity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _matter_overdensity(field, data):
          if not hasattr(data.ds, "cosmological_simulation") or \
            not data.ds.cosmological_simulation:
              raise NeedsConfiguration("cosmological_simulation", 1)
          co = data.ds.cosmology
          # critical_density(z) ~ omega_lambda + omega_matter * (1 + z)^3
          # mean density(z) ~ omega_matter * (1 + z)^3
          return data[ftype, "matter_density"] / data.ds.omega_matter / \
            co.critical_density(0.0) / \
            (1.0 + data.ds.current_redshift)**3
  

('gas', 'mean_molecular_weight')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _mean_molecular_weight(field, data):
          return (data[ftype, "density"] / (mh * data[ftype, "number_density"]))
  

('gas', 'metal_density')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'metal_mass')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _metal_mass(field, data):
          return data[ftype, "metal_density"] * data["index", "cell_volume"]
  

('gas', 'metallicity')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{Z}_\odot`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _metallicity(field, data):
          tr = data[ftype, "metal_density"] / data[ftype, "density"]
          tr /= metallicity_sun
          return data.apply_units(tr, "Zsun")
  

('gas', 'number_density')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{cm}^{3}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'overdensity')
^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _overdensity(field, data):
          if not hasattr(data.ds, "cosmological_simulation") or \
            not data.ds.cosmological_simulation:
              raise NeedsConfiguration("cosmological_simulation", 1)
          co = data.ds.cosmology
          return data[ftype, "matter_density"] / \
            co.critical_density(data.ds.current_redshift)
  

('gas', 'plasma_beta')
^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _plasma_beta(field,data):
          """This assumes that your front end has provided Bx, By, Bz in
          units of Gauss. If you use MKS, make sure to write your own
          PlasmaBeta field to deal with non-unitary \mu_0.
          """
          return data[ftype,'pressure']/data[ftype,'magnetic_energy']
  

('gas', 'pressure')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm} \cdot \rm{s}^{2}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'pressure_gradient_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{2} \cdot \rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('gas', 'pressure_gradient_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{2} \cdot \rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def func(field, data):
              ds = div_fac * data["index", "dx"]
              f  = data[grad_field][slice_3dr]/ds[slice_3d]
              f -= data[grad_field][slice_3dl]/ds[slice_3d]
              new_field = data.ds.arr(np.zeros_like(data[grad_field], dtype=np.float64),
                                      f.units)
              new_field[slice_3d] = f
              return new_field
  

('gas', 'pressure_gradient_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{2} \cdot \rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def func(field, data):
              ds = div_fac * data["index", "dx"]
              f  = data[grad_field][slice_3dr]/ds[slice_3d]
              f -= data[grad_field][slice_3dl]/ds[slice_3d]
              new_field = data.ds.arr(np.zeros_like(data[grad_field], dtype=np.float64),
                                      f.units)
              new_field[slice_3d] = f
              return new_field
  

('gas', 'pressure_gradient_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{g}}{\rm{cm}^{2} \cdot \rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def func(field, data):
              ds = div_fac * data["index", "dx"]
              f  = data[grad_field][slice_3dr]/ds[slice_3d]
              f -= data[grad_field][slice_3dl]/ds[slice_3d]
              new_field = data.ds.arr(np.zeros_like(data[grad_field], dtype=np.float64),
                                      f.units)
              new_field[slice_3d] = f
              return new_field
  

('gas', 'radial_mach_number')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _radial_mach_number(field, data):
          """ Radial component of M{|v|/c_sound} """
          tr = data[ftype, "radial_velocity"] / data[ftype, "sound_speed"]
          return np.abs(tr)
  

('gas', 'radial_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _radial(field, data):
          normal = data.get_field_parameter("normal")
          vectors = obtain_rv_vec(data, (xn, yn, zn),
                                  "bulk_%s" % basename)
          theta = data['index', 'spherical_theta']
          phi   = data['index', 'spherical_phi']
          rv = get_sph_r_component(vectors, theta, phi, normal)
          # Now, anywhere that radius is in fact zero, we want to zero out our
          # return values.
          rv[np.isnan(theta)] = 0.0
          return rv
  

('gas', 'radial_velocity_absolute')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _radial(field, data):
          normal = data.get_field_parameter("normal")
          vectors = obtain_rv_vec(data, (xn, yn, zn),
                                  "bulk_%s" % basename)
          theta = data['index', 'spherical_theta']
          phi   = data['index', 'spherical_phi']
          rv = get_sph_r_component(vectors, theta, phi, normal)
          # Now, anywhere that radius is in fact zero, we want to zero out our
          # return values.
          rv[np.isnan(theta)] = 0.0
          return rv
  

('gas', 'radiation_acceleration_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{cm}}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'radiation_acceleration_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{cm}}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'radiation_acceleration_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{cm}}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'shear')
^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _shear(field, data):
          """
          Shear is defined as [(dvx/dy + dvy/dx)^2 + (dvz/dy + dvy/dz)^2 +
                               (dvx/dz + dvz/dx)^2 ]^(0.5)
          where dvx/dy = [vx(j-1) - vx(j+1)]/[2dy]
          and is in units of s^(-1)
          (it's just like vorticity except add the derivative pairs instead
           of subtracting them)
          """
          
          if data.ds.dimensionality > 1:
              dvydx = (data[ftype, "velocity_y"][sl_right,sl_center,sl_center] -
                      data[ftype, "velocity_y"][sl_left,sl_center,sl_center]) \
                      / (div_fac*just_one(data["index", "dx"]))
              dvxdy = (data[ftype, "velocity_x"][sl_center,sl_right,sl_center] -
                      data[ftype, "velocity_x"][sl_center,sl_left,sl_center]) \
                      / (div_fac*just_one(data["index", "dy"]))
              f  = (dvydx + dvxdy)**2.0
              del dvydx, dvxdy
          if data.ds.dimensionality > 2:
              dvzdy = (data[ftype, "velocity_z"][sl_center,sl_right,sl_center] -
                      data[ftype, "velocity_z"][sl_center,sl_left,sl_center]) \
                      / (div_fac*just_one(data["index", "dy"]))
              dvydz = (data[ftype, "velocity_y"][sl_center,sl_center,sl_right] -
                      data[ftype, "velocity_y"][sl_center,sl_center,sl_left]) \
                      / (div_fac*just_one(data["index", "dz"]))
              f += (dvzdy + dvydz)**2.0
              del dvzdy, dvydz
              dvxdz = (data[ftype, "velocity_x"][sl_center,sl_center,sl_right] -
                      data[ftype, "velocity_x"][sl_center,sl_center,sl_left]) \
                      / (div_fac*just_one(data["index", "dz"]))
              dvzdx = (data[ftype, "velocity_z"][sl_right,sl_center,sl_center] -
                      data[ftype, "velocity_z"][sl_left,sl_center,sl_center]) \
                      / (div_fac*just_one(data["index", "dx"]))
              f += (dvxdz + dvzdx)**2.0
              del dvxdz, dvzdx
          np.sqrt(f, out=f)
          new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_x"]), f.units)
          new_field[sl_center, sl_center, sl_center] = f
          return new_field
  

('gas', 'shear_criterion')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _shear_criterion(field, data):
          """
          Divide by c_s to leave shear in units of cm**-1, which 
          can be compared against the inverse of the local cell size (1/dx) 
          to determine if refinement should occur.
          """
          
          return data[ftype, "shear"] / data[ftype, "sound_speed"]
  

('gas', 'shear_mach')
^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _shear_mach(field, data):
          """
          Dimensionless shear (shear_mach) is defined nearly the same as shear, 
          except that it is scaled by the local dx/dy/dz and the local sound speed.
          So it results in a unitless quantity that is effectively measuring 
          shear in mach number.  
  
          In order to avoid discontinuities created by multiplying by dx/dy/dz at
          grid refinement boundaries, we also multiply by 2**GridLevel.
  
          Shear (Mach) = [(dvx + dvy)^2 + (dvz + dvy)^2 +
                          (dvx + dvz)^2  ]^(0.5) / c_sound
          """
          
          if data.ds.dimensionality > 1:
              dvydx = (data[ftype, "velocity_y"][sl_right,sl_center,sl_center] -
                       data[ftype, "velocity_y"][sl_left,sl_center,sl_center]) \
                      / div_fac
              dvxdy = (data[ftype, "velocity_x"][sl_center,sl_right,sl_center] -
                       data[ftype, "velocity_x"][sl_center,sl_left,sl_center]) \
                      / div_fac
              f  = (dvydx + dvxdy)**2.0
              del dvydx, dvxdy
          if data.ds.dimensionality > 2:
              dvzdy = (data[ftype, "velocity_z"][sl_center,sl_right,sl_center] -
                       data[ftype, "velocity_z"][sl_center,sl_left,sl_center]) \
                      / div_fac
              dvydz = (data[ftype, "velocity_y"][sl_center,sl_center,sl_right] -
                       data[ftype, "velocity_y"][sl_center,sl_center,sl_left]) \
                      / div_fac
              f += (dvzdy + dvydz)**2.0
              del dvzdy, dvydz
              dvxdz = (data[ftype, "velocity_x"][sl_center,sl_center,sl_right] -
                       data[ftype, "velocity_x"][sl_center,sl_center,sl_left]) \
                      / div_fac
              dvzdx = (data[ftype, "velocity_z"][sl_right,sl_center,sl_center] -
                       data[ftype, "velocity_z"][sl_left,sl_center,sl_center]) \
                      / div_fac
              f += (dvxdz + dvzdx)**2.0
              del dvxdz, dvzdx
          f *= (2.0**data["index", "grid_level"][sl_center, sl_center, sl_center] /
                data[ftype, "sound_speed"][sl_center, sl_center, sl_center])**2.0
          np.sqrt(f, out=f)
          new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_x"]), f.units)
          new_field[sl_center, sl_center, sl_center] = f
          return new_field
  

('gas', 'sound_speed')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _sound_speed(field, data):
          tr = data.ds.gamma * data[ftype, "pressure"] / data[ftype, "density"]
          return np.sqrt(tr)
  

('gas', 'specific_angular_momentum_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('gas', 'specific_angular_momentum_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _specific_angular_momentum_x(field, data):
          xv, yv, zv = obtain_velocities(data, ftype)
          rv = obtain_rvec(data)
          rv = np.rollaxis(rv, 0, len(rv.shape))
          rv = data.ds.arr(rv, input_units = data["index", "x"].units)
          return yv * rv[...,2] - zv * rv[...,1]
  

('gas', 'specific_angular_momentum_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _specific_angular_momentum_y(field, data):
          xv, yv, zv = obtain_velocities(data, ftype)
          rv = obtain_rvec(data)
          rv = np.rollaxis(rv, 0, len(rv.shape))
          rv = data.ds.arr(rv, input_units = data["index", "x"].units)
          return - (xv * rv[...,2] - zv * rv[...,0])
  

('gas', 'specific_angular_momentum_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _specific_angular_momentum_z(field, data):
          xv, yv, zv = obtain_velocities(data, ftype)
          rv = obtain_rvec(data)
          rv = np.rollaxis(rv, 0, len(rv.shape))
          rv = data.ds.arr(rv, input_units = data["index", "x"].units)
          return xv * rv[...,1] - yv * rv[...,0]
  

('gas', 'sz_kinetic')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _sz_kinetic(field, data):
          scale = 0.88 * sigma_thompson / mh / clight
          vel_axis = data.get_field_parameter("axis")
          if vel_axis > 2:
              raise NeedsParameter(["axis"])
          vel = data[ftype, "velocity_%s" % ({0: "x", 1: "y", 2: "z"}[vel_axis])]
          return scale * vel * data[ftype, "density"]
  

('gas', 'szy')
^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _szy(field, data):
          scale = 0.88 / mh * kboltz / (me * clight*clight) * sigma_thompson
          return scale * data[ftype, "density"] * data[ftype, "temperature"]
  

('gas', 'tangential_over_velocity_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _tangential_over_magnitude(field, data):
          tr = data[ftype, "tangential_%s" % basename] / \
               data[ftype, "%s_magnitude" % basename]
          return np.abs(tr)
  

('gas', 'tangential_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _tangential(field, data):
          return np.sqrt(data[ftype, "%s_magnitude" % basename]**2.0
                       - data[ftype, "radial_%s" % basename]**2.0)
  

('gas', 'temperature')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{K}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'thermal_energy')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{\rm{cm}^{2}}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'velocity_divergence')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _divergence(field, data):
          ds = div_fac * just_one(data["index", "dx"])
          f  = data[xn][sl_right,1:-1,1:-1]/ds
          f -= data[xn][sl_left ,1:-1,1:-1]/ds
          ds = div_fac * just_one(data["index", "dy"])
          f += data[yn][1:-1,sl_right,1:-1]/ds
          f -= data[yn][1:-1,sl_left ,1:-1]/ds
          ds = div_fac * just_one(data["index", "dz"])
          f += data[zn][1:-1,1:-1,sl_right]/ds
          f -= data[zn][1:-1,1:-1,sl_left ]/ds
          new_field = data.ds.arr(np.zeros(data[xn].shape, dtype=np.float64),
                                  f.units)        
          new_field[1:-1,1:-1,1:-1] = f
          return new_field
  

('gas', 'velocity_divergence_absolute')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _divergence_abs(field, data):
          return np.abs(data[ftype, "%s_divergence" % basename])
  

('gas', 'velocity_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('gas', 'velocity_x')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'velocity_y')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'velocity_z')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: False

**Field Source**

No source available.

('gas', 'vorticity_growth_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_growth_magnitude(field, data):
          result = np.sqrt(data[ftype, "vorticity_growth_x"]**2 +
                           data[ftype, "vorticity_growth_y"]**2 +
                           data[ftype, "vorticity_growth_z"]**2)
          dot = data.ds.arr(np.zeros(result.shape), "")
          for ax in "xyz":
              dot += (data[ftype, "vorticity_%s" % ax] *
                      data[ftype, "vorticity_growth_%s" % ax]).to_ndarray()
          result = np.sign(dot) * result
          return result
  

('gas', 'vorticity_growth_magnitude_absolute')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_growth_magnitude_absolute(field, data):
          return np.sqrt(data[ftype, "vorticity_growth_x"]**2 +
                         data[ftype, "vorticity_growth_y"]**2 +
                         data[ftype, "vorticity_growth_z"]**2)
  

('gas', 'vorticity_growth_timescale')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_growth_timescale(field, data):
          domegax_dt = data[ftype, "vorticity_x"] / data[ftype, "vorticity_growth_x"]
          domegay_dt = data[ftype, "vorticity_y"] / data[ftype, "vorticity_growth_y"]
          domegaz_dt = data[ftype, "vorticity_z"] / data[ftype, "vorticity_growth_z"]
          return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
  

('gas', 'vorticity_growth_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_growth_x(field, data):
          return -data[ftype, "vorticity_stretching_x"] - \
            data[ftype, "baroclinic_vorticity_x"]
  

('gas', 'vorticity_growth_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_growth_y(field, data):
          return -data[ftype, "vorticity_stretching_y"] - \
            data[ftype, "baroclinic_vorticity_y"]
  

('gas', 'vorticity_growth_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_growth_z(field, data):
          return -data[ftype, "vorticity_stretching_z"] - \
            data[ftype, "baroclinic_vorticity_z"]
  

('gas', 'vorticity_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('gas', 'vorticity_radiation_pressure_growth_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_radiation_pressure_growth_magnitude(field, data):
          result = np.sqrt(data[ftype, "vorticity_radiation_pressure_growth_x"]**2 +
                           data[ftype, "vorticity_radiation_pressure_growth_y"]**2 +
                           data[ftype, "vorticity_radiation_pressure_growth_z"]**2)
          dot = data.ds.arr(np.zeros(result.shape), "")
          for ax in "xyz":
              dot += (data[ftype, "vorticity_%s" % ax] *
                      data[ftype, "vorticity_growth_%s" % ax]).to_ndarray()
          result = np.sign(dot) * result
          return result
  

('gas', 'vorticity_radiation_pressure_growth_magnitude_absolute')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_radiation_pressure_growth_magnitude_absolute(field, data):
          return np.sqrt(data[ftype, "vorticity_radiation_pressure_growth_x"]**2 +
                         data[ftype, "vorticity_radiation_pressure_growth_y"]**2 +
                         data[ftype, "vorticity_radiation_pressure_growth_z"]**2)
  

('gas', 'vorticity_radiation_pressure_growth_timescale')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_radiation_pressure_growth_timescale(field, data):
          domegax_dt = data[ftype, "vorticity_x"] / \
            data[ftype, "vorticity_radiation_pressure_growth_x"]
          domegay_dt = data[ftype, "vorticity_y"] / \
            data[ftype, "vorticity_radiation_pressure_growth_y"]
          domegaz_dt = data[ftype, "vorticity_z"] / \
            data[ftype, "vorticity_radiation_pressure_growth_z"]
          return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
  

('gas', 'vorticity_radiation_pressure_growth_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_radiation_pressure_growth_x(field, data):
          return -data[ftype, "vorticity_stretching_x"] - \
            data[ftype, "baroclinic_vorticity_x"] \
            -data[ftype, "vorticity_radiation_pressure_x"]
  

('gas', 'vorticity_radiation_pressure_growth_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_radiation_pressure_growth_y(field, data):
          return -data[ftype, "vorticity_stretching_y"] - \
            data[ftype, "baroclinic_vorticity_y"] \
            -data[ftype, "vorticity_radiation_pressure_y"]
  

('gas', 'vorticity_radiation_pressure_growth_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_radiation_pressure_growth_z(field, data):
          return -data[ftype, "vorticity_stretching_z"] - \
            data[ftype, "baroclinic_vorticity_z"] \
            -data[ftype, "vorticity_radiation_pressure_z"]
  

('gas', 'vorticity_radiation_pressure_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('gas', 'vorticity_radiation_pressure_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_radiation_pressure_x(field, data):
          rho = data[ftype, "density"].astype(np.float64)
          return (data[ftype, "radiation_acceleration_y"] *
                  data[ftype, "density_gradient_z"] -
                  data[ftype, "radiation_acceleration_z"] *
                  data[ftype, "density_gradient_y"]) / rho
  

('gas', 'vorticity_radiation_pressure_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_radiation_pressure_y(field, data):
          rho = data[ftype, "density"].astype(np.float64)
          return (data[ftype, "radiation_acceleration_z"] *
                  data[ftype, "density_gradient_x"] -
                  data[ftype, "radiation_acceleration_x"] *
                  data[ftype, "density_gradient_z"]) / rho
  

('gas', 'vorticity_radiation_pressure_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_radiation_pressure_z(field, data):
          rho = data[ftype, "density"].astype(np.float64)
          return (data[ftype, "radiation_acceleration_x"] *
                  data[ftype, "density_gradient_y"] -
                  data[ftype, "radiation_acceleration_y"] *
                  data[ftype, "density_gradient_x"]) / rho
  

('gas', 'vorticity_squared')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _squared(field, data):
          squared  = data[xn] * data[xn]
          squared += data[yn] * data[yn]
          squared += data[zn] * data[zn]
          return squared
  

('gas', 'vorticity_stretching_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('gas', 'vorticity_stretching_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_stretching_x(field, data):
          return data[ftype, "velocity_divergence"] * data[ftype, "vorticity_x"]
  

('gas', 'vorticity_stretching_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_stretching_y(field, data):
          return data[ftype, "velocity_divergence"] * data[ftype, "vorticity_y"]
  

('gas', 'vorticity_stretching_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\frac{1}{\rm{s}^{2}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_stretching_z(field, data):
          return data[ftype, "velocity_divergence"] * data[ftype, "vorticity_z"]
  

('gas', 'vorticity_x')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_x(field, data):
          f  = (data[ftype, "velocity_z"][sl_center,sl_right,sl_center] -
                data[ftype, "velocity_z"][sl_center,sl_left,sl_center]) \
                / (div_fac*just_one(data["index", "dy"]).in_cgs())
          f -= (data[ftype, "velocity_y"][sl_center,sl_center,sl_right] -
                data[ftype, "velocity_y"][sl_center,sl_center,sl_left]) \
                / (div_fac*just_one(data["index", "dz"].in_cgs()))
          new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_z"],
                                                dtype=np.float64),
                                  f.units)
          new_field[sl_center, sl_center, sl_center] = f
          return new_field
  

('gas', 'vorticity_y')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_y(field, data):
          f  = (data[ftype, "velocity_x"][sl_center,sl_center,sl_right] -
                data[ftype, "velocity_x"][sl_center,sl_center,sl_left]) \
                / (div_fac*just_one(data["index", "dz"]))
          f -= (data[ftype, "velocity_z"][sl_right,sl_center,sl_center] -
                data[ftype, "velocity_z"][sl_left,sl_center,sl_center]) \
                / (div_fac*just_one(data["index", "dx"]))
          new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_z"],
                                                dtype=np.float64),
                                  f.units)
          new_field[sl_center, sl_center, sl_center] = f
          return new_field
  

('gas', 'vorticity_z')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _vorticity_z(field, data):
          f  = (data[ftype, "velocity_y"][sl_right,sl_center,sl_center] -
                data[ftype, "velocity_y"][sl_left,sl_center,sl_center]) \
                / (div_fac*just_one(data["index", "dx"]))
          f -= (data[ftype, "velocity_x"][sl_center,sl_right,sl_center] -
                data[ftype, "velocity_x"][sl_center,sl_left,sl_center]) \
                / (div_fac*just_one(data["index", "dy"]))
          new_field = data.ds.arr(np.zeros_like(data[ftype, "velocity_z"],
                                                dtype=np.float64),
                                  f.units)
          new_field[sl_center, sl_center, sl_center] = f
          return new_field
  

('gas', 'weak_lensing_convergence')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _weak_lensing_convergence(field, data):
          if not hasattr(data.ds, "cosmological_simulation") or \
            not data.ds.cosmological_simulation:
              raise NeedsConfiguration("cosmological_simulation", 1)
          co = data.ds.cosmology
          observer_redshift = data.get_field_parameter('observer_redshift')
          source_redshift = data.get_field_parameter('source_redshift')
          
          # observer to lens
          dl = co.angular_diameter_distance(observer_redshift, data.ds.current_redshift)
          # observer to source
          ds = co.angular_diameter_distance(observer_redshift, source_redshift)
          # lens to source
          dls = co.angular_diameter_distance(data.ds.current_redshift, source_redshift)
  
          # removed the factor of 1 / a to account for the fact that we are projecting 
          # with a proper distance.
          return (1.5 * (co.hubble_constant / speed_of_light_cgs)**2 * (dl * dls / ds) * \
            data[ftype, "matter_overdensity"]).in_units("1/cm")
  

('gas', 'xray_emissivity')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _xray_emissivity(field, data):
          # old scaling coefficient was 2.168e60
          return data.ds.arr(data[ftype, "density"].to_ndarray().astype(np.float64)**2
                             * data[ftype, "temperature"].to_ndarray()**0.5,
                             "") # add correct units here
  

('index', 'cell_volume')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{code}~\rm{length}^{3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

          def _cell_volume(field, data):
              rv  = data["index", "dx"].copy(order='K')
              rv *= data["index", "dy"]
              rv *= data["index", "dz"]
              return rv
  

('index', 'cylindrical_r')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _cylindrical_r(field, data):
          normal = data.get_field_parameter("normal")
          coords = get_periodic_rvec(data)
          return data.ds.arr(get_cyl_r(coords, normal), "code_length").in_cgs()
  

('index', 'cylindrical_theta')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _cylindrical_theta(field, data):
          normal = data.get_field_parameter("normal")
          coords = get_periodic_rvec(data)
          return get_cyl_theta(coords, normal)
  

('index', 'cylindrical_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _cylindrical_z(field, data):
          normal = data.get_field_parameter("normal")
          coords = get_periodic_rvec(data)
          return data.ds.arr(get_cyl_z(coords, normal), "code_length").in_cgs()
  

('index', 'disk_angle')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _disk_angle(field, data):
          return data["index", "spherical_theta"]
  

('index', 'dx')
^^^^^^^^^^^^^^^

   * Units: :math:`\rm{code}~\rm{length}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dds(field, data):
          rv = data.ds.arr(data.fwidth[...,axi].copy(), units)
          return data._reshape_vals(rv)
  

('index', 'dy')
^^^^^^^^^^^^^^^

   * Units: :math:`\rm{code}~\rm{length}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dds(field, data):
          rv = data.ds.arr(data.fwidth[...,axi].copy(), units)
          return data._reshape_vals(rv)
  

('index', 'dz')
^^^^^^^^^^^^^^^

   * Units: :math:`\rm{code}~\rm{length}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dds(field, data):
          rv = data.ds.arr(data.fwidth[...,axi].copy(), units)
          return data._reshape_vals(rv)
  

('index', 'grid_indices')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _grid_indices(field, data):
          return np.ones(data["index", "ones"].shape)*(data.id-data._id_offset)
  

('index', 'grid_level')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _grid_level(field, data):
          return np.ones(data.ActiveDimensions)*(data.Level)
  

('index', 'height')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _height(field, data):
          return data["index", "cylindrical_z"]
  

('index', 'ones')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _ones(field, data):
          arr = np.ones(data.ires.shape, dtype="float64")
          if data._spatial:
              return data._reshape_vals(arr)
          return data.apply_units(arr, field.units)
  

('index', 'ones_over_dx')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`1 / \rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _ones_over_dx(field, data):
          return np.ones(data["index", "ones"].shape,
                         dtype="float64")/data["index", "dx"]
  

('index', 'radius')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _radius(field, data):
          return get_radius(data, "")
  

('index', 'spherical_phi')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _spherical_phi(field, data):
          normal = data.get_field_parameter("normal")
          coords = get_periodic_rvec(data)
          return get_sph_phi(coords, normal)
  

('index', 'spherical_r')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _spherical_r(field, data):
          coords = get_periodic_rvec(data)
          return data.ds.arr(get_sph_r(coords), "code_length").in_cgs()
  

('index', 'spherical_theta')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _spherical_theta(field, data):
          normal = data.get_field_parameter("normal")
          coords = get_periodic_rvec(data)
          return get_sph_theta(coords, normal)
  

('index', 'x')
^^^^^^^^^^^^^^

   * Units: :math:`\rm{code}~\rm{length}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _coords(field, data):
          rv = data.ds.arr(data.fcoords[...,axi].copy(), units)
          return data._reshape_vals(rv)
  

('index', 'y')
^^^^^^^^^^^^^^

   * Units: :math:`\rm{code}~\rm{length}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _coords(field, data):
          rv = data.ds.arr(data.fcoords[...,axi].copy(), units)
          return data._reshape_vals(rv)
  

('index', 'z')
^^^^^^^^^^^^^^

   * Units: :math:`\rm{code}~\rm{length}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _coords(field, data):
          rv = data.ds.arr(data.fcoords[...,axi].copy(), units)
          return data._reshape_vals(rv)
  

('index', 'zeros')
^^^^^^^^^^^^^^^^^^

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _zeros(field, data):
          arr = np.zeros(data["index", "ones"].shape, dtype='float64')
          return data.apply_units(arr, field.units)
  

('io', 'mesh_id')
^^^^^^^^^^^^^^^^^

   * Particle Type: True

**Field Source**

.. code-block:: python

      def particle_mesh_ids(field, data):
          pos = data[ptype, coord_name]
          ids = np.zeros(pos.shape[0], dtype="float64") - 1
          # This is float64 in name only.  It will be properly cast inside the
          # deposit operation.
          #_ids = ids.view("float64")
          data.deposit(pos, [ids], method = "mesh_id")
          return data.apply_units(ids, "")
  

('io', 'particle_angular_momentum_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('io', 'particle_angular_momentum_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_angular_momentum_x(field, data):
          return data[ptype, "particle_mass"] * \
                 data[ptype, "particle_specific_angular_momentum_x"]
  

('io', 'particle_angular_momentum_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_angular_momentum_y(field, data):
          return data[ptype, "particle_mass"] * \
                 data[ptype, "particle_specific_angular_momentum_y"]
  

('io', 'particle_angular_momentum_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} \cdot \rm{g} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_angular_momentum_z(field, data):
          return data[ptype, "particle_mass"] * \
                 data[ptype, "particle_specific_angular_momentum_z"]
  

('io', 'particle_ones')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

**Field Source**

.. code-block:: python

      def particle_ones(field, data):
          v = np.ones(data[ptype, mass_name].shape, dtype="float64")
          return data.apply_units(v, field.units)
  

('io', 'particle_position')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{code}~\rm{length}`
   * Particle Type: True

**Field Source**

.. code-block:: python

          def particle_vectors(field, data):
              v = [data[_ptype, name].in_units(field.units)
                    for name in names]
              c = np.column_stack(v)
              return data.apply_units(c, field.units)
  

('io', 'particle_radial_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_velocity_radius(field, data):
          """
          Radial component of the particles' velocity vectors in spherical coords
          based on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          vel = svel
          vel = YTArray([data[ptype, vel % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          vel = vel - np.reshape(bv, (3, 1))
          sphr = get_sph_r_component(vel, theta, phi, normal)
          return sphr
  

('io', 'particle_radius')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_radius(field, data):
          return get_radius(data, "particle_position_")
  

('io', 'particle_specific_angular_momentum')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_specific_angular_momentum(field, data):
          """
          Calculate the angular of a particle velocity.  Returns a vector for each
          particle.
          """
          if data.has_field_parameter("bulk_velocity"):
              bv = data.get_field_parameter("bulk_velocity")
          else: bv = np.zeros(3, dtype=np.float64)
          xv = data[ptype, svel % 'x'] - bv[0]
          yv = data[ptype, svel % 'y'] - bv[1]
          zv = data[ptype, svel % 'z'] - bv[2]
          center = data.get_field_parameter('center')
          coords = YTArray([data[ptype, spos % 'x'],
                             data[ptype, spos % 'y'],
                             data[ptype, spos % 'z']], dtype=np.float64)
          new_shape = tuple([3] + [1]*(len(coords.shape)-1))
          r_vec = coords - np.reshape(center,new_shape)
          v_vec = YTArray([xv,yv,zv], dtype=np.float64)
          return np.cross(r_vec, v_vec, axis=0)
  

('io', 'particle_specific_angular_momentum_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _magnitude(field, data):
          mag  = data[xn] * data[xn]
          mag += data[yn] * data[yn]
          mag += data[zn] * data[zn]
          return np.sqrt(mag)
  

('io', 'particle_specific_angular_momentum_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_specific_angular_momentum_x(field, data):
          if data.has_field_parameter("bulk_velocity"):
              bv = data.get_field_parameter("bulk_velocity")
          else: bv = np.zeros(3, dtype=np.float64)
          center = data.get_field_parameter('center')
          y = data[ptype, spos % "y"] - center[1]
          z = data[ptype, spos % "z"] - center[2]
          yv = data[ptype, svel % "y"] - bv[1]
          zv = data[ptype, svel % "z"] - bv[2]
          return yv*z - zv*y
  

('io', 'particle_specific_angular_momentum_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_specific_angular_momentum_y(field, data):
          if data.has_field_parameter("bulk_velocity"):
              bv = data.get_field_parameter("bulk_velocity")
          else: bv = np.zeros(3, dtype=np.float64)
          center = data.get_field_parameter('center')
          x = data[ptype, spos % "x"] - center[0]
          z = data[ptype, spos % "z"] - center[2]
          xv = data[ptype, svel % "x"] - bv[0]
          zv = data[ptype, svel % "z"] - bv[2]
          return -(xv*z - zv*x)
  

('io', 'particle_specific_angular_momentum_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}^{2} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_specific_angular_momentum_z(field, data):
          if data.has_field_parameter("bulk_velocity"):
              bv = data.get_field_parameter("bulk_velocity")
          else: bv = np.zeros(3, dtype=np.float64)
          center = data.get_field_parameter('center')
          x = data[ptype, spos % "x"] - center[0]
          y = data[ptype, spos % "y"] - center[1]
          xv = data[ptype, svel % "x"] - bv[0]
          yv = data[ptype, svel % "y"] - bv[1]
          return xv*y - yv*x
  

('io', 'particle_spherical_position_phi')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_position_phi(field, data):
          """
          Phi component of the particles' position vectors in spherical coords
          on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          sphp = get_sph_phi_component(pos, phi, normal)
          return sphp
  

('io', 'particle_spherical_position_radius')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_position_radius(field, data):
          """
          Radial component of the particles' position vectors in spherical coords
          on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          sphr = get_sph_r_component(pos, theta, phi, normal)
          return sphr
  

('io', 'particle_spherical_position_theta')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_position_theta(field, data):
          """
          Theta component of the particles' position vectors in spherical coords
          on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          spht = get_sph_theta_component(pos, theta, phi, normal)
          return spht
  

('io', 'particle_spherical_velocity_phi')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_velocity_phi(field, data):
          """
          Phi component of the particles' velocity vectors in spherical coords
          based on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = YTArray([data[ptype, spos % ax] for ax in "xyz"])
          vel = YTArray([data[ptype, svel % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          vel = vel - np.reshape(bv, (3, 1))
          sphp = get_sph_phi_component(vel, phi, normal)
          return sphp
  

('io', 'particle_spherical_velocity_radius')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_velocity_radius(field, data):
          """
          Radial component of the particles' velocity vectors in spherical coords
          based on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          vel = svel
          vel = YTArray([data[ptype, vel % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          vel = vel - np.reshape(bv, (3, 1))
          sphr = get_sph_r_component(vel, theta, phi, normal)
          return sphr
  

('io', 'particle_spherical_velocity_theta')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_spherical_velocity_theta(field, data):
          """
          Theta component of the particles' velocity vectors in spherical coords
          based on the provided field parameters for 'normal', 'center', and 
          'bulk_velocity', 
          """
          normal = data.get_field_parameter('normal')
          center = data.get_field_parameter('center')
          bv = data.get_field_parameter("bulk_velocity")
          pos = spos
          pos = YTArray([data[ptype, pos % ax] for ax in "xyz"])
          vel = svel
          vel = YTArray([data[ptype, vel % ax] for ax in "xyz"])
          theta = get_sph_theta(pos, center)
          phi = get_sph_phi(pos, center)
          pos = pos - np.reshape(center, (3, 1))
          vel = vel - np.reshape(bv, (3, 1))
          spht = get_sph_theta_component(vel, theta, phi, normal)
          return spht
  

('io', 'particle_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

          def particle_vectors(field, data):
              v = [data[_ptype, name].in_units(field.units)
                    for name in names]
              c = np.column_stack(v)
              return data.apply_units(c, field.units)
  

('io', 'particle_velocity_magnitude')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\rm{cm} / \rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _particle_velocity_magnitude(field, data):
          """ M{|v|} """
          bulk_velocity = data.get_field_parameter("bulk_velocity")
          if bulk_velocity is None:
              bulk_velocity = np.zeros(3)
          return np.sqrt((data[ptype, svel % 'x'] - bulk_velocity[0])**2
                       + (data[ptype, svel % 'y'] - bulk_velocity[1])**2
                       + (data[ptype, svel % 'z'] - bulk_velocity[2])**2 )
  

.. _ART_specific_fields:

ART-Specific Fields
-------------------

('art', 'Density')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('art', 'TotalEnergy')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``total_energy``
   * Particle Type: False

('art', 'XMomentumDensity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_x``
   * Particle Type: False

('art', 'YMomentumDensity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_y``
   * Particle Type: False

('art', 'ZMomentumDensity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_z``
   * Particle Type: False

('art', 'Pressure')
^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``pressure``
   * Particle Type: False

('art', 'Gamma')
^^^^^^^^^^^^^^^^

   * Aliased to: ``g`` ``a`` ``m`` ``m`` ``a``
   * Particle Type: False

('art', 'GasEnergy')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``thermal_energy``
   * Particle Type: False

('art', 'MetalDensitySNII')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``metal_ii_density``
   * Particle Type: False

('art', 'MetalDensitySNIa')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``metal_ia_density``
   * Particle Type: False

('art', 'PotentialNew')
^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``potential``
   * Particle Type: False

('art', 'PotentialOld')
^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``gas_potential``
   * Particle Type: False

('io', 'particle_position_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_velocity_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Particle Type: True

('io', 'particle_velocity_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Particle Type: True

('io', 'particle_velocity_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Particle Type: True

('io', 'particle_mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

('io', 'particle_index')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

('io', 'particle_species')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_type``
   * Particle Type: True

('io', 'particle_creation_time')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_mass_initial')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

('io', 'particle_metallicity1')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

('io', 'particle_metallicity2')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

.. _ARTIO_specific_fields:

ARTIO-Specific Fields
---------------------

('artio', 'HVAR_GAS_DENSITY')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('artio', 'HVAR_GAS_ENERGY')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``total_energy``
   * Particle Type: False

('artio', 'HVAR_INTERNAL_ENERGY')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``thermal_energy``
   * Particle Type: False

('artio', 'HVAR_PRESSURE')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``pressure``
   * Particle Type: False

('artio', 'HVAR_MOMENTUM_X')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_x``
   * Particle Type: False

('artio', 'HVAR_MOMENTUM_Y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_y``
   * Particle Type: False

('artio', 'HVAR_MOMENTUM_Z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_z``
   * Particle Type: False

('artio', 'HVAR_GAMMA')
^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``gamma``
   * Particle Type: False

('artio', 'HVAR_METAL_DENSITY_Ia')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``metal_ia_density``
   * Particle Type: False

('artio', 'HVAR_METAL_DENSITY_II')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``metal_ii_density``
   * Particle Type: False

('artio', 'VAR_POTENTIAL')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``potential``
   * Particle Type: False

('artio', 'VAR_POTENTIAL_HYDRO')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``gas_potential``
   * Particle Type: False

('io', 'POSITION_X')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position_x``
   * Particle Type: True

('io', 'POSITION_Y')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position_y``
   * Particle Type: True

('io', 'POSITION_Z')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position_z``
   * Particle Type: True

('io', 'VELOCITY_X')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity_x``
   * Particle Type: True

('io', 'VELOCITY_Y')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity_y``
   * Particle Type: True

('io', 'VELOCITY_Z')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity_z``
   * Particle Type: True

('io', 'MASS')
^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

('io', 'PID')
^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

('io', 'SPECIES')
^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_type``
   * Particle Type: True

('io', 'BIRTH_TIME')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Aliased to: ``creation_time``
   * Particle Type: True

('io', 'INITIAL_MASS')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``initial_mass``
   * Particle Type: True

('io', 'METALLICITY_SNIa')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``metallicity_snia``
   * Particle Type: True

('io', 'METALLICITY_SNII')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``metallicity_snii``
   * Particle Type: True

.. _Athena_specific_fields:

Athena-Specific Fields
----------------------

('athena', 'density')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('athena', 'cell_centered_B_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code_magnetic}}`
   * Aliased to: ``magnetic_field_x``
   * Particle Type: False

('athena', 'cell_centered_B_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code_magnetic}}`
   * Aliased to: ``magnetic_field_y``
   * Particle Type: False

('athena', 'cell_centered_B_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code_magnetic}}`
   * Aliased to: ``magnetic_field_z``
   * Particle Type: False

.. _Boxlib_specific_fields:

Boxlib-Specific Fields
----------------------

('boxlib', 'density')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('boxlib', 'eden')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length} \cdot \rm{code}~\rm{time}^{2}}}`
   * Aliased to: ``energy_density``
   * Particle Type: False

('boxlib', 'xmom')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_x``
   * Particle Type: False

('boxlib', 'ymom')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_y``
   * Particle Type: False

('boxlib', 'zmom')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_z``
   * Particle Type: False

('boxlib', 'temperature')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: False

('boxlib', 'Temp')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: False

('boxlib', 'x_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('boxlib', 'y_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('boxlib', 'z_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('boxlib', 'xvel')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('boxlib', 'yvel')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('boxlib', 'zvel')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('io', 'particle_mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

('io', 'particle_position_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_momentum_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Particle Type: True

('io', 'particle_momentum_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Particle Type: True

('io', 'particle_momentum_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Particle Type: True

('io', 'particle_angmomen_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}^{2} / \rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_angmomen_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}^{2} / \rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_angmomen_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}^{2} / \rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_id')
^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

('io', 'particle_mdot')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} / \rm{code}~\rm{time}}`
   * Particle Type: True

.. _Castro_specific_fields:

Castro-Specific Fields
----------------------

('boxlib', 'density')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{g}}{\rm{cm}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('boxlib', 'xmom')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{g}}{\rm{cm}^{2} \cdot \rm{s}}}`
   * Aliased to: ``momentum_x``
   * Particle Type: False

('boxlib', 'ymom')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{g}}{\rm{cm}^{2} \cdot \rm{s}}}`
   * Aliased to: ``momentum_y``
   * Particle Type: False

('boxlib', 'zmom')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{g}}{\rm{cm}^{2} \cdot \rm{s}}}`
   * Aliased to: ``momentum_z``
   * Particle Type: False

('boxlib', 'x_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('boxlib', 'y_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('boxlib', 'z_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('boxlib', 'rho_E')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Aliased to: ``energy_density``
   * Particle Type: False

('boxlib', 'rho_e')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'Temp')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: False

('boxlib', 'grav_x')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{cm}}{\rm{s}^{2}}}`
   * Particle Type: False

('boxlib', 'grav_y')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{cm}}{\rm{s}^{2}}}`
   * Particle Type: False

('boxlib', 'grav_z')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{cm}}{\rm{s}^{2}}}`
   * Particle Type: False

('boxlib', 'pressure')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{dyne}}{\rm{cm}^{2}}}`
   * Particle Type: False

('boxlib', 'kineng')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'soundspeed')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``sound_speed``
   * Particle Type: False

('boxlib', 'Machnumber')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``mach_number``
   * Particle Type: False

('boxlib', 'entropy')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{K} \cdot \rm{g}}}`
   * Aliased to: ``entropy``
   * Particle Type: False

('boxlib', 'magvort')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{1 / \rm{s}}`
   * Aliased to: ``vorticity_magnitude``
   * Particle Type: False

('boxlib', 'divu')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{1 / \rm{s}}`
   * Particle Type: False

('boxlib', 'eint_E')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{erg} / \rm{g}}`
   * Particle Type: False

('boxlib', 'eint_e')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{erg} / \rm{g}}`
   * Particle Type: False

('boxlib', 'magvel')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_magnitude``
   * Particle Type: False

('boxlib', 'radvel')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Particle Type: False

('boxlib', 'magmom')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} \cdot \rm{g} / \rm{s}}`
   * Aliased to: ``momentum_magnitude``
   * Particle Type: False

('boxlib', 'maggrav')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{cm}}{\rm{s}^{2}}}`
   * Particle Type: False

('boxlib', 'phiGrav')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{erg} / \rm{g}}`
   * Particle Type: False

.. _Maestro_specific_fields:

Maestro-Specific Fields
-----------------------

('boxlib', 'density')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{g}}{\rm{cm}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('boxlib', 'x_vel')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('boxlib', 'y_vel')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('boxlib', 'z_vel')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('boxlib', 'magvel')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_magnitude``
   * Particle Type: False

('boxlib', 'radial_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Particle Type: False

('boxlib', 'tfromp')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Particle Type: False

('boxlib', 'tfromh')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Particle Type: False

('boxlib', 'Machnumber')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``mach_number``
   * Particle Type: False

('boxlib', 'S')
^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{1 / \rm{s}}`
   * Particle Type: False

('boxlib', 'ad_excess')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

('boxlib', 'deltaT')
^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

('boxlib', 'deltagamma')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

('boxlib', 'deltap')
^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

('boxlib', 'divw0')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{1 / \rm{s}}`
   * Particle Type: False

('boxlib', 'entropy')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{K} \cdot \rm{g}}}`
   * Aliased to: ``entropy``
   * Particle Type: False

('boxlib', 'entropypert')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

('boxlib', 'enucdot')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{g} \cdot \rm{s}}}`
   * Particle Type: False

('boxlib', 'gpi_x')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{dyne}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'gpi_y')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{dyne}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'gpi_z')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{dyne}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'h')
^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{erg} / \rm{g}}`
   * Particle Type: False

('boxlib', 'h0')
^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{erg} / \rm{g}}`
   * Particle Type: False

('boxlib', 'momentum')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} \cdot \rm{g} / \rm{s}}`
   * Aliased to: ``momentum_magnitude``
   * Particle Type: False

('boxlib', 'p0')
^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'p0pluspi')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'pi')
^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'pioverp0')
^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

('boxlib', 'rho0')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{g}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'rhoh')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Aliased to: ``enthalpy_density``
   * Particle Type: False

('boxlib', 'rhoh0')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'rhohpert')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'rhopert')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{g}}{\rm{cm}^{3}}}`
   * Particle Type: False

('boxlib', 'soundspeed')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``sound_speed``
   * Particle Type: False

('boxlib', 'sponge')
^^^^^^^^^^^^^^^^^^^^

   * Particle Type: False

('boxlib', 'tpert')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Particle Type: False

('boxlib', 'vort')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{1 / \rm{s}}`
   * Aliased to: ``vorticity_magnitude``
   * Particle Type: False

('boxlib', 'w0_x')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Particle Type: False

('boxlib', 'w0_y')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Particle Type: False

('boxlib', 'w0_z')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Particle Type: False

.. _Orion_specific_fields:

Orion-Specific Fields
---------------------

('boxlib', 'density')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('boxlib', 'eden')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length} \cdot \rm{code}~\rm{time}^{2}}}`
   * Aliased to: ``energy_density``
   * Particle Type: False

('boxlib', 'xmom')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_x``
   * Particle Type: False

('boxlib', 'ymom')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_y``
   * Particle Type: False

('boxlib', 'zmom')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Aliased to: ``momentum_z``
   * Particle Type: False

('boxlib', 'temperature')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: False

('boxlib', 'Temp')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: False

('boxlib', 'x_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('boxlib', 'y_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('boxlib', 'z_velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('boxlib', 'xvel')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('boxlib', 'yvel')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('boxlib', 'zvel')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('io', 'particle_mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

('io', 'particle_position_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_momentum_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Particle Type: True

('io', 'particle_momentum_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Particle Type: True

('io', 'particle_momentum_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{2} \cdot \rm{code}~\rm{time}}}`
   * Particle Type: True

('io', 'particle_angmomen_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}^{2} / \rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_angmomen_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}^{2} / \rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_angmomen_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}^{2} / \rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_id')
^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

('io', 'particle_mdot')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} / \rm{code}~\rm{time}}`
   * Particle Type: True

.. _Enzo_specific_fields:

Enzo-Specific Fields
--------------------

('enzo', 'Cooling_Time')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{s}}`
   * Aliased to: ``cooling_time``
   * Particle Type: False

('enzo', 'Dengo_Cooling_Rate')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{g} \cdot \rm{s}}}`
   * Particle Type: False

('enzo', 'Grackle_Cooling_Rate')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3} \cdot \rm{s}}}`
   * Particle Type: False

('enzo', 'HI_kph')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{1 / \rm{code}~\rm{time}}`
   * Particle Type: False

('enzo', 'HeI_kph')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{1 / \rm{code}~\rm{time}}`
   * Particle Type: False

('enzo', 'HeII_kph')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{1 / \rm{code}~\rm{time}}`
   * Particle Type: False

('enzo', 'H2I_kdiss')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{1 / \rm{code}~\rm{time}}`
   * Particle Type: False

('enzo', 'Bx')
^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code_magnetic}}`
   * Aliased to: ``magnetic_field_x``
   * Particle Type: False

('enzo', 'By')
^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code_magnetic}}`
   * Aliased to: ``magnetic_field_y``
   * Particle Type: False

('enzo', 'Bz')
^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code_magnetic}}`
   * Aliased to: ``magnetic_field_z``
   * Particle Type: False

('enzo', 'RadAccel1')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{length}}{\rm{code}~\rm{time}^{2}}}`
   * Aliased to: ``radiation_acceleration_x``
   * Particle Type: False

('enzo', 'RadAccel2')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{length}}{\rm{code}~\rm{time}^{2}}}`
   * Aliased to: ``radiation_acceleration_y``
   * Particle Type: False

('enzo', 'RadAccel3')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{length}}{\rm{code}~\rm{time}^{2}}}`
   * Aliased to: ``radiation_acceleration_z``
   * Particle Type: False

('enzo', 'Dark_Matter_Density')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``dark_matter_density``
   * Particle Type: False

('enzo', 'Temperature')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: False

('enzo', 'Dust_Temperature')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``dust_temperature``
   * Particle Type: False

('enzo', 'x-velocity')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('enzo', 'y-velocity')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('enzo', 'z-velocity')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('enzo', 'RaySegments')
^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``ray_segments``
   * Particle Type: False

('enzo', 'PhotoGamma')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{length}}{\rm{code}~\rm{time}^{2}}}`
   * Aliased to: ``photo_gamma``
   * Particle Type: False

('enzo', 'PotentialField')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``gravitational_potential``
   * Particle Type: False

('enzo', 'Density')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('enzo', 'Metal_Density')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``metal_density``
   * Particle Type: False

('enzo', 'SN_Colour')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Particle Type: False

('enzo', 'Electron_Density')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Particle Type: False

('io', 'particle_position_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_velocity_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity_x``
   * Particle Type: True

('io', 'particle_velocity_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity_y``
   * Particle Type: True

('io', 'particle_velocity_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity_z``
   * Particle Type: True

('io', 'creation_time')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'dynamical_time')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'metallicity_fraction')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{metallicity}}`
   * Particle Type: True

('io', 'metallicity')
^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

('io', 'particle_type')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

('io', 'particle_index')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

('io', 'particle_mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

('io', 'GridID')
^^^^^^^^^^^^^^^^

   * Particle Type: True

('io', 'identifier')
^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

('io', 'level')
^^^^^^^^^^^^^^^

   * Particle Type: True

.. _FLASH_specific_fields:

FLASH-Specific Fields
---------------------

('flash', 'velx')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('flash', 'vely')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('flash', 'velz')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('flash', 'dens')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('flash', 'temp')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{temperature}}`
   * Aliased to: ``temperature``
   * Particle Type: False

('flash', 'pres')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length} \cdot \rm{code}~\rm{time}^{2}}}`
   * Aliased to: ``pressure``
   * Particle Type: False

('flash', 'gpot')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{length}^{2}}{\rm{code}~\rm{time}^{2}}}`
   * Aliased to: ``gravitational_potential``
   * Particle Type: False

('flash', 'gpol')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{length}^{2}}{\rm{code}~\rm{time}^{2}}}`
   * Particle Type: False

('flash', 'tion')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{temperature}}`
   * Particle Type: False

('flash', 'tele')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{temperature}}`
   * Particle Type: False

('flash', 'trad')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{temperature}}`
   * Particle Type: False

('flash', 'pion')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length} \cdot \rm{code}~\rm{time}^{2}}}`
   * Particle Type: False

('flash', 'pele')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length} \cdot \rm{code}~\rm{time}^{2}}}`
   * Particle Type: False

('flash', 'prad')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length} \cdot \rm{code}~\rm{time}^{2}}}`
   * Particle Type: False

('flash', 'eion')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{time}^{2}} \cdot \rm{code}~\rm{length}^{2}}`
   * Particle Type: False

('flash', 'eele')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{time}^{2}} \cdot \rm{code}~\rm{length}^{2}}`
   * Particle Type: False

('flash', 'erad')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{time}^{2}} \cdot \rm{code}~\rm{length}^{2}}`
   * Particle Type: False

('flash', 'pden')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Particle Type: False

('flash', 'depo')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{length}^{2}}{\rm{code}~\rm{time}^{2}}}`
   * Particle Type: False

('flash', 'ye')
^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'magp')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length} \cdot \rm{code}~\rm{time}^{2}}}`
   * Particle Type: False

('flash', 'divb')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} \cdot \rm{code_magnetic}}`
   * Particle Type: False

('flash', 'game')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'gamc')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'flam')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'absr')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'emis')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'cond')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'dfcf')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'fllm')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'pipe')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'tite')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'dbgs')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'cham')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'targ')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'sumy')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'mgdc')
^^^^^^^^^^^^^^^^^

   * Particle Type: False

('flash', 'magx')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code_magnetic}}`
   * Aliased to: ``magnetic_field_x``
   * Particle Type: False

('flash', 'magy')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code_magnetic}}`
   * Aliased to: ``magnetic_field_y``
   * Particle Type: False

('flash', 'magz')
^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code_magnetic}}`
   * Aliased to: ``magnetic_field_z``
   * Particle Type: False

('io', 'particle_posx')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position_x``
   * Particle Type: True

('io', 'particle_posy')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position_y``
   * Particle Type: True

('io', 'particle_posz')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position_z``
   * Particle Type: True

('io', 'particle_velx')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Aliased to: ``particle_velocity_x``
   * Particle Type: True

('io', 'particle_vely')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Aliased to: ``particle_velocity_y``
   * Particle Type: True

('io', 'particle_velz')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Aliased to: ``particle_velocity_z``
   * Particle Type: True

('io', 'particle_tag')
^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

('io', 'particle_mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

.. _GDF_specific_fields:

GDF-Specific Fields
-------------------

('gdf', 'density')
^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{g}}{\rm{cm}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('gdf', 'specific_energy')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{erg} / \rm{g}}`
   * Aliased to: ``thermal_energy``
   * Particle Type: False

('gdf', 'pressure')
^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{erg}}{\rm{cm}^{3}}}`
   * Aliased to: ``pressure``
   * Particle Type: False

('gdf', 'temperature')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: False

('gdf', 'velocity_x')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('gdf', 'velocity_y')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('gdf', 'velocity_z')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('gdf', 'mag_field_x')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{gauss}}`
   * Aliased to: ``magnetic_field_x``
   * Particle Type: False

('gdf', 'mag_field_y')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{gauss}}`
   * Aliased to: ``magnetic_field_y``
   * Particle Type: False

('gdf', 'mag_field_z')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{gauss}}`
   * Aliased to: ``magnetic_field_z``
   * Particle Type: False

.. _HaloCatalog_specific_fields:

HaloCatalog-Specific Fields
---------------------------

(particle_type, 'particle_identifier')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'particle_position_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm}}`
   * Particle Type: True

(particle_type, 'particle_position_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm}}`
   * Particle Type: True

(particle_type, 'particle_position_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm}}`
   * Particle Type: True

(particle_type, 'particle_velocity_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_velocity_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_velocity_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_mass')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{g}}`
   * Particle Type: True

(particle_type, 'virial_radius')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{cm}}`
   * Particle Type: True

.. _OWLSSubfind_specific_fields:

OWLSSubfind-Specific Fields
---------------------------

(particle_type, 'CenterOfMass_0')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Aliased to: ``particle_position_x``
   * Particle Type: True

(particle_type, 'CenterOfMass_1')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Aliased to: ``particle_position_y``
   * Particle Type: True

(particle_type, 'CenterOfMass_2')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Aliased to: ``particle_position_z``
   * Particle Type: True

(particle_type, 'CenterOfMassVelocity_0')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{100000.0 \cdot \rm{m}/(1+z) / \rm{s}}`
   * Aliased to: ``particle_velocity_x``
   * Particle Type: True

(particle_type, 'CenterOfMassVelocity_1')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{100000.0 \cdot \rm{m}/(1+z) / \rm{s}}`
   * Aliased to: ``particle_velocity_y``
   * Particle Type: True

(particle_type, 'CenterOfMassVelocity_2')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{100000.0 \cdot \rm{m}/(1+z) / \rm{s}}`
   * Aliased to: ``particle_velocity_z``
   * Particle Type: True

(particle_type, 'Mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

(particle_type, 'Halo_M_Crit200')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``Virial Mass``
   * Particle Type: True

(particle_type, 'Halo_M_Crit2500')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

(particle_type, 'Halo_M_Crit500')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

(particle_type, 'Halo_M_Mean200')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

(particle_type, 'Halo_M_Mean2500')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

(particle_type, 'Halo_M_Mean500')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

(particle_type, 'Halo_M_TopHat200')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

(particle_type, 'Halo_R_Crit200')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Aliased to: ``Virial Radius``
   * Particle Type: True

(particle_type, 'Halo_R_Crit2500')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'Halo_R_Crit500')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'Halo_R_Mean200')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'Halo_R_Mean2500')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'Halo_R_Mean500')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'Halo_R_TopHat200')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'BH_Mass')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

(particle_type, 'Stars/Mass')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

(particle_type, 'BH_Mdot')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} / \rm{code}~\rm{time}}`
   * Particle Type: True

(particle_type, 'StarFormationRate')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} / \rm{code}~\rm{time}}`
   * Particle Type: True

.. _Rockstar_specific_fields:

Rockstar-Specific Fields
------------------------

(particle_type, 'particle_identifier')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'particle_position_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'particle_position_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'particle_position_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'particle_velocity_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_velocity_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_velocity_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_corevel_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_corevel_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_corevel_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_bulkvel_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_bulkvel_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_bulkvel_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'particle_mass')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{M}_\odot / \rm{h}}`
   * Particle Type: True

(particle_type, 'virial_radius')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'child_r')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'vmax_r')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'mgrav')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'vmax')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'rvmax')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'rs')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'klypin_rs')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{pc}/(1+z) / \rm{h}}`
   * Particle Type: True

(particle_type, 'vrms')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{km} / \rm{s}}`
   * Particle Type: True

(particle_type, 'Jx')
^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'Jy')
^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'Jz')
^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'energy')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'spin')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'alt_m1')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{M}_\odot / \rm{h}}`
   * Particle Type: True

(particle_type, 'alt_m2')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{M}_\odot / \rm{h}}`
   * Particle Type: True

(particle_type, 'alt_m3')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{M}_\odot / \rm{h}}`
   * Particle Type: True

(particle_type, 'alt_m4')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{M}_\odot / \rm{h}}`
   * Particle Type: True

(particle_type, 'Xoff')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'Voff')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'b_to_a')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'c_to_a')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'Ax')
^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'Ay')
^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'Az')
^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'b_to_a2')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'c_to_a2')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'A2x')
^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'A2y')
^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'A2z')
^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'bullock_spin')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'kin_to_pot')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'm_pe_b')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'm_pe_d')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'num_p')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'num_child_particles')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'p_start')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'desc')
^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'flags')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'n_core')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'min_pos_err')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'min_vel_err')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

(particle_type, 'min_bulkvel_err')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

.. _RAMSES_specific_fields:

RAMSES-Specific Fields
----------------------

('ramses', 'Density')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: False

('ramses', 'x-velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Aliased to: ``velocity_x``
   * Particle Type: False

('ramses', 'y-velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Aliased to: ``velocity_y``
   * Particle Type: False

('ramses', 'z-velocity')
^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Aliased to: ``velocity_z``
   * Particle Type: False

('ramses', 'Pressure')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length} \cdot \rm{code}~\rm{time}^{2}}}`
   * Aliased to: ``pressure``
   * Particle Type: False

('ramses', 'Metallicity')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``metallicity``
   * Particle Type: False

('io', 'particle_position_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_position_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

('io', 'particle_velocity_x')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_velocity_y')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_velocity_z')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length} / \rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Particle Type: True

('io', 'particle_identifier')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

('io', 'particle_refinement_level')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

('io', 'particle_age')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Particle Type: True

('io', 'particle_metallicity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Particle Type: True

.. _Gadget_specific_fields:

Gadget-Specific Fields
----------------------

(particle_type, 'Mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

(particle_type, 'Masses')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

(particle_type, 'Coordinates')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position``
   * Particle Type: True

(particle_type, 'Velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity``
   * Particle Type: True

(particle_type, 'Velocities')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity``
   * Particle Type: True

(particle_type, 'ParticleIDs')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

(particle_type, 'InternalEnergy')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``thermal_energy``
   * Particle Type: True

(particle_type, 'SmoothingLength')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``smoothing_length``
   * Particle Type: True

(particle_type, 'Density')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: True

(particle_type, 'MaximumTemperature')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Particle Type: True

(particle_type, 'Temperature')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: True

(particle_type, 'Epsilon')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

(particle_type, 'Metals')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{metallicity}}`
   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Metallicity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{metallicity}}`
   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Phi')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

(particle_type, 'FormationTime')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Aliased to: ``creation_time``
   * Particle Type: True

(particle_type, 'Metallicity_00')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Metallicity_01')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``He_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_02')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``C_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_03')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``N_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_04')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``O_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_05')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Ne_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_06')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Mg_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_07')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Si_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_08')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``S_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_09')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Ca_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_10')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Fe_fraction``
   * Particle Type: True

.. _GadgetHDF5_specific_fields:

GadgetHDF5-Specific Fields
--------------------------

(particle_type, 'Mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

(particle_type, 'Masses')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

(particle_type, 'Coordinates')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position``
   * Particle Type: True

(particle_type, 'Velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity``
   * Particle Type: True

(particle_type, 'Velocities')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity``
   * Particle Type: True

(particle_type, 'ParticleIDs')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

(particle_type, 'InternalEnergy')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``thermal_energy``
   * Particle Type: True

(particle_type, 'SmoothingLength')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``smoothing_length``
   * Particle Type: True

(particle_type, 'Density')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: True

(particle_type, 'MaximumTemperature')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Particle Type: True

(particle_type, 'Temperature')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: True

(particle_type, 'Epsilon')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

(particle_type, 'Metals')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{metallicity}}`
   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Metallicity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{metallicity}}`
   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Phi')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

(particle_type, 'FormationTime')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Aliased to: ``creation_time``
   * Particle Type: True

(particle_type, 'Metallicity_00')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Metallicity_01')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``He_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_02')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``C_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_03')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``N_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_04')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``O_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_05')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Ne_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_06')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Mg_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_07')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Si_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_08')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``S_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_09')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Ca_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_10')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Fe_fraction``
   * Particle Type: True

.. _OWLS_specific_fields:

OWLS-Specific Fields
--------------------

(particle_type, 'Mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

(particle_type, 'Masses')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

(particle_type, 'Coordinates')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position``
   * Particle Type: True

(particle_type, 'Velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity``
   * Particle Type: True

(particle_type, 'Velocities')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity``
   * Particle Type: True

(particle_type, 'ParticleIDs')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

(particle_type, 'InternalEnergy')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``thermal_energy``
   * Particle Type: True

(particle_type, 'SmoothingLength')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``smoothing_length``
   * Particle Type: True

(particle_type, 'Density')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: True

(particle_type, 'MaximumTemperature')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Particle Type: True

(particle_type, 'Temperature')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: True

(particle_type, 'Epsilon')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

(particle_type, 'Metals')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{metallicity}}`
   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Metallicity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{metallicity}}`
   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Phi')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

(particle_type, 'FormationTime')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Aliased to: ``creation_time``
   * Particle Type: True

(particle_type, 'Metallicity_00')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Metallicity_01')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``He_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_02')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``C_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_03')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``N_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_04')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``O_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_05')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Ne_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_06')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Mg_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_07')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Si_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_08')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``S_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_09')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Ca_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_10')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Fe_fraction``
   * Particle Type: True

.. _Tipsy_specific_fields:

Tipsy-Specific Fields
---------------------

(particle_type, 'Mass')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

(particle_type, 'Masses')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``particle_mass``
   * Particle Type: True

(particle_type, 'Coordinates')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``particle_position``
   * Particle Type: True

(particle_type, 'Velocity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity``
   * Particle Type: True

(particle_type, 'Velocities')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``particle_velocity``
   * Particle Type: True

(particle_type, 'ParticleIDs')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``particle_index``
   * Particle Type: True

(particle_type, 'InternalEnergy')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``thermal_energy``
   * Particle Type: True

(particle_type, 'SmoothingLength')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Aliased to: ``smoothing_length``
   * Particle Type: True

(particle_type, 'Density')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\frac{\rm{code}~\rm{mass}}{\rm{code}~\rm{length}^{3}}}`
   * Aliased to: ``density``
   * Particle Type: True

(particle_type, 'MaximumTemperature')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Particle Type: True

(particle_type, 'Temperature')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{K}}`
   * Aliased to: ``temperature``
   * Particle Type: True

(particle_type, 'Epsilon')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

(particle_type, 'Metals')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{metallicity}}`
   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Metallicity')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{metallicity}}`
   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Phi')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{length}}`
   * Particle Type: True

(particle_type, 'FormationTime')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Aliased to: ``creation_time``
   * Particle Type: True

(particle_type, 'Metallicity_00')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``metallicity``
   * Particle Type: True

(particle_type, 'Metallicity_01')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``He_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_02')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``C_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_03')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``N_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_04')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``O_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_05')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Ne_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_06')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Mg_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_07')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Si_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_08')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``S_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_09')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Ca_fraction``
   * Particle Type: True

(particle_type, 'Metallicity_10')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Aliased to: ``Fe_fraction``
   * Particle Type: True

(particle_type, 'uDotFB')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``uDotFB``
   * Particle Type: True

(particle_type, 'acc')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity} / \rm{code}~\rm{time}}`
   * Aliased to: ``acc``
   * Particle Type: True

(particle_type, 'c')
^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity}}`
   * Aliased to: ``c``
   * Particle Type: True

(particle_type, 'uDotDiff')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``uDotDiff``
   * Particle Type: True

(particle_type, 'FeMassFrac')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{}`
   * Aliased to: ``FeMassFrac``
   * Particle Type: True

(particle_type, 'timeform')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Aliased to: ``timeform``
   * Particle Type: True

(particle_type, 'uDot')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``uDot``
   * Particle Type: True

(particle_type, 'accg')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{velocity} / \rm{code}~\rm{time}}`
   * Aliased to: ``accg``
   * Particle Type: True

(particle_type, 'uDotHydro')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``uDotHydro``
   * Particle Type: True

(particle_type, 'uDotPdV')
^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``uDotPdV``
   * Particle Type: True

(particle_type, 'uDotAV')
^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass} \cdot \rm{code}~\rm{velocity}^{2}}`
   * Aliased to: ``uDotAV``
   * Particle Type: True

(particle_type, 'HI')
^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{}`
   * Aliased to: ``HI``
   * Particle Type: True

(particle_type, 'coolontime')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{time}}`
   * Aliased to: ``coolontime``
   * Particle Type: True

(particle_type, 'HeII')
^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{}`
   * Aliased to: ``HeII``
   * Particle Type: True

(particle_type, 'HII')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{}`
   * Aliased to: ``HII``
   * Particle Type: True

(particle_type, 'massform')
^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{\rm{code}~\rm{mass}}`
   * Aliased to: ``massform``
   * Particle Type: True

(particle_type, 'OxMassFrac')
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{}`
   * Aliased to: ``OxMassFrac``
   * Particle Type: True

(particle_type, 'HeI')
^^^^^^^^^^^^^^^^^^^^^^

   * Units: :math:`\mathrm{}`
   * Aliased to: ``HeI``
   * Particle Type: True



Index of Fields
---------------

.. contents:: 
   :depth: 3
   :backlinks: none


