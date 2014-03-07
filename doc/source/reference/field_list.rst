
.. _field-list:

Field List
==========

This is a list of all fields available in ``yt``.  It has been organized by the
type of code that each field is supported by.  "Universal" fields are available
everywhere, "Enzo" fields in Enzo datasets, "Orion" fields in Orion datasets,
and so on.

Try using the ``pf.field_list`` and ``pf.derived_field_list`` to view the
native and derived fields available for your dataset respectively. For example
to display the native fields in alphabetical order:

.. notebook-cell::

  from yt.mods import *
  pf = load("Enzo_64/DD0043/data0043")
  for i in sorted(pf.field_list):
    print i

.. note:: Universal fields will be overridden by a code-specific field.

.. rubric:: Table of Contents

.. contents::
   :depth: 2
   :local:
   :backlinks: none

.. _universal-field-list:

Universal Field List
--------------------

AbsDivV
+++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _AbsDivV(field, data):
      return np.abs(data['DivV'])
  

**Convert Function Source**

No source available.

AngularMomentumX
++++++++++++++++

   * Units: :math:`\rm{g}\/\rm{cm}^2/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _AngularMomentumX(field, data):
      return data["CellMass"] * data["SpecificAngularMomentumX"]
  

**Convert Function Source**

No source available.

AngularMomentumY
++++++++++++++++

   * Units: :math:`\rm{g}\/\rm{cm}^2/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _AngularMomentumY(field, data):
      return data["CellMass"] * data["SpecificAngularMomentumY"]
  

**Convert Function Source**

No source available.

AngularMomentumZ
++++++++++++++++

   * Units: :math:`\rm{g}\/\rm{cm}^2/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _AngularMomentumZ(field, data):
      return data["CellMass"] * data["SpecificAngularMomentumZ"]
  

**Convert Function Source**

No source available.

AveragedDensity
+++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _AveragedDensity(field, data):
      nx, ny, nz = data["density"].shape
      new_field = np.zeros((nx-2,ny-2,nz-2), dtype='float64')
      weight_field = np.zeros((nx-2,ny-2,nz-2), dtype='float64')
      i_i, j_i, k_i = np.mgrid[0:3,0:3,0:3]
      for i,j,k in zip(i_i.ravel(),j_i.ravel(),k_i.ravel()):
          sl = [slice(i,nx-(2-i)),slice(j,ny-(2-j)),slice(k,nz-(2-k))]
          new_field += data["density"][sl] * data["CellMass"][sl]
          weight_field += data["CellMass"][sl]
      # Now some fancy footwork
      new_field2 = np.zeros((nx,ny,nz))
      new_field2[1:-1,1:-1,1:-1] = new_field/weight_field
      return new_field2
  

**Convert Function Source**

No source available.

BMagnitude
++++++++++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _BMagnitude(field,data):
      """This assumes that your front end has provided Bx, By, Bz in
      units of Gauss. If you use MKS, make sure to write your own
      BMagnitude field to deal with non-unitary \mu_0.
      """
      return np.sqrt((data["Bx"]**2 + data["By"]**2 + data["Bz"]**2))
  

**Convert Function Source**

No source available.

BPoloidal
+++++++++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _BPoloidal(field,data):
      normal = data.get_field_parameter("normal")
  
      Bfields = np.array([data['Bx'], data['By'], data['Bz']])
  
      theta = data['sph_theta']
      phi   = data['sph_phi']
  
      return get_sph_theta_component(Bfields, theta, phi, normal)
  

**Convert Function Source**

No source available.

BRadial
+++++++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _BRadial(field,data):
      normal = data.get_field_parameter("normal")
  
      Bfields = np.array([data['Bx'], data['By'], data['Bz']])
  
      theta = data['sph_theta']
      phi   = data['sph_phi']
  
      return get_sph_r_component(Bfields, theta, phi, normal)
  

**Convert Function Source**

No source available.

BToroidal
+++++++++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _BToroidal(field,data):
      normal = data.get_field_parameter("normal")
  
      Bfields = np.array([data['Bx'], data['By'], data['Bz']])
  
      phi   = data['sph_phi']
  
      return get_sph_phi_component(Bfields, phi, normal)
  

**Convert Function Source**

No source available.

BaroclinicVorticityMagnitude
++++++++++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _BaroclinicVorticityMagnitude(field, data):
      return np.sqrt(data["BaroclinicVorticityX"]**2 +
                     data["BaroclinicVorticityY"]**2 +
                     data["BaroclinicVorticityZ"]**2)
  

**Convert Function Source**

No source available.

BaroclinicVorticityX
++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _BaroclinicVorticityX(field, data):
      rho2 = data["density"].astype('float64')**2
      return (data["gradPressureY"] * data["gradDensityZ"] -
              data["gradPressureZ"] * data["gradDensityY"]) / rho2
  

**Convert Function Source**

No source available.

BaroclinicVorticityY
++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _BaroclinicVorticityY(field, data):
      rho2 = data["density"].astype('float64')**2
      return (data["gradPressureZ"] * data["gradDensityX"] -
              data["gradPressureX"] * data["gradDensityZ"]) / rho2
  

**Convert Function Source**

No source available.

BaroclinicVorticityZ
++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _BaroclinicVorticityZ(field, data):
      rho2 = data["density"].astype('float64')**2
      return (data["gradPressureX"] * data["gradDensityY"] -
              data["gradPressureY"] * data["gradDensityX"]) / rho2
  

**Convert Function Source**

No source available.

Baryon_Overdensity
++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Baryon_Overdensity(field, data):
      if data.pf.has_key('omega_baryon_now'):
          omega_baryon_now = data.pf['omega_baryon_now']
      else:
          omega_baryon_now = 0.0441
      return data['density'] / (omega_baryon_now * rho_crit_now * 
                                (data.pf.hubble_constant**2) * 
                                ((1+data.pf.current_redshift)**3))
  

**Convert Function Source**

No source available.

CellMass
++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CellMass(field, data):
      return data["density"] * data["CellVolume"]
  

**Convert Function Source**

No source available.

CellMassCode
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CellMassCode(field, data):
      return data["density"] * data["CellVolumeCode"]
  

**Convert Function Source**

.. code-block:: python

  def _convertCellMassCode(data):
      return 1.0/data.convert("density")
  

CellMassMsun
++++++++++++

   * Units: :math:`M_{\odot}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CellMass(field, data):
      return data["density"] * data["CellVolume"]
  

**Convert Function Source**

.. code-block:: python

  def _convertCellMassMsun(data):
      return 5.027854e-34 # g^-1
  

CellVolume
++++++++++

   * Units: :math:`\rm{cm}^3`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CellVolume(field, data):
      if data['dx'].size == 1:
          try:
              return data['dx'] * data['dy'] * data['dz'] * \
                  np.ones(data.ActiveDimensions, dtype='float64')
          except AttributeError:
              return data['dx'] * data['dy'] * data['dz']
      return data["dx"] * data["dy"] * data["dz"]
  

**Convert Function Source**

.. code-block:: python

  def _ConvertCellVolumeCGS(data):
      return data.convert("cm")**3.0
  

CellVolumeCode
++++++++++++++

   * Units: :math:`\rm{BoxVolume}^3`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CellVolume(field, data):
      if data['dx'].size == 1:
          try:
              return data['dx'] * data['dy'] * data['dz'] * \
                  np.ones(data.ActiveDimensions, dtype='float64')
          except AttributeError:
              return data['dx'] * data['dy'] * data['dz']
      return data["dx"] * data["dy"] * data["dz"]
  

**Convert Function Source**

No source available.

CellVolumeMpc
+++++++++++++

   * Units: :math:`\rm{Mpc}^3`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CellVolume(field, data):
      if data['dx'].size == 1:
          try:
              return data['dx'] * data['dy'] * data['dz'] * \
                  np.ones(data.ActiveDimensions, dtype='float64')
          except AttributeError:
              return data['dx'] * data['dy'] * data['dz']
      return data["dx"] * data["dy"] * data["dz"]
  

**Convert Function Source**

.. code-block:: python

  def _ConvertCellVolumeMpc(data):
      return data.convert("mpc")**3.0
  

CellsPerBin
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Ones(field, data):
      return np.ones(data.ActiveDimensions, dtype='float64')
  

**Convert Function Source**

No source available.

ChandraEmissivity
+++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _ChandraEmissivity(field, data):
      logT0 = np.log10(data["Temperature"]) - 7
      return ((data["NumberDensity"].astype('float64')**2.0) \
              *(10**(-0.0103*logT0**8 \
                     +0.0417*logT0**7 \
                     -0.0636*logT0**6 \
                     +0.1149*logT0**5 \
                     -0.3151*logT0**4 \
                     +0.6655*logT0**3 \
                     -1.1256*logT0**2 \
                     +1.0026*logT0**1 \
                     -0.6984*logT0) \
                +data["Metallicity"]*10**(0.0305*logT0**11 \
                                          -0.0045*logT0**10 \
                                          -0.3620*logT0**9 \
                                          +0.0513*logT0**8 \
                                          +1.6669*logT0**7 \
                                          -0.3854*logT0**6 \
                                          -3.3604*logT0**5 \
                                          +0.4728*logT0**4 \
                                          +4.5774*logT0**3 \
                                          -2.3661*logT0**2 \
                                          -1.6667*logT0**1 \
                                          -0.2193*logT0)))
  

**Convert Function Source**

.. code-block:: python

  def _convertChandraEmissivity(data):
      return 1.0 #1.0e-23*0.76**2
  

ComovingDensity
+++++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^3`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _ComovingDensity(field, data):
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data["density"]/ef
  

**Convert Function Source**

No source available.

Contours
++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Contours(field, data):
      return -np.ones_like(data["Ones"])
  

**Convert Function Source**

No source available.

CourantTimeStep
+++++++++++++++

   * Units: :math:`$\rm{s}$`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CourantTimeStep(field, data):
      t1 = data['dx'] / (
          data["SoundSpeed"] + \
          abs(data["x-velocity"]))
      t2 = data['dy'] / (
          data["SoundSpeed"] + \
          abs(data["y-velocity"]))
      t3 = data['dz'] / (
          data["SoundSpeed"] + \
          abs(data["z-velocity"]))
      return np.minimum(np.minimum(t1,t2),t3)
  

**Convert Function Source**

.. code-block:: python

  def _convertCourantTimeStep(data):
      # SoundSpeed and z-velocity are in cm/s, dx is in code
      return data.convert("cm")
  

CuttingPlaneBx
++++++++++++++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CuttingPlaneBx(field, data):
      x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                             for ax in 'xyz']
      b_vec = np.array([data["B%s" % ax] for ax in 'xyz'])
      return np.dot(x_vec, b_vec)
  

**Convert Function Source**

No source available.

CuttingPlaneBy
++++++++++++++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CuttingPlaneBy(field, data):
      x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                             for ax in 'xyz']
      b_vec = np.array([data["B%s" % ax] for ax in 'xyz'])
      return np.dot(y_vec, b_vec)
  

**Convert Function Source**

No source available.

CuttingPlaneVelocityX
+++++++++++++++++++++

   * Units: :math:`\rm{km}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CuttingPlaneVelocityX(field, data):
      x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                             for ax in 'xyz']
      bulk_velocity = data.get_field_parameter("bulk_velocity")
      if bulk_velocity == None:
          bulk_velocity = np.zeros(3)
      v_vec = np.array([data["%s-velocity" % ax] for ax in 'xyz']) \
                  - bulk_velocity[...,np.newaxis]
      return np.dot(x_vec, v_vec)
  

**Convert Function Source**

No source available.

CuttingPlaneVelocityY
+++++++++++++++++++++

   * Units: :math:`\rm{km}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _CuttingPlaneVelocityY(field, data):
      x_vec, y_vec, z_vec = [data.get_field_parameter("cp_%s_vec" % (ax))
                             for ax in 'xyz']
      bulk_velocity = data.get_field_parameter("bulk_velocity")
      if bulk_velocity == None:
          bulk_velocity = np.zeros(3)
      v_vec = np.array([data["%s-velocity" % ax] for ax in 'xyz']) \
                  - bulk_velocity[...,np.newaxis]
      return np.dot(y_vec, v_vec)
  

**Convert Function Source**

No source available.

DensityPerturbation
+++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _DensityPerturbation(field, data):
      rho_bar = rho_crit_now * data.pf.omega_matter * \
          data.pf.hubble_constant**2 * \
          (1.0 + data.pf.current_redshift)**3
      return ((data['Matter_Density'] - rho_bar) / rho_bar)
  

**Convert Function Source**

No source available.

DiskAngle
+++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _DiskAngle(field, data):
      return data['sph_theta']
  

**Convert Function Source**

No source available.

DivV
++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _DivV(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      ds = div_fac * data['dx'].flat[0]
      f  = data["x-velocity"][sl_right,1:-1,1:-1]/ds
      f -= data["x-velocity"][sl_left ,1:-1,1:-1]/ds
      if data.pf.dimensionality > 1:
          ds = div_fac * data['dy'].flat[0]
          f += data["y-velocity"][1:-1,sl_right,1:-1]/ds
          f -= data["y-velocity"][1:-1,sl_left ,1:-1]/ds
      if data.pf.dimensionality > 2:
          ds = div_fac * data['dz'].flat[0]
          f += data["z-velocity"][1:-1,1:-1,sl_right]/ds
          f -= data["z-velocity"][1:-1,1:-1,sl_left ]/ds
      new_field = np.zeros(data["x-velocity"].shape, dtype='float64')
      new_field[1:-1,1:-1,1:-1] = f
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertDivV(data):
      return data.convert("cm")**-1.0
  

DynamicalTime
+++++++++++++

   * Units: :math:`\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _DynamicalTime(field, data):
      """
      The formulation for the dynamical time is:
      M{sqrt(3pi/(16*G*rho))} or M{sqrt(3pi/(16G))*rho^-(1/2)}
      Note that we return in our natural units already
      """
      return (3.0*np.pi/(16*G*data["density"]))**(1./2.)
  

**Convert Function Source**

No source available.

Entropy
+++++++

   * Units: :math:`\rm{ergs}\ \rm{cm}^{3\gamma-3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Entropy(field, data):
      if data.has_field_parameter("mu"):
          mw = mh*data.get_field_parameter("mu")
      else :
          mw = mh
      try:
          gammam1 = data.pf["Gamma"] - 1.0
      except:
          gammam1 = 5./3. - 1.0
      return kboltz * data["Temperature"] / \
             ((data["density"]/mw)**gammam1)
  

**Convert Function Source**

No source available.

GridIndices
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _GridIndices(field, data):
      return np.ones(data["Ones"].shape)*(data.id-data._id_offset)
  

**Convert Function Source**

No source available.

GridLevel
+++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _GridLevel(field, data):
      return np.ones(data.ActiveDimensions)*(data.Level)
  

**Convert Function Source**

No source available.

Height
++++++

   * Units: :math:`cm`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Height(field, data):
      return data['cyl_z']
  

**Convert Function Source**

.. code-block:: python

  def _convertHeight(data):
      return data.convert("cm")
  

HeightAU
++++++++

   * Units: :math:`AU`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Height(field, data):
      return data['cyl_z']
  

**Convert Function Source**

.. code-block:: python

  def _convertHeightAU(data):
      return data.convert("au")
  

JeansMassMsun
+++++++++++++

   * Units: :math:`\rm{M_{\odot}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _JeansMassMsun(field,data):
      MJ_constant = (((5*kboltz)/(G*mh))**(1.5)) * \
      (3/(4*3.1415926535897931))**(0.5) / 1.989e33
  
      return (MJ_constant *
              ((data["Temperature"]/data["MeanMolecularWeight"])**(1.5)) *
              (data["density"]**(-0.5)))
  

**Convert Function Source**

No source available.

MachNumber
++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _MachNumber(field, data):
      """M{|v|/t_sound}"""
      return data["VelocityMagnitude"] / data["SoundSpeed"]
  

**Convert Function Source**

No source available.

MagneticEnergy
++++++++++++++

   * Units: :math:`\rm{ergs}\/\rm{cm}^{-3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _MagneticEnergy(field,data):
      """This assumes that your front end has provided Bx, By, Bz in
      units of Gauss. If you use MKS, make sure to write your own
      MagneticEnergy field to deal with non-unitary \mu_0.
      """
      return (data["Bx"]**2 + data["By"]**2 + data["Bz"]**2)/(8*np.pi)
  

**Convert Function Source**

No source available.

MagneticPressure
++++++++++++++++

   * Units: :math:`\rm{ergs}\/\rm{cm}^{-3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _MagneticPressure(field,data):
      return data['MagneticEnergy']
  

**Convert Function Source**

No source available.

Matter_Density
++++++++++++++

   * Units: :math:`\rm{g}/\rm{cm^3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Matter_Density(field,data):
      return (data['density'] + data['particle_density'])
  

**Convert Function Source**

No source available.

MeanMolecularWeight
+++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _MeanMolecularWeight(field,data):
      return (data["density"] / (mh *data["NumberDensity"]))
  

**Convert Function Source**

No source available.

Ones
++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Ones(field, data):
      return np.ones(data.ActiveDimensions, dtype='float64')
  

**Convert Function Source**

No source available.

OnesOverDx
++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _OnesOverDx(field, data):
      return np.ones(data["Ones"].shape,
                     dtype=data["density"].dtype)/data['dx']
  

**Convert Function Source**

No source available.

Overdensity
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Matter_Density(field,data):
      return (data['density'] + data['particle_density'])
  

**Convert Function Source**

.. code-block:: python

  def _Convert_Overdensity(data):
      return 1.0 / (rho_crit_now * data.pf.hubble_constant**2 * 
                  (1+data.pf.current_redshift)**3)
  

ParticleAngularMomentumX
++++++++++++++++++++++++

   * Units: :math:`\rm{g}\/\rm{cm}^2/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleAngularMomentumX(field, data):
      return data["CellMass"] * data["ParticleSpecificAngularMomentumX"]
  

**Convert Function Source**

No source available.

ParticleAngularMomentumY
++++++++++++++++++++++++

   * Units: :math:`\rm{g}\/\rm{cm}^2/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleAngularMomentumY(field, data):
      return data["CellMass"] * data["ParticleSpecificAngularMomentumY"]
  

**Convert Function Source**

No source available.

ParticleAngularMomentumZ
++++++++++++++++++++++++

   * Units: :math:`\rm{g}\/\rm{cm}^2/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleAngularMomentumZ(field, data):
      return data["CellMass"] * data["ParticleSpecificAngularMomentumZ"]
  

**Convert Function Source**

No source available.

ParticleMass
++++++++++++

   * Units: :math:`UNDEFINED`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def NullFunc(field, data):
      return
  

**Convert Function Source**

No source available.

ParticleRadius
++++++++++++++

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleRadius(field, data):
      return get_radius(data, "particle_position_")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiusCGS(data):
      return data.convert("cm")
  

ParticleRadiusAU
++++++++++++++++

   * Units: :math:`\rm{AU}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleRadius(field, data):
      return get_radius(data, "particle_position_")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiusAU(data):
      return data.convert("au")
  

ParticleRadiusCode
++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleRadius(field, data):
      return get_radius(data, "particle_position_")
  

**Convert Function Source**

No source available.

ParticleRadiusMpc
+++++++++++++++++

   * Units: :math:`\rm{Mpc}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleRadius(field, data):
      return get_radius(data, "particle_position_")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiusMpc(data):
      return data.convert("mpc")
  

ParticleRadiuskpc
+++++++++++++++++

   * Units: :math:`\rm{kpc}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleRadius(field, data):
      return get_radius(data, "particle_position_")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiuskpc(data):
      return data.convert("kpc")
  

ParticleRadiuskpch
++++++++++++++++++

   * Units: :math:`\rm{kpc}/\rm{h}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleRadius(field, data):
      return get_radius(data, "particle_position_")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiuskpch(data):
      return data.convert("kpch")
  

ParticleRadiuspc
++++++++++++++++

   * Units: :math:`\rm{pc}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleRadius(field, data):
      return get_radius(data, "particle_position_")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiuspc(data):
      return data.convert("pc")
  

ParticleSpecificAngularMomentumX
++++++++++++++++++++++++++++++++

   * Units: :math:`\rm{cm}^2/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleSpecificAngularMomentumX(field, data):
      if data.has_field_parameter("bulk_velocity"):
          bv = data.get_field_parameter("bulk_velocity")
      else: bv = np.zeros(3, dtype='float64')
      center = data.get_field_parameter('center')
      y = data["particle_position_y"] - center[1]
      z = data["particle_position_z"] - center[2]
      yv = data["particle_velocity_y"] - bv[1]
      zv = data["particle_velocity_z"] - bv[2]
      return yv*z - zv*y
  

**Convert Function Source**

.. code-block:: python

  def _convertSpecificAngularMomentum(data):
      return data.convert("cm")
  

ParticleSpecificAngularMomentumXKMSMPC
++++++++++++++++++++++++++++++++++++++

   * Units: :math:`\rm{cm}^2/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleSpecificAngularMomentumX(field, data):
      if data.has_field_parameter("bulk_velocity"):
          bv = data.get_field_parameter("bulk_velocity")
      else: bv = np.zeros(3, dtype='float64')
      center = data.get_field_parameter('center')
      y = data["particle_position_y"] - center[1]
      z = data["particle_position_z"] - center[2]
      yv = data["particle_velocity_y"] - bv[1]
      zv = data["particle_velocity_z"] - bv[2]
      return yv*z - zv*y
  

**Convert Function Source**

.. code-block:: python

  def _convertSpecificAngularMomentumKMSMPC(data):
      return km_per_cm*data.convert("mpc")
  

ParticleSpecificAngularMomentumY
++++++++++++++++++++++++++++++++

   * Units: :math:`\rm{cm}^2/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleSpecificAngularMomentumY(field, data):
      if data.has_field_parameter("bulk_velocity"):
          bv = data.get_field_parameter("bulk_velocity")
      else: bv = np.zeros(3, dtype='float64')
      center = data.get_field_parameter('center')
      x = data["particle_position_x"] - center[0]
      z = data["particle_position_z"] - center[2]
      xv = data["particle_velocity_x"] - bv[0]
      zv = data["particle_velocity_z"] - bv[2]
      return -(xv*z - zv*x)
  

**Convert Function Source**

.. code-block:: python

  def _convertSpecificAngularMomentum(data):
      return data.convert("cm")
  

ParticleSpecificAngularMomentumYKMSMPC
++++++++++++++++++++++++++++++++++++++

   * Units: :math:`\rm{cm}^2/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleSpecificAngularMomentumY(field, data):
      if data.has_field_parameter("bulk_velocity"):
          bv = data.get_field_parameter("bulk_velocity")
      else: bv = np.zeros(3, dtype='float64')
      center = data.get_field_parameter('center')
      x = data["particle_position_x"] - center[0]
      z = data["particle_position_z"] - center[2]
      xv = data["particle_velocity_x"] - bv[0]
      zv = data["particle_velocity_z"] - bv[2]
      return -(xv*z - zv*x)
  

**Convert Function Source**

.. code-block:: python

  def _convertSpecificAngularMomentumKMSMPC(data):
      return km_per_cm*data.convert("mpc")
  

ParticleSpecificAngularMomentumZ
++++++++++++++++++++++++++++++++

   * Units: :math:`\rm{cm}^2/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleSpecificAngularMomentumZ(field, data):
      if data.has_field_parameter("bulk_velocity"):
          bv = data.get_field_parameter("bulk_velocity")
      else: bv = np.zeros(3, dtype='float64')
      center = data.get_field_parameter('center')
      x = data["particle_position_x"] - center[0]
      y = data["particle_position_y"] - center[1]
      xv = data["particle_velocity_x"] - bv[0]
      yv = data["particle_velocity_y"] - bv[1]
      return xv*y - yv*x
  

**Convert Function Source**

.. code-block:: python

  def _convertSpecificAngularMomentum(data):
      return data.convert("cm")
  

ParticleSpecificAngularMomentumZKMSMPC
++++++++++++++++++++++++++++++++++++++

   * Units: :math:`\rm{cm}^2/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleSpecificAngularMomentumZ(field, data):
      if data.has_field_parameter("bulk_velocity"):
          bv = data.get_field_parameter("bulk_velocity")
      else: bv = np.zeros(3, dtype='float64')
      center = data.get_field_parameter('center')
      x = data["particle_position_x"] - center[0]
      y = data["particle_position_y"] - center[1]
      xv = data["particle_velocity_x"] - bv[0]
      yv = data["particle_velocity_y"] - bv[1]
      return xv*y - yv*x
  

**Convert Function Source**

.. code-block:: python

  def _convertSpecificAngularMomentumKMSMPC(data):
      return km_per_cm*data.convert("mpc")
  

ParticleVelocityMagnitude
+++++++++++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleVelocityMagnitude(field, data):
      """M{|v|}"""
      bulk_velocity = data.get_field_parameter("bulk_velocity")
      if bulk_velocity == None:
          bulk_velocity = np.zeros(3)
      return ( (data["particle_velocity_x"]-bulk_velocity[0])**2.0 + \
               (data["particle_velocity_y"]-bulk_velocity[1])**2.0 + \
               (data["particle_velocity_z"]-bulk_velocity[2])**2.0 )**(1.0/2.0)
  

**Convert Function Source**

No source available.

PlasmaBeta
++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _PlasmaBeta(field,data):
      """This assumes that your front end has provided Bx, By, Bz in
      units of Gauss. If you use MKS, make sure to write your own
      PlasmaBeta field to deal with non-unitary \mu_0.
      """
      return data['Pressure']/data['MagneticEnergy']
  

**Convert Function Source**

No source available.

Pressure
++++++++

   * Units: :math:`\rm{dyne}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Pressure(field, data):
      """M{(Gamma-1.0)*rho*E}"""
      return (data.pf["Gamma"] - 1.0) * \
             data["density"] * data["ThermalEnergy"]
  

**Convert Function Source**

No source available.

RadialMachNumber
++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _RadialMachNumber(field, data):
      """M{|v|/t_sound}"""
      return np.abs(data["RadialVelocity"]) / data["SoundSpeed"]
  

**Convert Function Source**

No source available.

RadialVelocity
++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _RadialVelocity(field, data):
      normal = data.get_field_parameter("normal")
      velocities = obtain_rv_vec(data)    
      theta = data['sph_theta']
      phi   = data['sph_phi']
  
      return get_sph_r_component(velocities, theta, phi, normal)
  

**Convert Function Source**

No source available.

RadialVelocityABS
+++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _RadialVelocityABS(field, data):
      return np.abs(_RadialVelocity(field, data))
  

**Convert Function Source**

No source available.

RadialVelocityKMS
+++++++++++++++++

   * Units: :math:`\rm{km}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _RadialVelocity(field, data):
      normal = data.get_field_parameter("normal")
      velocities = obtain_rv_vec(data)    
      theta = data['sph_theta']
      phi   = data['sph_phi']
  
      return get_sph_r_component(velocities, theta, phi, normal)
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadialVelocityKMS(data):
      return km_per_cm
  

RadialVelocityKMSABS
++++++++++++++++++++

   * Units: :math:`\rm{km}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _RadialVelocityABS(field, data):
      return np.abs(_RadialVelocity(field, data))
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadialVelocityKMS(data):
      return km_per_cm
  

Radius
++++++

   * Units: :math:`\rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Radius(field, data):
      return get_radius(data, "")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiusCGS(data):
      return data.convert("cm")
  

RadiusAU
++++++++

   * Units: :math:`\rm{AU}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Radius(field, data):
      return get_radius(data, "")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiusAU(data):
      return data.convert("au")
  

RadiusCode
++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Radius(field, data):
      return get_radius(data, "")
  

**Convert Function Source**

No source available.

RadiusMpc
+++++++++

   * Units: :math:`\rm{Mpc}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Radius(field, data):
      return get_radius(data, "")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiusMpc(data):
      return data.convert("mpc")
  

Radiuskpc
+++++++++

   * Units: :math:`\rm{kpc}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Radius(field, data):
      return get_radius(data, "")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiuskpc(data):
      return data.convert("kpc")
  

Radiuspc
++++++++

   * Units: :math:`\rm{pc}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Radius(field, data):
      return get_radius(data, "")
  

**Convert Function Source**

.. code-block:: python

  def _ConvertRadiuspc(data):
      return data.convert("pc")
  

SZKinetic
+++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SZKinetic(field, data):
      vel_axis = data.get_field_parameter('axis')
      if vel_axis > 2:
          raise NeedsParameter(['axis'])
      vel = data["%s-velocity" % ({0:'x',1:'y',2:'z'}[vel_axis])]
      return (vel*data["density"])
  

**Convert Function Source**

.. code-block:: python

  def _convertSZKinetic(data):
      return 0.88*((sigma_thompson/mh)/clight)
  

SZY
+++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SZY(field, data):
      return (data["density"]*data["Temperature"])
  

**Convert Function Source**

.. code-block:: python

  def _convertSZY(data):
      conv = (0.88/mh) * (kboltz)/(me * clight*clight) * sigma_thompson
      return conv
  

Shear
+++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Shear(field, data):
      """
      Shear is defined as [(dvx/dy + dvy/dx)^2 + (dvz/dy + dvy/dz)^2 +
                           (dvx/dz + dvz/dx)^2 ]^(0.5)
      where dvx/dy = [vx(j-1) - vx(j+1)]/[2dy]
      and is in units of s^(-1)
      (it's just like vorticity except add the derivative pairs instead
       of subtracting them)
      """
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["x-velocity"].shape)
      if data.pf.dimensionality > 1:
          dvydx = (data["y-velocity"][sl_right,1:-1,1:-1] -
                  data["y-velocity"][sl_left,1:-1,1:-1]) \
                  / (div_fac*data["dx"].flat[0])
          dvxdy = (data["x-velocity"][1:-1,sl_right,1:-1] -
                  data["x-velocity"][1:-1,sl_left,1:-1]) \
                  / (div_fac*data["dy"].flat[0])
          new_field[1:-1,1:-1,1:-1] += (dvydx + dvxdy)**2.0
          del dvydx, dvxdy
      if data.pf.dimensionality > 2:
          dvzdy = (data["z-velocity"][1:-1,sl_right,1:-1] -
                  data["z-velocity"][1:-1,sl_left,1:-1]) \
                  / (div_fac*data["dy"].flat[0])
          dvydz = (data["y-velocity"][1:-1,1:-1,sl_right] -
                  data["y-velocity"][1:-1,1:-1,sl_left]) \
                  / (div_fac*data["dz"].flat[0])
          new_field[1:-1,1:-1,1:-1] += (dvzdy + dvydz)**2.0
          del dvzdy, dvydz
          dvxdz = (data["x-velocity"][1:-1,1:-1,sl_right] -
                  data["x-velocity"][1:-1,1:-1,sl_left]) \
                  / (div_fac*data["dz"].flat[0])
          dvzdx = (data["z-velocity"][sl_right,1:-1,1:-1] -
                  data["z-velocity"][sl_left,1:-1,1:-1]) \
                  / (div_fac*data["dx"].flat[0])
          new_field[1:-1,1:-1,1:-1] += (dvxdz + dvzdx)**2.0
          del dvxdz, dvzdx
      new_field = new_field**0.5
      new_field = np.abs(new_field)
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertShear(data):
      return data.convert("cm")**-1.0
  

ShearCriterion
++++++++++++++

   * Units: :math:`\rm{cm}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _ShearCriterion(field, data):
      """
      Shear is defined as [(dvx/dy + dvy/dx)^2 + (dvz/dy + dvy/dz)^2 +
                           (dvx/dz + dvz/dx)^2 ]^(0.5)
      where dvx/dy = [vx(j-1) - vx(j+1)]/[2dy]
      and is in units of s^(-1)
      (it's just like vorticity except add the derivative pairs instead
       of subtracting them)
  
      Divide by c_s to leave Shear in units of cm**-1, which 
      can be compared against the inverse of the local cell size (1/dx) 
      to determine if refinement should occur.
      """
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["x-velocity"].shape)
      if data.pf.dimensionality > 1:
          dvydx = (data["y-velocity"][sl_right,1:-1,1:-1] -
                  data["y-velocity"][sl_left,1:-1,1:-1]) \
                  / (div_fac*data["dx"].flat[0])
          dvxdy = (data["x-velocity"][1:-1,sl_right,1:-1] -
                  data["x-velocity"][1:-1,sl_left,1:-1]) \
                  / (div_fac*data["dy"].flat[0])
          new_field[1:-1,1:-1,1:-1] += (dvydx + dvxdy)**2.0
          del dvydx, dvxdy
      if data.pf.dimensionality > 2:
          dvzdy = (data["z-velocity"][1:-1,sl_right,1:-1] -
                  data["z-velocity"][1:-1,sl_left,1:-1]) \
                  / (div_fac*data["dy"].flat[0])
          dvydz = (data["y-velocity"][1:-1,1:-1,sl_right] -
                  data["y-velocity"][1:-1,1:-1,sl_left]) \
                  / (div_fac*data["dz"].flat[0])
          new_field[1:-1,1:-1,1:-1] += (dvzdy + dvydz)**2.0
          del dvzdy, dvydz
          dvxdz = (data["x-velocity"][1:-1,1:-1,sl_right] -
                  data["x-velocity"][1:-1,1:-1,sl_left]) \
                  / (div_fac*data["dz"].flat[0])
          dvzdx = (data["z-velocity"][sl_right,1:-1,1:-1] -
                  data["z-velocity"][sl_left,1:-1,1:-1]) \
                  / (div_fac*data["dx"].flat[0])
          new_field[1:-1,1:-1,1:-1] += (dvxdz + dvzdx)**2.0
          del dvxdz, dvzdx
      new_field /= data["SoundSpeed"]**2.0
      new_field = new_field**(0.5)
      new_field = np.abs(new_field)
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertShearCriterion(data):
      return data.convert("cm")**-1.0
  

ShearMach
+++++++++

   * Units: :math:`\rm{Mach}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _ShearMach(field, data):
      """
      Dimensionless Shear (ShearMach) is defined nearly the same as shear, 
      except that it is scaled by the local dx/dy/dz and the local sound speed.
      So it results in a unitless quantity that is effectively measuring 
      shear in mach number.  
  
      In order to avoid discontinuities created by multiplying by dx/dy/dz at
      grid refinement boundaries, we also multiply by 2**GridLevel.
  
      Shear (Mach) = [(dvx + dvy)^2 + (dvz + dvy)^2 +
                      (dvx + dvz)^2  ]^(0.5) / c_sound
      """
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["x-velocity"].shape)
      if data.pf.dimensionality > 1:
          dvydx = (data["y-velocity"][sl_right,1:-1,1:-1] -
                  data["y-velocity"][sl_left,1:-1,1:-1]) \
                  / (div_fac)
          dvxdy = (data["x-velocity"][1:-1,sl_right,1:-1] -
                  data["x-velocity"][1:-1,sl_left,1:-1]) \
                  / (div_fac)
          new_field[1:-1,1:-1,1:-1] += (dvydx + dvxdy)**2.0
          del dvydx, dvxdy
      if data.pf.dimensionality > 2:
          dvzdy = (data["z-velocity"][1:-1,sl_right,1:-1] -
                  data["z-velocity"][1:-1,sl_left,1:-1]) \
                  / (div_fac)
          dvydz = (data["y-velocity"][1:-1,1:-1,sl_right] -
                  data["y-velocity"][1:-1,1:-1,sl_left]) \
                  / (div_fac)
          new_field[1:-1,1:-1,1:-1] += (dvzdy + dvydz)**2.0
          del dvzdy, dvydz
          dvxdz = (data["x-velocity"][1:-1,1:-1,sl_right] -
                  data["x-velocity"][1:-1,1:-1,sl_left]) \
                  / (div_fac)
          dvzdx = (data["z-velocity"][sl_right,1:-1,1:-1] -
                  data["z-velocity"][sl_left,1:-1,1:-1]) \
                  / (div_fac)
          new_field[1:-1,1:-1,1:-1] += (dvxdz + dvzdx)**2.0
          del dvxdz, dvzdx
      new_field *= ((2.0**data.level)/data["SoundSpeed"])**2.0
      new_field = new_field**0.5
      new_field = np.abs(new_field)
      return new_field
  

**Convert Function Source**

No source available.

SoundSpeed
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SoundSpeed(field, data):
      if data.pf["EOSType"] == 1:
          return np.ones(data["density"].shape, dtype='float64') * \
                  data.pf["EOSSoundSpeed"]
      return ( data.pf["Gamma"]*data["Pressure"] / \
               data["density"] )**(1.0/2.0)
  

**Convert Function Source**

No source available.

SpecificAngularMomentumX
++++++++++++++++++++++++

   * Units: :math:`\rm{cm}^2/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpecificAngularMomentumX(field, data):
      xv, yv, zv = obtain_velocities(data)
      rv = obtain_rvec(data)
      return yv*rv[2,:] - zv*rv[1,:]
  

**Convert Function Source**

.. code-block:: python

  def _convertSpecificAngularMomentum(data):
      return data.convert("cm")
  

SpecificAngularMomentumY
++++++++++++++++++++++++

   * Units: :math:`\rm{cm}^2/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpecificAngularMomentumY(field, data):
      xv, yv, zv = obtain_velocities(data)
      rv = obtain_rvec(data)
      return -(xv*rv[2,:] - zv*rv[0,:])
  

**Convert Function Source**

.. code-block:: python

  def _convertSpecificAngularMomentum(data):
      return data.convert("cm")
  

SpecificAngularMomentumZ
++++++++++++++++++++++++

   * Units: :math:`\rm{cm}^2/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpecificAngularMomentumZ(field, data):
      xv, yv, zv = obtain_velocities(data)
      rv = obtain_rvec(data)
      return xv*rv[1,:] - yv*rv[0,:]
  

**Convert Function Source**

.. code-block:: python

  def _convertSpecificAngularMomentum(data):
      return data.convert("cm")
  

StarMassMsun
++++++++++++

   * Units: :math:`M_{\odot}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _StarMass(field,data):
      return data["star_density"] * data["CellVolume"]
  

**Convert Function Source**

.. code-block:: python

  def _convertCellMassMsun(data):
      return 5.027854e-34 # g^-1
  

TangentialOverVelocityMagnitude
+++++++++++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _TangentialOverVelocityMagnitude(field, data):
      return np.abs(data["TangentialVelocity"])/np.abs(data["VelocityMagnitude"])
  

**Convert Function Source**

No source available.

TangentialVelocity
++++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _TangentialVelocity(field, data):
      return np.sqrt(data["VelocityMagnitude"]**2.0
                   - data["RadialVelocity"]**2.0)
  

**Convert Function Source**

No source available.

TempkeV
+++++++

   * Units: :math:`\rm{keV}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _TempkeV(field, data):
      return data["Temperature"] * keV_per_K
  

**Convert Function Source**

No source available.

TotalMass
+++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _TotalMass(field,data):
      return (data["density"]+data["particle_density"]) * data["CellVolume"]
  

**Convert Function Source**

No source available.

TotalMassMsun
+++++++++++++

   * Units: :math:`M_{\odot}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _TotalMass(field,data):
      return (data["density"]+data["particle_density"]) * data["CellVolume"]
  

**Convert Function Source**

.. code-block:: python

  def _convertCellMassMsun(data):
      return 5.027854e-34 # g^-1
  

VelocityMagnitude
+++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VelocityMagnitude(field, data):
      """M{|v|}"""
      velocities = obtain_rv_vec(data)
      return np.sqrt(np.sum(velocities**2,axis=0))
  

**Convert Function Source**

No source available.

VorticityGrowthMagnitude
++++++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityGrowthMagnitude(field, data):
      result = np.sqrt(data["VorticityGrowthX"]**2 +
                       data["VorticityGrowthY"]**2 +
                       data["VorticityGrowthZ"]**2)
      dot = np.zeros(result.shape)
      for ax in "XYZ":
          dot += data["Vorticity%s" % ax] * data["VorticityGrowth%s" % ax]
      result = np.sign(dot) * result
      return result
  

**Convert Function Source**

No source available.

VorticityGrowthMagnitudeABS
+++++++++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityGrowthMagnitudeABS(field, data):
      return np.sqrt(data["VorticityGrowthX"]**2 +
                     data["VorticityGrowthY"]**2 +
                     data["VorticityGrowthZ"]**2)
  

**Convert Function Source**

No source available.

VorticityGrowthTimescale
++++++++++++++++++++++++

   * Units: :math:`\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityGrowthTimescale(field, data):
      domegax_dt = data["VorticityX"] / data["VorticityGrowthX"]
      domegay_dt = data["VorticityY"] / data["VorticityGrowthY"]
      domegaz_dt = data["VorticityZ"] / data["VorticityGrowthZ"]
      return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
  

**Convert Function Source**

No source available.

VorticityGrowthX
++++++++++++++++

   * Units: :math:`\rm{s}^{-2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityGrowthX(field, data):
      return -data["VorticityStretchingX"] - data["BaroclinicVorticityX"]
  

**Convert Function Source**

No source available.

VorticityGrowthY
++++++++++++++++

   * Units: :math:`\rm{s}^{-2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityGrowthY(field, data):
      return -data["VorticityStretchingY"] - data["BaroclinicVorticityY"]
  

**Convert Function Source**

No source available.

VorticityGrowthZ
++++++++++++++++

   * Units: :math:`\rm{s}^{-2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityGrowthZ(field, data):
      return -data["VorticityStretchingZ"] - data["BaroclinicVorticityZ"]
  

**Convert Function Source**

No source available.

VorticityMagnitude
++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityMagnitude(field, data):
      return np.sqrt(data["VorticityX"]**2 +
                     data["VorticityY"]**2 +
                     data["VorticityZ"]**2)
  

**Convert Function Source**

No source available.

VorticityRPGrowthMagnitude
++++++++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityGrowthMagnitude(field, data):
      result = np.sqrt(data["VorticityGrowthX"]**2 +
                       data["VorticityGrowthY"]**2 +
                       data["VorticityGrowthZ"]**2)
      dot = np.zeros(result.shape)
      for ax in "XYZ":
          dot += data["Vorticity%s" % ax] * data["VorticityGrowth%s" % ax]
      result = np.sign(dot) * result
      return result
  

**Convert Function Source**

No source available.

VorticityRPGrowthMagnitudeABS
+++++++++++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityRPGrowthMagnitudeABS(field, data):
      return np.sqrt(data["VorticityRPGrowthX"]**2 +
                     data["VorticityRPGrowthY"]**2 +
                     data["VorticityRPGrowthZ"]**2)
  

**Convert Function Source**

No source available.

VorticityRPGrowthTimescale
++++++++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityRPGrowthTimescale(field, data):
      domegax_dt = data["VorticityX"] / data["VorticityRPGrowthX"]
      domegay_dt = data["VorticityY"] / data["VorticityRPGrowthY"]
      domegaz_dt = data["VorticityZ"] / data["VorticityRPGrowthZ"]
      return np.sqrt(domegax_dt**2 + domegay_dt**2 + domegaz_dt**2)
  

**Convert Function Source**

No source available.

VorticityRPGrowthX
++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityRPGrowthX(field, data):
      return -data["VorticityStretchingX"] - data["BaroclinicVorticityX"] \
             -data["VorticityRadPressureX"]
  

**Convert Function Source**

No source available.

VorticityRPGrowthY
++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityRPGrowthY(field, data):
      return -data["VorticityStretchingY"] - data["BaroclinicVorticityY"] \
             -data["VorticityRadPressureY"]
  

**Convert Function Source**

No source available.

VorticityRPGrowthZ
++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityRPGrowthZ(field, data):
      return -data["VorticityStretchingZ"] - data["BaroclinicVorticityZ"] \
             -data["VorticityRadPressureZ"]
  

**Convert Function Source**

No source available.

VorticityRadPressureMagnitude
+++++++++++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityRadPressureMagnitude(field, data):
      return np.sqrt(data["VorticityRadPressureX"]**2 +
                     data["VorticityRadPressureY"]**2 +
                     data["VorticityRadPressureZ"]**2)
  

**Convert Function Source**

No source available.

VorticityRadPressureX
+++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityRadPressureX(field, data):
      rho = data["density"].astype('float64')
      return (data["RadAccel2"] * data["gradDensityZ"] -
              data["RadAccel3"] * data["gradDensityY"]) / rho
  

**Convert Function Source**

.. code-block:: python

  def _convertRadAccel(data):
      return data.convert("x-velocity")/data.convert("Time")
  

VorticityRadPressureY
+++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityRadPressureY(field, data):
      rho = data["density"].astype('float64')
      return (data["RadAccel3"] * data["gradDensityX"] -
              data["RadAccel1"] * data["gradDensityZ"]) / rho
  

**Convert Function Source**

.. code-block:: python

  def _convertRadAccel(data):
      return data.convert("x-velocity")/data.convert("Time")
  

VorticityRadPressureZ
+++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityRadPressureZ(field, data):
      rho = data["density"].astype('float64')
      return (data["RadAccel1"] * data["gradDensityY"] -
              data["RadAccel2"] * data["gradDensityX"]) / rho
  

**Convert Function Source**

.. code-block:: python

  def _convertRadAccel(data):
      return data.convert("x-velocity")/data.convert("Time")
  

VorticitySquared
++++++++++++++++

   * Units: :math:`\rm{s}^{-2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticitySquared(field, data):
      mylog.debug("Generating vorticity on %s", data)
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["x-velocity"].shape)
      dvzdy = (data["z-velocity"][1:-1,sl_right,1:-1] -
               data["z-velocity"][1:-1,sl_left,1:-1]) \
               / (div_fac*data["dy"].flat[0])
      dvydz = (data["y-velocity"][1:-1,1:-1,sl_right] -
               data["y-velocity"][1:-1,1:-1,sl_left]) \
               / (div_fac*data["dz"].flat[0])
      new_field[1:-1,1:-1,1:-1] += (dvzdy - dvydz)**2.0
      del dvzdy, dvydz
      dvxdz = (data["x-velocity"][1:-1,1:-1,sl_right] -
               data["x-velocity"][1:-1,1:-1,sl_left]) \
               / (div_fac*data["dz"].flat[0])
      dvzdx = (data["z-velocity"][sl_right,1:-1,1:-1] -
               data["z-velocity"][sl_left,1:-1,1:-1]) \
               / (div_fac*data["dx"].flat[0])
      new_field[1:-1,1:-1,1:-1] += (dvxdz - dvzdx)**2.0
      del dvxdz, dvzdx
      dvydx = (data["y-velocity"][sl_right,1:-1,1:-1] -
               data["y-velocity"][sl_left,1:-1,1:-1]) \
               / (div_fac*data["dx"].flat[0])
      dvxdy = (data["x-velocity"][1:-1,sl_right,1:-1] -
               data["x-velocity"][1:-1,sl_left,1:-1]) \
               / (div_fac*data["dy"].flat[0])
      new_field[1:-1,1:-1,1:-1] += (dvydx - dvxdy)**2.0
      del dvydx, dvxdy
      new_field = np.abs(new_field)
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertVorticitySquared(data):
      return data.convert("cm")**-2.0
  

VorticityStretchingMagnitude
++++++++++++++++++++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityStretchingMagnitude(field, data):
      return np.sqrt(data["VorticityStretchingX"]**2 +
                     data["VorticityStretchingY"]**2 +
                     data["VorticityStretchingZ"]**2)
  

**Convert Function Source**

No source available.

VorticityStretchingX
++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityStretchingX(field, data):
      return data["DivV"] * data["VorticityX"]
  

**Convert Function Source**

No source available.

VorticityStretchingY
++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityStretchingY(field, data):
      return data["DivV"] * data["VorticityY"]
  

**Convert Function Source**

No source available.

VorticityStretchingZ
++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityStretchingZ(field, data):
      return data["DivV"] * data["VorticityZ"]
  

**Convert Function Source**

No source available.

VorticityX
++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityX(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["z-velocity"].shape, dtype='float64')
      new_field[1:-1,1:-1,1:-1] = (data["z-velocity"][1:-1,sl_right,1:-1] -
                                   data["z-velocity"][1:-1,sl_left,1:-1]) \
                                   / (div_fac*data["dy"].flat[0])
      new_field[1:-1,1:-1,1:-1] -= (data["y-velocity"][1:-1,1:-1,sl_right] -
                                    data["y-velocity"][1:-1,1:-1,sl_left]) \
                                    / (div_fac*data["dz"].flat[0])
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertVorticity(data):
      return 1.0/data.convert("cm")
  

VorticityY
++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityY(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["z-velocity"].shape, dtype='float64')
      new_field[1:-1,1:-1,1:-1] = (data["x-velocity"][1:-1,1:-1,sl_right] -
                                   data["x-velocity"][1:-1,1:-1,sl_left]) \
                                   / (div_fac*data["dz"].flat[0])
      new_field[1:-1,1:-1,1:-1] -= (data["z-velocity"][sl_right,1:-1,1:-1] -
                                    data["z-velocity"][sl_left,1:-1,1:-1]) \
                                    / (div_fac*data["dx"].flat[0])
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertVorticity(data):
      return 1.0/data.convert("cm")
  

VorticityZ
++++++++++

   * Units: :math:`\rm{s}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _VorticityZ(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["x-velocity"].shape, dtype='float64')
      new_field[1:-1,1:-1,1:-1] = (data["y-velocity"][sl_right,1:-1,1:-1] -
                                   data["y-velocity"][sl_left,1:-1,1:-1]) \
                                   / (div_fac*data["dx"].flat[0])
      new_field[1:-1,1:-1,1:-1] -= (data["x-velocity"][1:-1,sl_right,1:-1] -
                                    data["x-velocity"][1:-1,sl_left,1:-1]) \
                                    / (div_fac*data["dy"].flat[0])
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertVorticity(data):
      return 1.0/data.convert("cm")
  

WeakLensingConvergence
++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _DensityPerturbation(field, data):
      rho_bar = rho_crit_now * data.pf.omega_matter * \
          data.pf.hubble_constant**2 * \
          (1.0 + data.pf.current_redshift)**3
      return ((data['Matter_Density'] - rho_bar) / rho_bar)
  

**Convert Function Source**

.. code-block:: python

  def _convertConvergence(data):
      if not data.pf.parameters.has_key('cosmology_calculator'):
          data.pf.parameters['cosmology_calculator'] = Cosmology(
              HubbleConstantNow=(100.*data.pf.hubble_constant),
              OmegaMatterNow=data.pf.omega_matter, OmegaLambdaNow=data.pf.omega_lambda)
      # observer to lens
      DL = data.pf.parameters['cosmology_calculator'].AngularDiameterDistance(
          data.pf.parameters['observer_redshift'], data.pf.current_redshift)
      # observer to source
      DS = data.pf.parameters['cosmology_calculator'].AngularDiameterDistance(
          data.pf.parameters['observer_redshift'], data.pf.parameters['lensing_source_redshift'])
      # lens to source
      DLS = data.pf.parameters['cosmology_calculator'].AngularDiameterDistance(
          data.pf.current_redshift, data.pf.parameters['lensing_source_redshift'])
      return (((DL * DLS) / DS) * (1.5e14 * data.pf.omega_matter * 
                                  (data.pf.hubble_constant / speed_of_light_cgs)**2 *
                                  (1 + data.pf.current_redshift)))
  

XRayEmissivity
++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _XRayEmissivity(field, data):
      return ((data["density"].astype('float64')**2.0) \
              *data["Temperature"]**0.5)
  

**Convert Function Source**

.. code-block:: python

  def _convertXRayEmissivity(data):
      return 2.168e60
  

Zeros
+++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Zeros(field, data):
      return np.zeros(data.ActiveDimensions, dtype='float64')
  

**Convert Function Source**

No source available.

cyl_R
+++++

   * Units: :math:`\rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_R(field, data):
      center = data.get_field_parameter("center")
      normal = data.get_field_parameter("normal")
        
      coords = obtain_rvec(data)
  
      return get_cyl_r(coords, normal)
  

**Convert Function Source**

.. code-block:: python

  def _Convert_cyl_R_CGS(data):
     return data.convert("cm")
  

cyl_RCode
+++++++++

   * Units: :math:`Radius (code)`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_R(field, data):
      center = data.get_field_parameter("center")
      normal = data.get_field_parameter("normal")
        
      coords = obtain_rvec(data)
  
      return get_cyl_r(coords, normal)
  

**Convert Function Source**

No source available.

cyl_RadialVelocity
++++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_RadialVelocity(field, data):
      normal = data.get_field_parameter("normal")
      velocities = obtain_rv_vec(data)
  
      theta = data['cyl_theta']
  
      return get_cyl_r_component(velocities, theta, normal)
  

**Convert Function Source**

No source available.

cyl_RadialVelocityABS
+++++++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_RadialVelocityABS(field, data):
      return np.abs(_cyl_RadialVelocity(field, data))
  

**Convert Function Source**

No source available.

cyl_RadialVelocityKMS
+++++++++++++++++++++

   * Units: :math:`\rm{km}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_RadialVelocity(field, data):
      normal = data.get_field_parameter("normal")
      velocities = obtain_rv_vec(data)
  
      theta = data['cyl_theta']
  
      return get_cyl_r_component(velocities, theta, normal)
  

**Convert Function Source**

.. code-block:: python

  def _Convert_cyl_RadialVelocityKMS(data):
      return km_per_cm
  

cyl_RadialVelocityKMSABS
++++++++++++++++++++++++

   * Units: :math:`\rm{km}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_RadialVelocityABS(field, data):
      return np.abs(_cyl_RadialVelocity(field, data))
  

**Convert Function Source**

.. code-block:: python

  def _Convert_cyl_RadialVelocityKMS(data):
      return km_per_cm
  

cyl_TangentialVelocity
++++++++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_TangentialVelocity(field, data):
      normal = data.get_field_parameter("normal")
      velocities = obtain_rv_vec(data)
      theta = data['cyl_theta']
  
      return get_cyl_theta_component(velocities, theta, normal)
  

**Convert Function Source**

No source available.

cyl_TangentialVelocityABS
+++++++++++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_TangentialVelocityABS(field, data):
      return np.abs(_cyl_TangentialVelocity(field, data))
  

**Convert Function Source**

No source available.

cyl_TangentialVelocityKMS
+++++++++++++++++++++++++

   * Units: :math:`\rm{km}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_TangentialVelocity(field, data):
      normal = data.get_field_parameter("normal")
      velocities = obtain_rv_vec(data)
      theta = data['cyl_theta']
  
      return get_cyl_theta_component(velocities, theta, normal)
  

**Convert Function Source**

.. code-block:: python

  def _Convert_cyl_TangentialVelocityKMS(data):
      return km_per_cm
  

cyl_TangentialVelocityKMSABS
++++++++++++++++++++++++++++

   * Units: :math:`\rm{km}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_TangentialVelocityABS(field, data):
      return np.abs(_cyl_TangentialVelocity(field, data))
  

**Convert Function Source**

.. code-block:: python

  def _Convert_cyl_TangentialVelocityKMS(data):
      return km_per_cm
  

cyl_theta
+++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_theta(field, data):
      center = data.get_field_parameter("center")
      normal = data.get_field_parameter("normal")
      
      coords = obtain_rvec(data)
  
      return get_cyl_theta(coords, normal)
  

**Convert Function Source**

No source available.

cyl_z
+++++

   * Units: :math:`\rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cyl_z(field, data):
      center = data.get_field_parameter("center")
      normal = data.get_field_parameter("normal")
      
      coords = obtain_rvec(data)
  
      return get_cyl_z(coords, normal)
  

**Convert Function Source**

.. code-block:: python

  def _Convert_cyl_z_CGS(data):
     return data.convert("cm")
  

dx
++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _dx(field, data):
      return data.dds[0]
      return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[0]
  

**Convert Function Source**

No source available.

dy
++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _dy(field, data):
      return data.dds[1]
      return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[1]
  

**Convert Function Source**

No source available.

dz
++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _dz(field, data):
      return data.dds[2]
      return np.ones(data.ActiveDimensions, dtype='float64') * data.dds[2]
  

**Convert Function Source**

No source available.

gradDensityMagnitude
++++++++++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{4}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _gradDensityMagnitude(field, data):
      return np.sqrt(data["gradDensityX"]**2 +
                     data["gradDensityY"]**2 +
                     data["gradDensityZ"]**2)
  

**Convert Function Source**

No source available.

gradDensityX
++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{4}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _gradDensityX(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["density"].shape, dtype='float64')
      ds = div_fac * data['dx'].flat[0]
      new_field[1:-1,1:-1,1:-1]  = data["density"][sl_right,1:-1,1:-1]/ds
      new_field[1:-1,1:-1,1:-1] -= data["density"][sl_left ,1:-1,1:-1]/ds
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertgradDensity(data):
      return 1.0/data.convert("cm")
  

gradDensityY
++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{4}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _gradDensityY(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["density"].shape, dtype='float64')
      ds = div_fac * data['dy'].flat[0]
      new_field[1:-1,1:-1,1:-1]  = data["density"][1:-1,sl_right,1:-1]/ds
      new_field[1:-1,1:-1,1:-1] -= data["density"][1:-1,sl_left ,1:-1]/ds
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertgradDensity(data):
      return 1.0/data.convert("cm")
  

gradDensityZ
++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{4}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _gradDensityZ(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["density"].shape, dtype='float64')
      ds = div_fac * data['dz'].flat[0]
      new_field[1:-1,1:-1,1:-1]  = data["density"][1:-1,1:-1,sl_right]/ds
      new_field[1:-1,1:-1,1:-1] -= data["density"][1:-1,1:-1,sl_left ]/ds
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertgradDensity(data):
      return 1.0/data.convert("cm")
  

gradPressureMagnitude
+++++++++++++++++++++

   * Units: :math:`\rm{dyne}/\rm{cm}^{3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _gradPressureMagnitude(field, data):
      return np.sqrt(data["gradPressureX"]**2 +
                     data["gradPressureY"]**2 +
                     data["gradPressureZ"]**2)
  

**Convert Function Source**

No source available.

gradPressureX
+++++++++++++

   * Units: :math:`\rm{dyne}/\rm{cm}^{3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _gradPressureX(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["Pressure"].shape, dtype='float64')
      ds = div_fac * data['dx'].flat[0]
      new_field[1:-1,1:-1,1:-1]  = data["Pressure"][sl_right,1:-1,1:-1]/ds
      new_field[1:-1,1:-1,1:-1] -= data["Pressure"][sl_left ,1:-1,1:-1]/ds
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertgradPressure(data):
      return 1.0/data.convert("cm")
  

gradPressureY
+++++++++++++

   * Units: :math:`\rm{dyne}/\rm{cm}^{3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _gradPressureY(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["Pressure"].shape, dtype='float64')
      ds = div_fac * data['dy'].flat[0]
      new_field[1:-1,1:-1,1:-1]  = data["Pressure"][1:-1,sl_right,1:-1]/ds
      new_field[1:-1,1:-1,1:-1] -= data["Pressure"][1:-1,sl_left ,1:-1]/ds
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertgradPressure(data):
      return 1.0/data.convert("cm")
  

gradPressureZ
+++++++++++++

   * Units: :math:`\rm{dyne}/\rm{cm}^{3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _gradPressureZ(field, data):
      # We need to set up stencils
      if data.pf["HydroMethod"] == 2:
          sl_left = slice(None,-2,None)
          sl_right = slice(1,-1,None)
          div_fac = 1.0
      else:
          sl_left = slice(None,-2,None)
          sl_right = slice(2,None,None)
          div_fac = 2.0
      new_field = np.zeros(data["Pressure"].shape, dtype='float64')
      ds = div_fac * data['dz'].flat[0]
      new_field[1:-1,1:-1,1:-1]  = data["Pressure"][1:-1,1:-1,sl_right]/ds
      new_field[1:-1,1:-1,1:-1] -= data["Pressure"][1:-1,1:-1,sl_left ]/ds
      return new_field
  

**Convert Function Source**

.. code-block:: python

  def _convertgradPressure(data):
      return 1.0/data.convert("cm")
  

particle_density
++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _pdensity(field, data):
      blank = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return blank
      CICDeposit_3(data["particle_position_x"].astype(np.float64),
                   data["particle_position_y"].astype(np.float64),
                   data["particle_position_z"].astype(np.float64),
                   data["ParticleMass"],
                   data["particle_position_x"].size,
                   blank, np.array(data.LeftEdge).astype(np.float64),
                   np.array(data.ActiveDimensions).astype(np.int32),
                   np.float64(data['dx']))
      np.divide(blank, data["CellVolume"], blank)
      return blank
  

**Convert Function Source**

No source available.

particle_position_x
+++++++++++++++++++

   * Units: :math:`UNDEFINED`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def NullFunc(field, data):
      return
  

**Convert Function Source**

No source available.

particle_position_y
+++++++++++++++++++

   * Units: :math:`UNDEFINED`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def NullFunc(field, data):
      return
  

**Convert Function Source**

No source available.

particle_position_z
+++++++++++++++++++

   * Units: :math:`UNDEFINED`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def NullFunc(field, data):
      return
  

**Convert Function Source**

No source available.

sph_phi
+++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _sph_phi(field, data):
      center = data.get_field_parameter("center")
      normal = data.get_field_parameter("normal")
      
      coords = obtain_rvec(data)
  
      return get_sph_phi(coords, normal)
  

**Convert Function Source**

No source available.

sph_r
+++++

   * Units: :math:`\rm{cm}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _sph_r(field, data):
      center = data.get_field_parameter("center")
        
      coords = obtain_rvec(data)
  
      return get_sph_r(coords)
  

**Convert Function Source**

.. code-block:: python

  def _Convert_sph_r_CGS(data):
     return data.convert("cm")
  

sph_theta
+++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _sph_theta(field, data):
      center = data.get_field_parameter("center")
      normal = data.get_field_parameter("normal")
      
      coords = obtain_rvec(data)
  
      return get_sph_theta(coords, normal)
  

**Convert Function Source**

No source available.

tempContours
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Contours(field, data):
      return -np.ones_like(data["Ones"])
  

**Convert Function Source**

No source available.

x
+

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _coordX(field, data):
      dim = data.ActiveDimensions[0]
      return (np.ones(data.ActiveDimensions, dtype='float64')
                     * np.arange(data.ActiveDimensions[0])[:,None,None]
              +0.5) * data['dx'] + data.LeftEdge[0]
  

**Convert Function Source**

No source available.

y
+

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _coordY(field, data):
      dim = data.ActiveDimensions[1]
      return (np.ones(data.ActiveDimensions, dtype='float64')
                     * np.arange(data.ActiveDimensions[1])[None,:,None]
              +0.5) * data['dy'] + data.LeftEdge[1]
  

**Convert Function Source**

No source available.

z
+

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _coordZ(field, data):
      dim = data.ActiveDimensions[2]
      return (np.ones(data.ActiveDimensions, dtype='float64')
                     * np.arange(data.ActiveDimensions[2])[None,None,:]
              +0.5) * data['dz'] + data.LeftEdge[2]
  

**Convert Function Source**

No source available.

.. _enzo-field-names:

Enzo-Specific Field List
------------------------

Bmag
++++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Bmag(field, data):
      """ magnitude of bvec
      """
      return np.sqrt(data['Bx']**2 + data['By']**2 + data['Bz']**2)
  

**Convert Function Source**

No source available.

Comoving_DII_Density
++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_DI_Density
+++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_Electron_Density
+++++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_H2II_Density
+++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_H2I_Density
++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_HDI_Density
++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_HII_Density
++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_HI_Density
+++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_HM_Density
+++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_HeIII_Density
++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_HeII_Density
+++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_HeI_Density
++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_MetalSNIa_Density
++++++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_Metal_Density
++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

Comoving_PreShock_Density
+++++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesComovingDensity(field, data):
      sp = field.name.split("_")[0] + "_Density"
      ef = (1.0 + data.pf.current_redshift)**3.0
      return data[sp] / ef
  

**Convert Function Source**

No source available.

DII_Fraction
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

DII_Mass
++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

DII_MassMsun
++++++++++++

   * Units: :math:`M_{\odot}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

.. code-block:: python

  def _convertCellMassMsun(data):
      return 5.027854e-34 # g^-1
  

DII_NumberDensity
+++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesNumberDensity(field, data):
      species = field.name.split("_")[0]
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / _speciesMass[species]
  

**Convert Function Source**

.. code-block:: python

  def _ConvertNumberDensity(data):
      return 1.0/mh
  

DI_Fraction
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

DI_Mass
+++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

Dark_Matter_Mass
++++++++++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Dark_Matter_Mass(field, data):
      return data['Dark_Matter_Density'] * data["CellVolume"]
  

**Convert Function Source**

No source available.

Dark_Matter_MassMsun
++++++++++++++++++++

   * Units: :math:`M_{\odot}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Dark_Matter_Mass(field, data):
      return data['Dark_Matter_Density'] * data["CellVolume"]
  

**Convert Function Source**

.. code-block:: python

  def _convertCellMassMsun(data):
      return 5.027854e-34 # g^-1
  

Electron_Fraction
+++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

Electron_Mass
+++++++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

Gas_Energy
++++++++++

   * Units: :math:`\rm{ergs}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Gas_Energy(field, data):
      return data["GasEnergy"] / _convertEnergy(data)
  

**Convert Function Source**

.. code-block:: python

  def _convertEnergy(data):
      return data.convert("x-velocity")**2.0
  

H2II_Fraction
+++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

H2II_Mass
+++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

H2I_Fraction
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

H2I_Mass
++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

HDI_Fraction
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

HDI_Mass
++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

HII_Fraction
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

HII_Mass
++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

HI_Fraction
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

HI_Mass
+++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

HM_Fraction
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

HM_Mass
+++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

H_NumberDensity
+++++++++++++++

   * Units: :math:`\rm{cm}^{-3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _H_NumberDensity(field, data):
      field_data = np.zeros(data["density"].shape,
                            dtype=data["density"].dtype)
      if data.pf.parameters["MultiSpecies"] == 0:
          field_data += data["density"] * \
            data.pf.parameters["HydrogenFractionByMass"]
      if data.pf.parameters["MultiSpecies"] > 0:
          field_data += data["HI_Density"]
          field_data += data["HII_Density"]
      if data.pf.parameters["MultiSpecies"] > 1:
          field_data += data["HM_Density"]
          field_data += data["H2I_Density"]
          field_data += data["H2II_Density"]
      if data.pf.parameters["MultiSpecies"] > 2:
          field_data += data["HDI_Density"] / 2.0
      return field_data
  

**Convert Function Source**

.. code-block:: python

  def _ConvertNumberDensity(data):
      return 1.0/mh
  

HeIII_Fraction
++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

HeIII_Mass
++++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

HeII_Fraction
+++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

HeII_Mass
+++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

HeI_Fraction
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

HeI_Mass
++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

IsStarParticle
++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _IsStarParticle(field, data):
      is_star = (data['creation_time'] > 0).astype('float64')
      return is_star
  

**Convert Function Source**

No source available.

KineticEnergy
+++++++++++++

   * Units: :math:`\rm{ergs}/\rm{cm^3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _KineticEnergy(field, data):
      return 0.5*data["density"] * ( data["x-velocity"]**2.0
                                     + data["y-velocity"]**2.0
                                     + data["z-velocity"]**2.0 )
  

**Convert Function Source**

No source available.

MetalSNIa_Fraction
++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

MetalSNIa_Mass
++++++++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

Metal_Fraction
++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

Metal_Mass
++++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

Metallicity
+++++++++++

   * Units: :math:`Z_{\rm{\odot}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Metallicity(field, data):
      return data["Metal_Fraction"]
  

**Convert Function Source**

.. code-block:: python

  def _ConvertMetallicity(data):
      return 49.0196 # 1 / 0.0204
  

Metallicity3
++++++++++++

   * Units: :math:`Z_{\rm{\odot}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Metallicity3(field, data):
      return data["SN_Colour"]/data["density"]
  

**Convert Function Source**

.. code-block:: python

  def _ConvertMetallicity(data):
      return 49.0196 # 1 / 0.0204
  

NumberDensity
+++++++++++++

   * Units: :math:`\rm{cm}^{-3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _NumberDensity(field, data):
      # We can assume that we at least have Density
      # We should actually be guaranteeing the presence of a .shape attribute,
      # but I am not currently implementing that
      fieldData = np.zeros(data["density"].shape,
                           dtype = data["density"].dtype)
      if data.pf["MultiSpecies"] == 0:
          if data.has_field_parameter("mu"):
              mu = data.get_field_parameter("mu")
          else:
              mu = 0.6
          fieldData += data["density"] / mu
      if data.pf["MultiSpecies"] > 0:
          fieldData += data["HI_Density"] / 1.0
          fieldData += data["HII_Density"] / 1.0
          fieldData += data["HeI_Density"] / 4.0
          fieldData += data["HeII_Density"] / 4.0
          fieldData += data["HeIII_Density"] / 4.0
          fieldData += data["Electron_Density"] / 1.0
      if data.pf["MultiSpecies"] > 1:
          fieldData += data["HM_Density"] / 1.0
          fieldData += data["H2I_Density"] / 2.0
          fieldData += data["H2II_Density"] / 2.0
      if data.pf["MultiSpecies"] > 2:
          fieldData += data["DI_Density"] / 2.0
          fieldData += data["DII_Density"] / 2.0
          fieldData += data["HDI_Density"] / 3.0
      return fieldData
  

**Convert Function Source**

.. code-block:: python

  def _ConvertNumberDensity(data):
      return 1.0/mh
  

ParticleAge
+++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleAge(field, data):
      current_time = data.pf.current_time
      return (current_time - data["creation_time"])
  

**Convert Function Source**

.. code-block:: python

  def _convertParticleAge(data):
      return data.convert("years")
  

ParticleMass
++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleMass(field, data):
      particles = data["particle_mass"].astype('float64') * \
                  just_one(data["CellVolumeCode"].ravel())
      # Note that we mandate grid-type here, so this is okay
      return particles
  

**Convert Function Source**

.. code-block:: python

  def _convertParticleMass(data):
      return data.convert("density")*(data.convert("cm")**3.0)
  

ParticleMassMsun
++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleMass(field, data):
      particles = data["particle_mass"].astype('float64') * \
                  just_one(data["CellVolumeCode"].ravel())
      # Note that we mandate grid-type here, so this is okay
      return particles
  

**Convert Function Source**

.. code-block:: python

  def _convertParticleMassMsun(data):
      return data.convert("density")*((data.convert("cm")**3.0)/1.989e33)
  

PreShock_Fraction
+++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesFraction(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] / data["density"]
  

**Convert Function Source**

No source available.

PreShock_Mass
+++++++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _SpeciesMass(field, data):
      sp = field.name.split("_")[0] + "_Density"
      return data[sp] * data["CellVolume"]
  

**Convert Function Source**

No source available.

RadiationAcceleration
+++++++++++++++++++++

   * Units: :math:`\rm{cm} \rm{s}^{-2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _RadiationAccelerationMagnitude(field, data):
      return ( data["RadAccel1"]**2 + data["RadAccel2"]**2 +
               data["RadAccel3"]**2 )**(1.0/2.0)
  

**Convert Function Source**

No source available.

StarAgeYears
++++++++++++

   * Units: :math:`\rm{yr}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _StarAge(field, data):
      star_age = np.zeros(data['StarCreationTimeYears'].shape)
      with_stars = data['StarCreationTimeYears'] > 0
      star_age[with_stars] = data.pf.time_units['years'] * \
          data.pf.current_time - \
          data['StarCreationTimeYears'][with_stars]
      return star_age
  

**Convert Function Source**

No source available.

StarCreationTimeYears
+++++++++++++++++++++

   * Units: :math:`\rm{yr}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _StarCreationTime(field, data):
      return data['star_creation_time']
  

**Convert Function Source**

.. code-block:: python

  def _ConvertEnzoTimeYears(data):
      return data.pf.time_units['years']
  

StarDynamicalTimeYears
++++++++++++++++++++++

   * Units: :math:`\rm{yr}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _StarDynamicalTime(field, data):
      return data['star_dynamical_time']
  

**Convert Function Source**

.. code-block:: python

  def _ConvertEnzoTimeYears(data):
      return data.pf.time_units['years']
  

StarMetallicity
+++++++++++++++

   * Units: :math:`Z_{\rm{\odot}}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _StarMetallicity(field, data):
      return data['star_metallicity_fraction']
  

**Convert Function Source**

.. code-block:: python

  def _ConvertMetallicity(data):
      return 49.0196 # 1 / 0.0204
  

ThermalEnergy
+++++++++++++

   * Units: :math:`\rm{ergs}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _ThermalEnergy(field, data):
      if data.pf["HydroMethod"] == 2:
          return data["TotalEnergy"]
      
      if data.pf["DualEnergyFormalism"]:
          return data["GasEnergy"]
  
      if data.pf["HydroMethod"] in (4,6):
          return data["TotalEnergy"] - 0.5*(
              data["x-velocity"]**2.0
              + data["y-velocity"]**2.0
              + data["z-velocity"]**2.0 ) \
              - data["MagneticEnergy"]/data["density"]
  
      return data["TotalEnergy"] - 0.5*(
          data["x-velocity"]**2.0
          + data["y-velocity"]**2.0
          + data["z-velocity"]**2.0 )
  

**Convert Function Source**

No source available.

TotalEnergy
+++++++++++

   * Units: :math:`\rm{ergs}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _TotalEnergy(field, data):
      return data["Total_Energy"] / _convertEnergy(data)
  

**Convert Function Source**

.. code-block:: python

  def _convertEnergy(data):
      return data.convert("x-velocity")**2.0
  

Total_Energy
++++++++++++

   * Units: :math:`\rm{ergs}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Total_Energy(field, data):
      return data["TotalEnergy"] / _convertEnergy(data)
  

**Convert Function Source**

.. code-block:: python

  def _convertEnergy(data):
      return data.convert("x-velocity")**2.0
  

cic_particle_velocity_x
+++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cic_particle_field(field, data):
      """
      Create a grid field for particle quantities weighted by particle mass, 
      using cloud-in-cell deposit.
      """
      particle_field = field.name[4:]
      top = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return top
      particle_field_data = data[particle_field] * data['particle_mass']
      amr_utils.CICDeposit_3(data["particle_position_x"].astype(np.float64),
                             data["particle_position_y"].astype(np.float64),
                             data["particle_position_z"].astype(np.float64),
                             particle_field_data,
                             data["particle_position_x"].size,
                             top, np.array(data.LeftEdge).astype(np.float64),
                             np.array(data.ActiveDimensions).astype(np.int32), 
                             np.float64(data['dx']))
      del particle_field_data
  
      bottom = np.zeros(data.ActiveDimensions, dtype='float64')
      amr_utils.CICDeposit_3(data["particle_position_x"].astype(np.float64),
                             data["particle_position_y"].astype(np.float64),
                             data["particle_position_z"].astype(np.float64),
                             data["particle_mass"],
                             data["particle_position_x"].size,
                             bottom, np.array(data.LeftEdge).astype(np.float64),
                             np.array(data.ActiveDimensions).astype(np.int32), 
                             np.float64(data['dx']))
      top[bottom == 0] = 0.0
      bnz = bottom.nonzero()
      top[bnz] /= bottom[bnz]
      return top
  

**Convert Function Source**

No source available.

cic_particle_velocity_y
+++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cic_particle_field(field, data):
      """
      Create a grid field for particle quantities weighted by particle mass, 
      using cloud-in-cell deposit.
      """
      particle_field = field.name[4:]
      top = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return top
      particle_field_data = data[particle_field] * data['particle_mass']
      amr_utils.CICDeposit_3(data["particle_position_x"].astype(np.float64),
                             data["particle_position_y"].astype(np.float64),
                             data["particle_position_z"].astype(np.float64),
                             particle_field_data,
                             data["particle_position_x"].size,
                             top, np.array(data.LeftEdge).astype(np.float64),
                             np.array(data.ActiveDimensions).astype(np.int32), 
                             np.float64(data['dx']))
      del particle_field_data
  
      bottom = np.zeros(data.ActiveDimensions, dtype='float64')
      amr_utils.CICDeposit_3(data["particle_position_x"].astype(np.float64),
                             data["particle_position_y"].astype(np.float64),
                             data["particle_position_z"].astype(np.float64),
                             data["particle_mass"],
                             data["particle_position_x"].size,
                             bottom, np.array(data.LeftEdge).astype(np.float64),
                             np.array(data.ActiveDimensions).astype(np.int32), 
                             np.float64(data['dx']))
      top[bottom == 0] = 0.0
      bnz = bottom.nonzero()
      top[bnz] /= bottom[bnz]
      return top
  

**Convert Function Source**

No source available.

cic_particle_velocity_z
+++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _cic_particle_field(field, data):
      """
      Create a grid field for particle quantities weighted by particle mass, 
      using cloud-in-cell deposit.
      """
      particle_field = field.name[4:]
      top = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return top
      particle_field_data = data[particle_field] * data['particle_mass']
      amr_utils.CICDeposit_3(data["particle_position_x"].astype(np.float64),
                             data["particle_position_y"].astype(np.float64),
                             data["particle_position_z"].astype(np.float64),
                             particle_field_data,
                             data["particle_position_x"].size,
                             top, np.array(data.LeftEdge).astype(np.float64),
                             np.array(data.ActiveDimensions).astype(np.int32), 
                             np.float64(data['dx']))
      del particle_field_data
  
      bottom = np.zeros(data.ActiveDimensions, dtype='float64')
      amr_utils.CICDeposit_3(data["particle_position_x"].astype(np.float64),
                             data["particle_position_y"].astype(np.float64),
                             data["particle_position_z"].astype(np.float64),
                             data["particle_mass"],
                             data["particle_position_x"].size,
                             bottom, np.array(data.LeftEdge).astype(np.float64),
                             np.array(data.ActiveDimensions).astype(np.int32), 
                             np.float64(data['dx']))
      top[bottom == 0] = 0.0
      bnz = bottom.nonzero()
      top[bnz] /= bottom[bnz]
      return top
  

**Convert Function Source**

No source available.

dm_density
++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^3`
   * Projected Units: :math:`\rm{g}/\rm{cm}^2`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _dmpdensity(field, data):
      blank = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return blank
      if 'creation_time' in data.pf.field_info:
          filter = data['creation_time'] <= 0.0
          if not filter.any(): return blank
          num = filter.sum()
      else:
          filter = Ellipsis
          num = data["particle_position_x"].size
      amr_utils.CICDeposit_3(data["particle_position_x"][filter].astype(np.float64),
                             data["particle_position_y"][filter].astype(np.float64),
                             data["particle_position_z"][filter].astype(np.float64),
                             data["particle_mass"][filter],
                             num,
                             blank, np.array(data.LeftEdge).astype(np.float64),
                             np.array(data.ActiveDimensions).astype(np.int32), 
                             np.float64(data['dx']))
      return blank
  

**Convert Function Source**

.. code-block:: python

  def _convertDensity(data):
      return data.convert("density")
  

particle_mass
+++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          try:
              return io._read_data_set(data, p_field).astype(dtype)
          except io._read_exception:
              pass
          # This is bad.  But it's the best idea I have right now.
          return data._read_data(p_field.replace("_"," ")).astype(dtype)
  

**Convert Function Source**

No source available.

star_creation_time
++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _star_field(field, data):
      """
      Create a grid field for star quantities, weighted by star mass.
      """
      particle_field = field.name[5:]
      top = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return top
      filter = data['creation_time'] > 0.0
      if not filter.any(): return top
      particle_field_data = data[particle_field][filter] * data['particle_mass'][filter]
      amr_utils.CICDeposit_3(data["particle_position_x"][filter].astype(np.float64),
                            data["particle_position_y"][filter].astype(np.float64),
                            data["particle_position_z"][filter].astype(np.float64),
                            particle_field_data,
                            np.int64(np.where(filter)[0].size),
                            top, np.array(data.LeftEdge).astype(np.float64),
                            np.array(data.ActiveDimensions).astype(np.int32), 
                            np.float64(data['dx']))
      del particle_field_data
  
      bottom = np.zeros(data.ActiveDimensions, dtype='float64')
      amr_utils.CICDeposit_3(data["particle_position_x"][filter].astype(np.float64),
                            data["particle_position_y"][filter].astype(np.float64),
                            data["particle_position_z"][filter].astype(np.float64),
                            data["particle_mass"][filter],
                            np.int64(np.where(filter)[0].size),
                            bottom, np.array(data.LeftEdge).astype(np.float64),
                            np.array(data.ActiveDimensions).astype(np.int32), 
                            np.float64(data['dx']))
      top[bottom == 0] = 0.0
      bnz = bottom.nonzero()
      top[bnz] /= bottom[bnz]
      return top
  

**Convert Function Source**

No source available.

star_density
++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^3`
   * Projected Units: :math:`\rm{g}/\rm{cm}^2`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _spdensity(field, data):
      blank = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return blank
      filter = data['creation_time'] > 0.0
      if not filter.any(): return blank
      amr_utils.CICDeposit_3(data["particle_position_x"][filter].astype(np.float64),
                             data["particle_position_y"][filter].astype(np.float64),
                             data["particle_position_z"][filter].astype(np.float64),
                             data["particle_mass"][filter],
                             np.int64(np.where(filter)[0].size),
                             blank, np.array(data.LeftEdge).astype(np.float64),
                             np.array(data.ActiveDimensions).astype(np.int32), 
                             np.float64(data['dx']))
      return blank
  

**Convert Function Source**

.. code-block:: python

  def _convertDensity(data):
      return data.convert("density")
  

star_dynamical_time
+++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _star_field(field, data):
      """
      Create a grid field for star quantities, weighted by star mass.
      """
      particle_field = field.name[5:]
      top = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return top
      filter = data['creation_time'] > 0.0
      if not filter.any(): return top
      particle_field_data = data[particle_field][filter] * data['particle_mass'][filter]
      amr_utils.CICDeposit_3(data["particle_position_x"][filter].astype(np.float64),
                            data["particle_position_y"][filter].astype(np.float64),
                            data["particle_position_z"][filter].astype(np.float64),
                            particle_field_data,
                            np.int64(np.where(filter)[0].size),
                            top, np.array(data.LeftEdge).astype(np.float64),
                            np.array(data.ActiveDimensions).astype(np.int32), 
                            np.float64(data['dx']))
      del particle_field_data
  
      bottom = np.zeros(data.ActiveDimensions, dtype='float64')
      amr_utils.CICDeposit_3(data["particle_position_x"][filter].astype(np.float64),
                            data["particle_position_y"][filter].astype(np.float64),
                            data["particle_position_z"][filter].astype(np.float64),
                            data["particle_mass"][filter],
                            np.int64(np.where(filter)[0].size),
                            bottom, np.array(data.LeftEdge).astype(np.float64),
                            np.array(data.ActiveDimensions).astype(np.int32), 
                            np.float64(data['dx']))
      top[bottom == 0] = 0.0
      bnz = bottom.nonzero()
      top[bnz] /= bottom[bnz]
      return top
  

**Convert Function Source**

No source available.

star_metallicity_fraction
+++++++++++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _star_field(field, data):
      """
      Create a grid field for star quantities, weighted by star mass.
      """
      particle_field = field.name[5:]
      top = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return top
      filter = data['creation_time'] > 0.0
      if not filter.any(): return top
      particle_field_data = data[particle_field][filter] * data['particle_mass'][filter]
      amr_utils.CICDeposit_3(data["particle_position_x"][filter].astype(np.float64),
                            data["particle_position_y"][filter].astype(np.float64),
                            data["particle_position_z"][filter].astype(np.float64),
                            particle_field_data,
                            np.int64(np.where(filter)[0].size),
                            top, np.array(data.LeftEdge).astype(np.float64),
                            np.array(data.ActiveDimensions).astype(np.int32), 
                            np.float64(data['dx']))
      del particle_field_data
  
      bottom = np.zeros(data.ActiveDimensions, dtype='float64')
      amr_utils.CICDeposit_3(data["particle_position_x"][filter].astype(np.float64),
                            data["particle_position_y"][filter].astype(np.float64),
                            data["particle_position_z"][filter].astype(np.float64),
                            data["particle_mass"][filter],
                            np.int64(np.where(filter)[0].size),
                            bottom, np.array(data.LeftEdge).astype(np.float64),
                            np.array(data.ActiveDimensions).astype(np.int32), 
                            np.float64(data['dx']))
      top[bottom == 0] = 0.0
      bnz = bottom.nonzero()
      top[bnz] /= bottom[bnz]
      return top
  

**Convert Function Source**

No source available.

tracer_number_density
+++++++++++++++++++++

   * Units: :math:`\rm{particles}/\rm{kpc}^3`
   * Projected Units: :math:`\rm{particles}/\rm{kpc}^2`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _tpdensity(field, data): 
      blank = np.zeros(data.ActiveDimensions, dtype='float64')
      if data["particle_position_x"].size == 0: return blank
      filter = data['particle_type'] == 3 # tracer particles
      if not filter.any(): return blank
      amr_utils.CICDeposit_3(data["particle_position_x"][filter].astype(np.float64),
                             data["particle_position_y"][filter].astype(np.float64),
                             data["particle_position_z"][filter].astype(np.float64),
                             np.ones(filter.sum(), dtype="float64"),
                             np.int64(np.where(filter)[0].size),
                             blank, np.array(data.LeftEdge).astype(np.float64),
                             np.array(data.ActiveDimensions).astype(np.int32), 
                             np.float64(data['dx']))
      blank /= data['CellVolume']
      return blank
  

**Convert Function Source**

.. code-block:: python

  def _convertCmToKpc(data):
      return 1/(kpc_per_cm)**3
  

Orion-Specific Field List
-------------------------

Density
+++++++

   * Units: :math:`\rm{g}/\rm{cm}^3`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

ParticleMass
++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleMass(field, data):
      particles = data["particle_mass"].astype('float64')
      return particles
  

**Convert Function Source**

No source available.

ParticleMassMsun
++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleMassMsun(field, data):
      particles = data["particle_mass"].astype('float64')
      return particles/1.989e33
  

**Convert Function Source**

No source available.

Pressure
++++++++

   * Units: :math:`\rm{dyne}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Pressure(field,data):
      """M{(Gamma-1.0)*e, where e is thermal energy density
         NB: this will need to be modified for radiation
      """
      return (data.pf["Gamma"] - 1.0)*data["ThermalEnergy"]
  

**Convert Function Source**

No source available.

Temperature
+++++++++++

   * Units: :math:`\rm{Kelvin}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Temperature(field,data):
      return (data.pf["Gamma"]-1.0)*data.pf["mu"]*mh*data["ThermalEnergy"]/(kboltz*data["density"])
  

**Convert Function Source**

No source available.

ThermalEnergy
+++++++++++++

   * Units: :math:`\rm{ergs}/\rm{cm^3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _ThermalEnergy(field, data):
      """generate thermal (gas energy). Dual Energy Formalism was
          implemented by Stella, but this isn't how it's called, so I'll
          leave that commented out for now.
      """
      #if data.pf["DualEnergyFormalism"]:
      #    return data["GasEnergy"]
      #else:
      return data["TotalEnergy"] - 0.5 * data["density"] * (
          data["x-velocity"]**2.0
          + data["y-velocity"]**2.0
          + data["z-velocity"]**2.0 )
  

**Convert Function Source**

No source available.

TotalEnergy
+++++++++++

   * Units: :math:`\rm{erg}/\rm{cm}^3`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

particle_angmomen_x
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_angmomen_y
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_angmomen_z
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_burnstate
++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_id
+++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_luminosity
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_mass
+++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_mdeut
++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_mdot
+++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_mlast
++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_momentum_x
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_momentum_y
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_momentum_z
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_n
++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_position_x
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_position_y
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_position_z
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_r
++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

x-momentum
++++++++++

   * Units: :math:`\rm{g}/\rm{cm^2\ s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

x-velocity
++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

y-momentum
++++++++++

   * Units: :math:`\rm{gm}/\rm{cm^2\ s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

y-velocity
++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

z-momentum
++++++++++

   * Units: :math:`\rm{g}/\rm{cm^2\ s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

z-velocity
++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

FLASH-Specific Field List
-------------------------

Bx
++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Bx(fields, data):
      factor = GetMagRescalingFactor(data.pf)
      return data['magx']*factor
  

**Convert Function Source**

No source available.

By
++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _By(fields, data):
      factor = GetMagRescalingFactor(data.pf)
      return data['magy']*factor
  

**Convert Function Source**

No source available.

Bz
++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Bz(fields, data):
      factor = GetMagRescalingFactor(data.pf)
      return data['magz']*factor
  

**Convert Function Source**

No source available.

DII_Density
+++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

DII_Fraction
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

DI_Density
++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

DI_Fraction
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

Density
+++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

DivB
++++

   * Units: :math:`\rm{Gauss}\/\rm{cm}^{-1}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _DivB(fields, data):
      factor = GetMagRescalingFactor(data.pf)
      return data['divb']*factor
  

**Convert Function Source**

No source available.

Electron_Density
++++++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

Electron_Fraction
+++++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

Flame_Density
+++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

Flame_Fraction
++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

GasEnergy
+++++++++

   * Units: :math:`\rm{ergs}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _GasEnergy(fields, data) :
      return data["ThermalEnergy"]
  

**Convert Function Source**

No source available.

Grav_Potential
++++++++++++++

   * Units: :math:`\rm{ergs}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

H2II_Density
++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

H2II_Fraction
+++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

H2I_Density
+++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

H2I_Fraction
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

HD_Density
++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

HD_Fraction
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

HII_Density
+++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

HII_Fraction
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

HI_Density
++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

HI_Fraction
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

HM_Density
++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

HM_Fraction
+++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

HeIII_Density
+++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

HeIII_Fraction
++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

HeII_Density
++++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

HeII_Fraction
+++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

HeI_Density
+++++++++++

   * Units: :math:`\rm{g}/\rm{cm}^{3}`
   * Projected Units: :math:`\rm{g}/\rm{cm}^{2}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _dens(field, data):
          return data[fname] * data['density']
  

**Convert Function Source**

No source available.

HeI_Fraction
++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

NumberDensity
+++++++++++++

   * Units: :math:`\rm{cm}^{-3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _NumberDensity(fields,data) :
      try:
          return data["nele"]+data["nion"]
      except:
          pass
      return data['pres']/(data['temp']*kboltz)
  

**Convert Function Source**

No source available.

ParticleMass
++++++++++++

   * Units: :math:`\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

ParticleMassMsun
++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleMassMsun(field, data):
      return data["ParticleMass"]
  

**Convert Function Source**

.. code-block:: python

  def _convertParticleMassMsun(data):
      return 1.0/1.989e33
  

Pressure
++++++++

   * Units: :math:`\rm{erg}/\rm{cm}^{3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

Temperature
+++++++++++

   * Units: :math:`\rm{K}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

ThermalEnergy
+++++++++++++

   * Units: :math:`\rm{ergs}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _ThermalEnergy(fields, data) :
      try:
          return data["eint"]
      except:
          pass
      try:
          return data["Pressure"] / (data.pf["Gamma"] - 1.0) / data["density"]
      except:
          pass
      if data.has_field_parameter("mu") :
          mu = data.get_field_parameter("mu")
      else:
          mu = 0.6
      return kboltz*data["density"]*data["Temperature"]/(mu*mh) / (data.pf["Gamma"] - 1.0)
  

**Convert Function Source**

No source available.

TotalEnergy
+++++++++++

   * Units: :math:`\rm{ergs}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _TotalEnergy(fields, data) :
      try:
          etot = data["ener"]
      except:
          etot = data["ThermalEnergy"] + 0.5 * (
              data["x-velocity"]**2.0 +
              data["y-velocity"]**2.0 +
              data["z-velocity"]**2.0)
      try:
          etot += data['magp']/data["density"]
      except:
          pass
      return etot
  

**Convert Function Source**

No source available.

abar
++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _abar(field, data):
      try:
          return 1.0 / data['sumy']
      except:
          pass
      return data['dens']*Na*kboltz*data['temp']/data['pres']
  

**Convert Function Source**

No source available.

edens
+++++

   * Units: :math:`\rm{cm}^{-3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _nele(field, data):
      return data['dens'] * data['ye'] * Na
  

**Convert Function Source**

No source available.

nele
++++

   * Units: :math:`\rm{cm}^{-3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _nele(field, data):
      return data['dens'] * data['ye'] * Na
  

**Convert Function Source**

No source available.

nion
++++

   * Units: :math:`\rm{cm}^{-3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _nion(field, data):
      return data['dens'] * data['sumy'] * Na
  

**Convert Function Source**

No source available.

particle_index
++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

particle_position_x
+++++++++++++++++++

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

particle_position_y
+++++++++++++++++++

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

particle_position_z
+++++++++++++++++++

   * Units: :math:`\rm{cm}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

particle_velocity_x
+++++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

particle_velocity_y
+++++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

particle_velocity_z
+++++++++++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: True

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

x-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

y-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

z-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

Athena-Specific Field List
--------------------------

Bx
++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Bx(field, data):
      return data['cell_centered_B_x']
  

**Convert Function Source**

.. code-block:: python

  def _convertBfield(data):
          return np.sqrt(4*np.pi*data.convert("density")*data.convert("x-velocity")**2)
  

By
++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _By(field, data):
      return data['cell_centered_B_y']
  

**Convert Function Source**

.. code-block:: python

  def _convertBfield(data):
          return np.sqrt(4*np.pi*data.convert("density")*data.convert("x-velocity")**2)
  

Bz
++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Bz(field, data):
      return data['cell_centered_B_z']
  

**Convert Function Source**

.. code-block:: python

  def _convertBfield(data):
          return np.sqrt(4*np.pi*data.convert("density")*data.convert("x-velocity")**2)
  

Density
+++++++

   * Units: :math:`\rm{g}/\rm{cm}^3`
   * Projected Units: :math:`\rm{g}/\rm{cm}^2`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _density(field, data) :
      return data["density"]
  

**Convert Function Source**

.. code-block:: python

  def _convertDensity(data) :
      return data.convert("density")
  

Gas_Energy
++++++++++

   * Units: :math:`\rm{erg}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _gasenergy(field, data) :
      if "pressure" in data.pf.field_info:
          return data["pressure"]/(data.pf["Gamma"]-1.0)/data["density"]
      else:
          eint = data["total_energy"] - 0.5*(data["momentum_x"]**2 +
                                             data["momentum_y"]**2 +
                                             data["momentum_z"]**2)/data["density"]
          if "cell_centered_B_x" in data.pf.field_info:
              eint -= 0.5*(data["cell_centered_B_x"]**2 +
                           data["cell_centered_B_y"]**2 +
                           data["cell_centered_B_z"]**2)
          return eint/data["density"]
  

**Convert Function Source**

.. code-block:: python

  def _convertEnergy(data) :
      return data.convert("x-velocity")**2
  

Pressure
++++++++

   * Units: :math:`\rm{erg}/\rm{cm}^3`
   * Projected Units: :math:`\rm{erg}/\rm{cm}^2`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _pressure(field, data) :
      if "pressure" in data.pf.field_info:
          return data["pressure"]
      else:
          eint = data["total_energy"] - 0.5*(data["momentum_x"]**2 +
                                             data["momentum_y"]**2 +
                                             data["momentum_z"]**2)/data["density"]
          if "cell_centered_B_x" in data.pf.field_info:
              eint -= 0.5*(data["cell_centered_B_x"]**2 +
                           data["cell_centered_B_y"]**2 +
                           data["cell_centered_B_z"]**2)
          return eint*(data.pf["Gamma"]-1.0)
  

**Convert Function Source**

.. code-block:: python

  def _convertPressure(data) :
      return data.convert("density")*data.convert("x-velocity")**2
  

Temperature
+++++++++++

   * Units: :math:`\rm{K}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _temperature(field, data):
      if data.has_field_parameter("mu"):
          mu = data.get_field_parameter("mu")
      else:
          mu = 0.6
      return mu*mh*data["Pressure"]/data["density"]/kboltz
  

**Convert Function Source**

No source available.

x-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _xvelocity(field, data):
      if "velocity_x" in data.pf.field_info:
          return data["velocity_x"]
      else:
          return data["momentum_x"]/data["density"]           
  

**Convert Function Source**

.. code-block:: python

  def _convertVelocity(data):
      return data.convert("x-velocity")
  

y-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _yvelocity(field, data):
      if "velocity_y" in data.pf.field_info:
          return data["velocity_y"]
      else:
          return data["momentum_y"]/data["density"]
  

**Convert Function Source**

.. code-block:: python

  def _convertVelocity(data):
      return data.convert("x-velocity")
  

z-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _zvelocity(field, data):
      if "velocity_z" in data.pf.field_info:
          return data["velocity_z"]
      else:
          return data["momentum_z"]/data["density"]
  

**Convert Function Source**

.. code-block:: python

  def _convertVelocity(data):
      return data.convert("x-velocity")
  

Nyx-Specific Field List
-----------------------

Density
+++++++

   * Units: :math:`\rm{g} / \rm{cm}^3`
   * Projected Units: :math:`\rm{g} / \rm{cm}^2`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

ParticleMassMsun
++++++++++++++++

   * Units: :math:`\rm{M_{\odot}}`
   * Particle Type: True

**Field Source**

.. code-block:: python

  def _particle_mass_m_sun(field, data):
      return data["particle_mass"]
  

**Convert Function Source**

.. code-block:: python

  def _convertParticleMassMsun(data):
      return (1/1.989e33)
  

Pressure
++++++++

   * Units: :math:`\rm{M_{\odot}} (\rm{km} / \rm{s})^2 / \rm{Mpc}^3`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _pressure(field, data):
      """
      Computed using
  
      $$ pressure = (\gamma - 1.0) * e$$
  
      where e is thermal energy density. Note that this will need to be modified
      when radiation is accounted for.
  
      """
      return (data.pf["Gamma"] - 1.0) * data["ThermalEnergy"]
  

**Convert Function Source**

No source available.

Temperature
+++++++++++

   * Units: :math:`\rm{Kelvin}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _temperature(field, data):
      return ((data.pf["Gamma"] - 1.0) * data.pf["mu"] * mh *
              data["ThermalEnergy"] / (kboltz * data["density"]))
  

**Convert Function Source**

No source available.

ThermalEnergy
+++++++++++++

   * Units: :math:`\rm{M_{\odot}} (\rm{km} / \rm{s})^2`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _thermal_energy(field, data):
      """
      Generate thermal (gas energy). Dual Energy Formalism was implemented by
      Stella, but this isn't how it's called, so I'll leave that commented out for
      now.
  
      """
      #if data.pf["DualEnergyFormalism"]:
      #    return data["Gas_Energy"]
      #else:
      return data["Total_Energy"] - 0.5 * data["density"] * (
                                            data["x-velocity"]**2.0
                                          + data["y-velocity"]**2.0
                                          + data["z-velocity"]**2.0 )
  

**Convert Function Source**

No source available.

x-velocity
++++++++++

   * Units: :math:`\rm{km} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _x_velocity(field, data):
      """ Generate x-velocity from x-momentum and density. """
      return data["x-momentum"] / data["density"]
  

**Convert Function Source**

No source available.

y-velocity
++++++++++

   * Units: :math:`\rm{km} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _y_velocity(field, data):
      """ Generate y-velocity from y-momentum and density. """
      return data["y-momentum"] / data["density"]
  

**Convert Function Source**

No source available.

z-velocity
++++++++++

   * Units: :math:`\rm{km} / \rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _z_velocity(field, data):
      """ Generate z-velocity from z-momentum and density. """
      return data["z-momentum"] / data["density"]
  

**Convert Function Source**

No source available.

Chombo-Specific Field List
--------------------------

Bx
++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Bx(field,data):
      return data["X-magnfield"]
  

**Convert Function Source**

No source available.

By
++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _By(field,data):
      return data["Y-magnfield"]
  

**Convert Function Source**

No source available.

Bz
++

   * Units: :math:`\rm{Gauss}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Bz(field,data):
      return data["Z-magnfield"]
  

**Convert Function Source**

No source available.

Density
+++++++

   * Units: :math:`\rm{g}/\rm{cm^3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Density(field,data):
      """A duplicate of the density field. This is needed because when you try 
      to instantiate a PlotCollection without passing in a center, the code
      will try to generate one for you using the "density" field, which gives an error 
      if it isn't defined.
  
      """
      return data["density"]
  

**Convert Function Source**

No source available.

MagneticEnergy
++++++++++++++

   * Particle Type: False

**Field Source**

.. code-block:: python

  def _MagneticEnergy(field,data):
      return (data["X-magnfield"]**2 +
              data["Y-magnfield"]**2 +
              data["Z-magnfield"]**2)/2.
  

**Convert Function Source**

No source available.

ParticleMass
++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleMass(field, data):
      particles = data["particle_mass"].astype('float64')
      return particles
  

**Convert Function Source**

No source available.

ParticleMassMsun
++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

  def _ParticleMassMsun(field, data):
      particles = data["particle_mass"].astype('float64')
      return particles/1.989e33
  

**Convert Function Source**

No source available.

particle_angmomen_x
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_angmomen_y
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_angmomen_z
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_burnstate
++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_id
+++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_luminosity
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_mass
+++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_mdeut
++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_mdot
+++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_mlast
++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_momentum_x
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_momentum_y
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_momentum_z
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_n
++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_position_x
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_position_y
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_position_z
+++++++++++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

particle_r
++++++++++

   * Particle Type: True

**Field Source**

.. code-block:: python

      def _Particles(field, data):
          io = data.index.io
          if not data.NumberOfParticles > 0:
              return np.array([], dtype=dtype)
          else:
              return io._read_particles(data, p_field).astype(dtype)
  

**Convert Function Source**

No source available.

x-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _xVelocity(field, data):
      """ Generate x-velocity from x-momentum and density. """
      return data["X-momentum"]/data["density"]
  

**Convert Function Source**

No source available.

y-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _yVelocity(field,data):
      """ Generate y-velocity from y-momentum and density. """
      #try:
      #    return data["xvel"]
      #except KeyError:
      return data["Y-momentum"]/data["density"]
  

**Convert Function Source**

No source available.

z-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _zVelocity(field,data):
      """ Generate z-velocity from z-momentum and density. """
      return data["Z-momentum"]/data["density"]
  

**Convert Function Source**

No source available.

Pluto-Specific Field List
--------------------------

Density
+++++++

   * Units: :math:`\rm{g}/\rm{cm^3}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Density(field,data):
      """A duplicate of the density field. This is needed because when you try 
      to instantiate a PlotCollection without passing in a center, the code
      will try to generate one for you using the "density" field, which gives an error 
      if it isn't defined.
  
      """
      return data["rho"]
  

**Convert Function Source**

No source available.

X-momentum
++++++++++

   * Units: :math:`\rm{g}/\rm{cm^2 s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Xmomentum(field, data):
      """ Generate x-momentum. """
      return data["vx1"]*data["density"]
  

**Convert Function Source**

No source available.

Y-momentum
++++++++++

   * Units: :math:`\rm{g}/\rm{cm^2 s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Ymomentum(field, data):
      """ Generate y-momentum  """
      return data["vx2"]*data["density"]
  

**Convert Function Source**

No source available.

Z-Momentum
++++++++++

   * Units: :math:`\rm{g}/\rm{cm^2 s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

  def _Zmomentum(field,data):
      """ Generate z-momentum"""
      return data["vx3"]*data["density"]
  

**Convert Function Source**

No source available.

Grid-Data-Format-Specific Field List
------------------------------------

Density
+++++++

   * Units: :math:`\rm{g}/\rm{cm}^3`
   * Projected Units: :math:`\rm{g}/\rm{cm}^2`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

Pressure
++++++++

   * Units: :math:`\rm{erg}/\rm{g}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

x-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

y-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

z-velocity
++++++++++

   * Units: :math:`\rm{cm}/\rm{s}`
   * Particle Type: False

**Field Source**

.. code-block:: python

      def _TranslationFunc(field, data):
          return data[field_name]
  

**Convert Function Source**

No source available.

Generic-Format (Stream) Field List
----------------------------------

