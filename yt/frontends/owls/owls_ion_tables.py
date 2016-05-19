"""
OWLS ion tables

A module to handle the HM01 UV background spectra and ionization data from the
OWLS photoionization equilibrium lookup tables.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from yt.utilities.on_demand_imports import _h5py as h5py
import yt.extern.six as six
import numpy as np




def h5rd(fname, path, dtype=None):
    """ Read Data. Return a dataset located at <path> in file <fname> as
    a numpy array.
    e.g. rd( fname, '/PartType0/Coordinates' ). """

    data = None
    fid = h5py.h5f.open(six.b(fname), h5py.h5f.ACC_RDONLY)
    dg = h5py.h5d.open(fid, path.encode('ascii'))
    if dtype is None:
       dtype = dg.dtype
    data = np.zeros(dg.shape, dtype=dtype)
    dg.read(h5py.h5s.ALL, h5py.h5s.ALL, data)
    fid.close()
    return data


class IonTableSpectrum:

    """ A class to handle the HM01 spectra in the OWLS ionization tables. """

    def __init__(self, ion_file):

        where = '/header/spectrum/gammahi'
        self.GH1 = h5rd(ion_file, where) # GH1[1/s]

        where = '/header/spectrum/logenergy_ryd'
        self.logryd = h5rd(ion_file, where) # E[ryd]

        where = '/header/spectrum/logflux'
        self.logflux = h5rd(ion_file, where) # J[ergs/s/Hz/Sr/cm^2]

        where = '/header/spectrum/redshift'
        self.z = h5rd(ion_file, where) # z


    def return_table_GH1_at_z(self,z):

        # find redshift indices
        #-----------------------------------------------------------------
        i_zlo = np.argmin( np.abs( self.z - z ) )
        if self.z[i_zlo] < z:
            i_zhi = i_zlo + 1
        else:
            i_zhi = i_zlo
            i_zlo = i_zlo - 1

        z_frac = (z - self.z[i_zlo]) / (self.z[i_zhi] - self.z[i_zlo])

        # find GH1 from table
        #-----------------------------------------------------------------
        logGH1_all = np.log10( self.GH1 )
        dlog_GH1 = logGH1_all[i_zhi] - logGH1_all[i_zlo]

        logGH1_table = logGH1_all[i_zlo] + z_frac * dlog_GH1
        GH1_table = 10.0**logGH1_table

        return GH1_table


class IonTableOWLS:

    """ A class to handle OWLS ionization tables. """

    DELTA_nH = 0.25
    DELTA_T = 0.1

    def __init__(self, ion_file):

        self.ion_file = ion_file

        # ionbal is indexed like [nH, T, z]
        # nH and T are log quantities
        #---------------------------------------------------------------
        self.nH = h5rd( ion_file, '/logd' )         # log nH [cm^-3]
        self.T = h5rd( ion_file, '/logt' )          # log T [K]
        self.z = h5rd( ion_file, '/redshift' )      # z

        # read the ionization fractions
        # linear values stored in file so take log here
        # ionbal is the ionization balance (i.e. fraction)
        #---------------------------------------------------------------
        self.ionbal = h5rd( ion_file, '/ionbal' ).astype(np.float64)
        self.ionbal_orig = self.ionbal.copy()

        ipositive = self.ionbal > 0.0
        izero = np.logical_not(ipositive)
        self.ionbal[izero] = self.ionbal[ipositive].min()

        self.ionbal = np.log10( self.ionbal )


        # load in background spectrum
        #---------------------------------------------------------------
        self.spectrum = IonTableSpectrum( ion_file )

        # calculate the spacing along each dimension
        #---------------------------------------------------------------
        self.dnH = self.nH[1:] - self.nH[0:-1]
        self.dT = self.T[1:] - self.T[0:-1]
        self.dz = self.z[1:] - self.z[0:-1]

        self.order_str = '[log nH, log T, z]'


    # sets iz and fz
    #-----------------------------------------------------
    def set_iz( self, z ):

        if z <= self.z[0]:
            self.iz = 0
            self.fz = 0.0
        elif z >= self.z[-1]:
            self.iz = len(self.z) - 2
            self.fz = 1.0
        else:
            for iz in range( len(self.z)-1 ):
                if z < self.z[iz+1]:
                    self.iz = iz
                    self.fz = ( z - self.z[iz] ) / self.dz[iz]
                    break



    # interpolate the table at a fixed redshift for the input
    # values of nH and T ( input should be log ).  A simple
    # tri-linear interpolation is used.
    #-----------------------------------------------------
    def interp( self, nH, T ):

        nH = np.array( nH )
        T  = np.array( T )

        if nH.size != T.size:
            raise ValueError(' owls_ion_tables: array size mismatch !!! ')

        # field discovery will have nH.size == 1 and T.size == 1
        # in that case we simply return 1.0

        if nH.size == 1 and T.size == 1:
            ionfrac = 1.0
            return ionfrac


        # find inH and fnH
        #-----------------------------------------------------
        x_nH = ( nH - self.nH[0] ) / self.DELTA_nH
        x_nH_clip = np.clip( x_nH, 0.0, self.nH.size-1.001 )
        fnH,inH = np.modf( x_nH_clip )
        inH = inH.astype( np.int32 )


        # find iT and fT
        #-----------------------------------------------------
        x_T = ( T - self.T[0] ) / self.DELTA_T
        x_T_clip = np.clip( x_T, 0.0, self.T.size-1.001 )
        fT,iT = np.modf( x_T_clip )
        iT = iT.astype( np.int32 )


        # short names for previously calculated iz and fz
        #-----------------------------------------------------
        iz = self.iz
        fz = self.fz


        # calculate interpolated value
        # use tri-linear interpolation on the log values
        #-----------------------------------------------------

        ionfrac = self.ionbal[inH,   iT,   iz  ] * (1-fnH) * (1-fT) * (1-fz) + \
                  self.ionbal[inH+1, iT,   iz  ] * (fnH)   * (1-fT) * (1-fz) + \
                  self.ionbal[inH,   iT+1, iz  ] * (1-fnH) * (fT)   * (1-fz) + \
                  self.ionbal[inH,   iT,   iz+1] * (1-fnH) * (1-fT) * (fz)   + \
                  self.ionbal[inH+1, iT,   iz+1] * (fnH)   * (1-fT) * (fz)   + \
                  self.ionbal[inH,   iT+1, iz+1] * (1-fnH) * (fT)   * (fz)   + \
                  self.ionbal[inH+1, iT+1, iz]   * (fnH)   * (fT)   * (1-fz) + \
                  self.ionbal[inH+1, iT+1, iz+1] * (fnH)   * (fT)   * (fz)

        return 10**ionfrac
