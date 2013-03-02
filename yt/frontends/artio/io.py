"""
ARTIO-specific IO

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt-project.org/
License:
  Copyright (C) 2007-2011 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
from collections import defaultdict
import numpy as np

from yt.utilities.io_handler import \
    BaseIOHandler
from yt.utilities.logger import ytLogger as mylog
import yt.utilities.fortran_utils as fpu
import cStringIO

class IOHandlerARTIO(BaseIOHandler):
    _data_style = "artio"

    def _read_fluid_selection(self, chunks, selector, fields, size):
	print 'reading in fluid data'
        tr = dict((ftuple, np.empty(size, dtype='float64')) for ftuple in fields)
        cp = 0
        for onechunk in chunks:
            for subset in onechunk.objs:
                print 'reading from ',fields, subset.domain.grid_fn
                rv = subset.fill(fields) 
                for fieldtype, fieldname in fields:
                    mylog.debug("Filling %s with %s (%0.3e %0.3e) (%s:%s)",
                        fieldname, subset.masked_cell_count, rv[fieldname].min(), 
                        rv[fieldname].max(), cp, cp+subset.masked_cell_count)
                    tr[(fieldtype, fieldname)][cp:cp+subset.masked_cell_count] = rv.pop(fieldname)
                cp += subset.masked_cell_count
        return tr

    def _read_particle_selection(self, chunks, selector, fields):
        # http://yt-project.org/doc/analyzing/particles.html
        # ->creation_time >0 used to indicate star particles
	print "reading particle data"
        #
        # FIX need an input for particle type (in fields?)
#        accessed_species = ['N-BODY','STAR']
        accessed_species = ['STAR']
        #
	print 'io.py particle fields ',fields
        cp = 0
        tr = dict((ftuple, np.empty(0, dtype='float64')) for ftuple in fields)
        for onechunk in chunks:
            for subset in onechunk.objs:
                print 'reading values from', subset.domain.part_fn
                rv = subset.fill_particles( accessed_species, selector, fields)
                #get size and ensure that all fields have the same size (np)
                subset_size = 0 
                for fieldtype, fieldname in fields:
                    if subset_size != 0 and subset_size != len(rv[fieldname]) :
                        print 'size varies between fields! exiting'
                        sys.exit(1)
                    subset_size = len(rv[fieldname])
                cp_newsize = cp+subset_size
                for fieldtype, fieldname in fields:
                    tr[(fieldtype,fieldname)].resize(cp_newsize) 
                    tr[(fieldtype,fieldname)][cp:cp_newsize] = rv.pop(fieldname)
		cp = cp_newsize
        return tr

#stolen from frontends/art/
#All of these functions are to convert from hydro time var to 
#proper time
sqrt = np.sqrt
sign = np.sign

def find_root(f,a,b,tol=1e-6):
    c = (a+b)/2.0
    last = -np.inf
    assert(sign(f(a)) != sign(f(b)))  
    while np.abs(f(c)-last) > tol:
        last=f(c)
        if sign(last)==sign(f(b)):
            b=c
        else:
            a=c
        c = (a+b)/2.0
    return c

def quad(fintegrand,xmin,xmax,n=1e4):
    spacings = np.logspace(np.log10(xmin),np.log10(xmax),n)
    integrand_arr = fintegrand(spacings)
    val = np.trapz(integrand_arr,dx=np.diff(spacings))
    return val

def a2b(at,Om0=0.27,Oml0=0.73,h=0.700):
    def f_a2b(x):
        val = 0.5*sqrt(Om0) / x**3.0
        val /= sqrt(Om0/x**3.0 +Oml0 +(1.0 - Om0-Oml0)/x**2.0)
        return val
    #val, err = si.quad(f_a2b,1,at)
    val = quad(f_a2b,1,at)
    return val

def b2a(bt,**kwargs):
    #converts code time into expansion factor 
    #if Om0 ==1and OmL == 0 then b2a is (1 / (1-td))**2
    #if bt < -190.0 or bt > -.10:  raise 'bt outside of range'
    f_b2a = lambda at: a2b(at,**kwargs)-bt
    return find_root(f_b2a,1e-4,1.1)
    #return so.brenth(f_b2a,1e-4,1.1)
    #return brent.brent(f_b2a)

def a2t(at,Om0=0.27,Oml0=0.73,h=0.700):
    integrand = lambda x : 1./(x*sqrt(Oml0+Om0*x**-3.0))
    #current_time,err = si.quad(integrand,0.0,at,epsabs=1e-6,epsrel=1e-6)
    current_time = quad(integrand,1e-4,at)
    #spacings = np.logspace(-5,np.log10(at),1e5)
    #integrand_arr = integrand(spacings)
    #current_time = np.trapz(integrand_arr,dx=np.diff(spacings))
    current_time *= 9.779/h
    return current_time

def b2t(tb,n = 1e2,logger=None,**kwargs):
    tb = np.array(tb)
    if type(tb) == type(1.1): 
        return a2t(b2a(tb))
    if tb.shape == (): 
        return a2t(b2a(tb))
    if len(tb) < n: n= len(tb)
    age_min = a2t(b2a(tb.max(),**kwargs),**kwargs)
    age_max = a2t(b2a(tb.min(),**kwargs),**kwargs)
    tbs  = -1.*np.logspace(np.log10(-tb.min()),
                          np.log10(-tb.max()),n)
    ages = []
    for i,tbi in enumerate(tbs):
        ages += a2t(b2a(tbi)),
        if logger: logger(i)
    ages = np.array(ages)
    fb2t = np.interp(tb,tbs,ages)
    #fb2t = interp1d(tbs,ages)
    return fb2t

def spread_ages(ages,logger=None,spread=1.0e7*365*24*3600):
    #stars are formed in lumps; spread out the ages linearly
    da= np.diff(ages)
    assert np.all(da<=0)
    #ages should always be decreasing, and ordered so
    agesd = np.zeros(ages.shape)
    idx, = np.where(da<0)
    idx+=1 #mark the right edges
    #spread this age evenly out to the next age
    lidx=0
    lage=0
    for i in idx:
        n = i-lidx #n stars affected
        rage = ages[i]
        lage = max(rage-spread,0.0)
        agesd[lidx:i]=np.linspace(lage,rage,n)
        lidx=i
        #lage=rage
        if logger: logger(i)
    #we didn't get the last iter
    i=ages.shape[0]-1
    n = i-lidx #n stars affected
    rage = ages[i]
    lage = max(rage-spread,0.0)
    agesd[lidx:i]=np.linspace(lage,rage,n)
    return agesd
