"""
Halo filters to be used with the HaloProfiler.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2011 Britton Smith.  All Rights Reserved.

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

from copy import deepcopy
import numpy as na

from yt.funcs import *

def VirialFilter(profile, overdensity_field='ActualOverdensity',
                 virial_overdensity=200., must_be_virialized=True,
                 virial_filters=[['TotalMassMsun', '>=','1e14']],
                 virial_quantities=['TotalMassMsun', 'RadiusMpc'],
                 virial_index=None, use_log=False):
    """
    Filter halos by virial quantities.
    Return values are a True or False whether the halo passed the filter, 
    along with a dictionary of virial quantities for the fields specified in 
    the virial_quantities keyword.  Thresholds for virial quantities are 
    given with the virial_filters keyword in the following way: 
    [field, condition, value].
    :param: overdensity_field (str): the field used for interpolation with the 
    specified critical value given with 'virial_overdensity'.  
    Default: 'ActualOverdensity'.
    :param: virial_overdensity (flt): the value used for interpolation.  
    Default: 200.[['TotalMassMsun', '>=','1e14']]
    :param: must_be_virialized (bool): if no values in the profile are above the 
    value of virial_overdensity, the halo does not pass the filter.  
    Default: True.
    :param: virial_filters (list): conditional filters based on virial quantities 
    given in the following way: [field, condition, value].  
    Default: [['TotalMassMsun', '>=','1e14']].
    :param: virial_quantities (list): fields for which interpolated values should 
    be calculated and returned.  Default: ['TotalMassMsun', 'RadiusMpc'].
    :param: virial_index (list): if given as a list, the index of the radial profile 
    which is used for interpolation is placed here.  Default: None.
    :param: use_log (bool): if True, interpolation is done in log space.  
    Default: False.
    """

    fields = deepcopy(virial_quantities)
    if virial_filters is None: virial_filters = []
    for vfilter in virial_filters:
        if not vfilter[0] in fields:
            fields.append(vfilter[0])
    
    overDensity = []
    temp_profile = dict((field, []) for field in fields)

    for q in range(len(profile[overdensity_field])):
        good = True
        if (profile[overdensity_field][q] != profile[overdensity_field][q]):
            good = False
            continue
        for field in fields:
            if (profile[field][q] != profile[field][q]):
                good = False
                break
        if good:
            overDensity.append(profile[overdensity_field][q])
            for field in fields:
                temp_profile[field].append(profile[field][q])

    if use_log:
        for field in temp_profile.keys():
            temp_profile[field] = na.log10(temp_profile[field])

    virial = dict((field, 0.0) for field in fields)

    if (not (na.array(overDensity) >= virial_overdensity).any()) and \
            must_be_virialized:
        mylog.error("This halo is not virialized!")
        return [False, {}]

    if (len(overDensity) < 2):
        mylog.error("Skipping halo with no valid points in profile.")
        return [False, {}]

    if (overDensity[1] <= virial_overdensity):
        index = 0
    elif (overDensity[-1] >= virial_overdensity):
        index = -2
    else:
        for q in (na.arange(len(overDensity),0,-1)-1):
            if (overDensity[q] < virial_overdensity) and (overDensity[q-1] >= virial_overdensity):
                index = q - 1
                break

    if type(virial_index) is list:
        virial_index.append(index)

    for field in fields:
        if (overDensity[index+1] - overDensity[index]) == 0:
            mylog.error("Overdensity profile has slope of zero.")
            return [False, {}]
        else:
            slope = (temp_profile[field][index+1] - temp_profile[field][index]) / \
                (overDensity[index+1] - overDensity[index])
            value = slope * (virial_overdensity - overDensity[index]) + \
                temp_profile[field][index]
            virial[field] = value

    if use_log:
        for field in virial.keys():
            virial[field] = na.power(10, virial[field])

    for vfilter in virial_filters:
        if eval("%s %s %s" % (virial[vfilter[0]],vfilter[1],vfilter[2])):
            mylog.debug("(%s %s %s) returned True for %s." % \
                            (vfilter[0],vfilter[1],vfilter[2],virial[vfilter[0]]))
            continue
        else:
            mylog.debug("(%s %s %s) returned False for %s." % \
                            (vfilter[0],vfilter[1],vfilter[2],virial[vfilter[0]]))
            return [False, {}]

    return [True, dict((("%s_%s" % (q, virial_overdensity)), virial[q])
                       for q in virial_quantities)]

