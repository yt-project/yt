"""
Halo filters to be used with the HaloProfiler.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008-2009 Britton Smith.  All Rights Reserved.

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

from yt.logger import lagosLogger as mylog
from copy import deepcopy
import numpy as na

def VirialFilter(profile,overdensity_field='ActualOverdensity',
                 virial_overdensity=200.,must_be_virialized=True,
                 virial_filters=[['TotalMassMsun','>=','1e14']],
                 virial_quantities=['TotalMassMsun','RadiusMpc'],
                 virial_index=None):
    """
    Filter halos by virial quantities.
    Return values are a True or False whether the halo passed the filter, 
    along with a dictionary of virial quantities for the fields specified in 
    the virial_quantities keyword.  Thresholds for virial quantities are 
    given with the virial_filters keyword in the following way: 
    [field, condition, value].
    """

    fields = deepcopy(virial_quantities)
    if virial_filters is None: virial_filters = []
    for vfilter in virial_filters:
        if not vfilter[0] in fields:
            fields.append(vfilter[0])
    
    overDensity = []
    temp_profile = {}
    for field in fields:
        temp_profile[field] = []

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

    virial = {}
    for field in fields:
        virial[field] = 0.0

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
        for q in (na.arange(len(overDensity)-2))+2:
            if (overDensity[q] < virial_overdensity):
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

    for vfilter in virial_filters:
        if eval("%s %s %s" % (virial[vfilter[0]],vfilter[1],vfilter[2])):
            mylog.debug("(%s %s %s) returned True for %s." % (vfilter[0],vfilter[1],vfilter[2],virial[vfilter[0]]))
            continue
        else:
            mylog.debug("(%s %s %s) returned False for %s." % (vfilter[0],vfilter[1],vfilter[2],virial[vfilter[0]]))
            return [False, {}]

    return [True, dict((q,virial[q]) for q in virial_quantities)]

