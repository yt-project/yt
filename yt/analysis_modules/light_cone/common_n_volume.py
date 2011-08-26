"""
Function to calculate volume in common between two n-cubes, with optional 
periodic boundary conditions.

Author: Britton Smith <brittons@origins.colorado.edu>
Affiliation: CASA/University of CO, Boulder
Homepage: http://yt-project.org/
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

import numpy as na

def commonNVolume(nCube1,nCube2,periodic=None):
    "Return the n-volume in common between the two n-cubes."

    # Check for proper args.
    if ((len(na.shape(nCube1)) != 2) or
        (na.shape(nCube1)[1] != 2) or
        (na.shape(nCube1) != na.shape(nCube2))):
        print "Arguments must be 2 (n,2) numpy array."
        return 0

    if ((periodic is not None) and
        (na.shape(nCube1) != na.shape(periodic))):
        print "periodic argument must be (n,2) numpy array."
        return 0

    nCommon = 1.0
    for q in range(na.shape(nCube1)[0]):
        if (periodic is None):
            nCommon *= commonSegment(nCube1[q],nCube2[q])
        else:
            nCommon *= commonSegment(nCube1[q],nCube2[q],periodic=periodic[q])

    return nCommon

def commonSegment(seg1,seg2,periodic=None):
    "Return the length of the common segment."

    # Check for proper args.
    if ((len(seg1) != 2) or (len(seg2) != 2)):
        print "Arguments must be arrays of size 2."
        return 0

    # If not periodic, then this is very easy.
    if periodic is None:
        seg1.sort()
        len1 = seg1[1] - seg1[0]
        seg2.sort()
        len2 = seg2[1] - seg2[0]

        common = 0.0

        add = seg1[1] - seg2[0]
        if ((add > 0) and (add <= max(len1,len2))):
            common += add
        add = seg2[1] - seg1[0]
        if ((add > 0) and (add <= max(len1,len2))):
            common += add
        common = min(common,len1,len2)
        return common

    # If periodic, it's a little more complicated.
    else:
        if len(periodic) != 2:
            print "periodic array must be of size 2."
            return 0

        seg1.sort()
        flen1 = seg1[1] - seg1[0]
        len1 = flen1 - int(flen1)
        seg2.sort()
        flen2 = seg2[1] - seg2[0]
        len2 = flen2 - int(flen2)

        periodic.sort()
        scale = periodic[1] - periodic[0]

        if (abs(int(flen1)-int(flen2)) >= scale):
            return min(flen1,flen2)

        # Adjust for periodicity
        seg1[0] = na.mod(seg1[0],scale) + periodic[0]
        seg1[1] = seg1[0] + len1
        if (seg1[1] > periodic[1]): seg1[1] -= scale
        seg2[0] = na.mod(seg2[0],scale) + periodic[0]
        seg2[1] = seg2[0] + len2
        if (seg2[1] > periodic[1]): seg2[1] -= scale

        # create list of non-periodic segments
        pseg1 = []
        if (seg1[0] >= seg1[1]):
            pseg1.append([seg1[0],periodic[1]])
            pseg1.append([periodic[0],seg1[1]])
        else:
            pseg1.append(seg1)
        pseg2 = []
        if (seg2[0] >= seg2[1]):
            pseg2.append([seg2[0],periodic[1]])
            pseg2.append([periodic[0],seg2[1]])
        else:
            pseg2.append(seg2)

        # Add up common segments.
        common = min(int(flen1),int(flen2))

        for subseg1 in pseg1:
            for subseg2 in pseg2:
                common += commonSegment(subseg1,subseg2)

        return common
