"""
The strings we need for weaving.

@author: U{Matthew Turk<http://www.stanford.edu/~mturk/>}
@organization: U{KIPAC<http://www-group.slac.stanford.edu/KIPAC/>}
@contact: U{mturk@slac.stanford.edu<mailto:mturk@slac.stanford.edu>}
@license:
  Copyright (C) 2007 Matthew Turk.  All Rights Reserved.

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

ProfileBinningWeighted = \
"""
int bin;
int element;
int mybin;
npy_float64 myweight;          // should be higher precision
npy_float64 weights[num_bins]; // should be higher precision

for (bin = 0; bin < num_bins; bin++) {
  weights[bin] = 0;
  profilevalues(bin) = 0;
}

for (element = 0; element < num_elements ; element++) {
  mybin = binindices(element);
  profilevalues(mybin) += weightvalues(mybin)*fieldvalues(element);
  weights[mybin] += weightvalues(mybin);
}

for (bin = 0 ; bin < num_bins ; bin++ ) {
  profilevalues(bin) = profilevalues(bin) / weights[bin];
}
"""


ProfileBinningAccumulation = \
"""
int bin;
int element;
int mybin;
npy_float64 myweight;          // should be higher precision
npy_float64 weights[num_bins]; // should be higher precision

for (bin = 0; bin < num_bins; bin++) {
  weights[bin] = 0;
  profilevalues(bin) = 0;
}

for (element = 0; element < num_elements ; element++) {
  mybin = binindices(element);
  profilevalues(mybin) += fieldvalues(element);
}

for (bin = 1 ; bin < num_bins ; bin++ ) {
  profilevalues(bin) += profilevalues(bin-1);
}
"""

ProjectionRefineCoarseData = \
"""
/*
fpoints,  finedata_x,  finedata_y, finedata_vals, finedata_wgt,
cpoints,  coarsedata_x,  coarsedata_y, coarsedata_vals, coarsedata_wgt,
refinementFactor,  totalRefined

#define MIPCOMB(A,B) ((A) > (B) ? (A) : (B))
#define SUMCOMB(A,B) (A + B)

long fi, ci;

int flagged[cpoints], tr = 0;

long double rf = refinementFactor;

npy_int64 coarseCell_x, coarseCell_y;

for (ci=0; ci<cpoints; ci++) flagged[ci]=0;

for (fi = 0; fi < fpoints; fi++) {
    coarseCell_x = floorl((finedata_x(fi))/rf);
    coarseCell_y = floorl((finedata_y(fi))/rf);
    for (ci = 0; ci < cpoints; ci++) {
        if ((coarseCell_x == coarsedata_x(ci)) &&
            (coarseCell_y == coarsedata_y(ci))) {
                tr += 1;
                finedata_vals(fi) = %(COMBTYPE)sCOMB(coarsedata_vals[ci], finedata_vals[fi]);
                finedata_wgt(fi)  = %(COMBTYPE)sCOMB(coarsedata_wgt[ci],  finedata_wgt[fi]);
                flagged[ci] = 1;
                break;  // Each fine cell maps to one and only one coarse
                        // cell
        }
    }
}
*totalRefined = tr;

for (ci=0; ci<cpoints; ci++) {
    if (flagged[ci]==1) {
        coarsedata_x(ci)=-1;
        coarsedata_y(ci)=-1;
        coarsedata_vals(ci)=0.0;
        coarsedata_wgt(ci)=0.0;
    }
}
*/

"""

DataCubeRefineCoarseData = \
"""

// For every cell in fieldData, there are (dx_c / dx_f)**3.0 cells that it maps
// to in our cube, were dx_f is the fine resolution and dx_c is the coarse
// resolution.

int bx_f, by_f, bz_f, xc, yc, zc, xf, yf, zf;

for (xc = 0; xc < nxc; xc++)
  for (yc = 0; yc < nyc; yc++)
    for (zc = 0; zc < nzc; zc++) {
        // Now we have a cell.  What are the left edges?
        bx_f = (int) floorl((leftEdgeCoarse(0)+dx_c*xc - cubeLeftEdge(0))/dx_f);
        by_f = (int) floorl((leftEdgeCoarse(1)+dy_c*yc - cubeLeftEdge(1))/dy_f);
        bz_f = (int) floorl((leftEdgeCoarse(2)+dz_c*zc - cubeLeftEdge(2))/dz_f);
        if ((bx_f + rf > nxf) || (bx_f + rf < 0)
        ||  (by_f + rf > nyf) || (by_f + rf < 0)
        ||  (bz_f + rf > nzf) || (bz_f + rf < 0)) continue;
        for (xf = bx_f; xf < bx_f + rf ; xf++) {
          if (xf < 0 || xf > nxf) continue;
          for (yf = by_f; yf < by_f + rf ; yf++) {
            if (yf < 0 || yf > nyf) continue;
            for (zf = bz_f; zf < bz_f + rf ; zf++) {
              if (zf < 0 || zf > nzf) continue;
              cubeData(xf,yf,zf) = fieldData(xc,yc,zc);
            }
          }
        }
    }

"""
