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

iterate_over_contours = r"""
int i,j,k,n;
for(n=0; n<np; n++){
  i=xi(n);j=yi(n);k=zi(n);
  process_neighbors(fd, i,j,k, mi,mj,mk);
}
"""

recursively_find_contours = r"""
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

int process_neighbors(blitz::Array<double,3> fd, int i, int j, int k,
                       int mi, int mj, int mk) {
  int off_i, off_j, off_k, new_cid, spawn_check, original;
  using namespace std;
  original = fd(i,j,k);
  do {
    spawn_check = 0;
    for (off_i=MAX(i-1,0);off_i<=MIN(i+1,mi-1);off_i++)
      for (off_j=MAX(j-1,0);off_j<=MIN(j+1,mj-1);off_j++)
        for (off_k=MAX(k-1,0);off_k<=MIN(k+1,mk-1);off_k++) {
          if(off_i==i&&off_j==j&&off_k==k) continue;
          if(fd(off_i,off_j,off_k)!=-1) {
              if(fd(off_i,off_j,off_k) > fd(i,j,k)){
                fd(i,j,k) = fd(off_i,off_j,off_k);
                spawn_check += 1;
              }
              if(fd(off_i,off_j,off_k) < fd(i,j,k)){
                fd(off_i,off_j,off_k) = fd(i,j,k);
                new_cid = process_neighbors(fd, off_i, off_j, off_k, mi,mj,mk);
                if (new_cid != fd(i,j,k)) spawn_check += 1;
                fd(i,j,k) = new_cid;
              }
            }
          }
    } while (spawn_check > 0);
    return fd(i,j,k);
}
"""


check_cell_distance = \
r"""
using namespace std;
int i, j, k;
long double cell_dist, rad_m, rad_n;
k=0;
for(i=0;i<mp;i++){
  for(j=0;j<np;j++){
    /*
   cell_dist = sqrtl(pow(mx(i)-nx(j),2) +
                     pow(my(i)-ny(j),2) +
                     pow(mz(i)-nz(j),2));
   rad_m = sqrtl(pow(mdx(i),2) +
                 pow(mdy(i),2) +
                 pow(mdz(i),2));
   rad_n = sqrtl(pow(ndx(j),2) +
                 pow(ndy(j),2) +
                 pow(ndz(j),2));
    */
   //if(cell_dist > 1.01 * (rad_n/2.0+rad_m/2.0)) continue;
   if(fabsl(mx(i)-nx(j))>(mdx(i)+ndx(j))/2.0) continue;
   if(fabsl(my(i)-ny(j))>(mdy(i)+ndy(j))/2.0) continue;
   if(fabsl(mz(i)-nz(j))>(mdz(i)+ndz(j))/2.0) continue;
   k++;
   break;
   cout << cell_dist << "\t" << 1.01*(rad_n/2.0+rad_m/2.0) << "\t";
   cout << grid_ids_m(i) << "\t" << grid_ids_n(j) << endl;
  }
}
n_bad(0) += k;
"""
