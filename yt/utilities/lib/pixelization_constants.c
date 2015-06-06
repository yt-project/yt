/*******************************************************************************
# Copyright (c) 2015, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
*******************************************************************************/
//
// Some Cython versions don't like module-level constants, so we'll put them
// here.
//

#include "pixelization_constants.h"

/*
 Six faces, two vectors for each, two indices for each vector.  The function
 below unrolls how these are defined.  Some info can be found at:
 http://www.mscsoftware.com/training_videos/patran/Reverb_help/index.html#page/Finite%20Element%20Modeling/elem_lib_topics.16.8.html
 This is [6][2][2] in shape.
 Here are the faces and their four edges each:
 F1    1   2   3   4
 F2    5   6   7   8
 F3    1   10  5   9
 F4    2   11  6   10
 F5    3   12  7   11
 F6    4   9   8   12

 The edges are then defined by:
 E1    1 2
 E2    2 6
 E3    6 5
 E4    5 1
 E5    4 3
 E6    3 7
 E7    7 8
 E8    8 4
 E9    1 4
 E10   2 3
 E11   6 7
 E12   5 8
 Now we unroll these here ...
 */
const npy_uint8 hex_face_defs[MAX_NUM_FACES][2][2] = {
   /* Note that the first of each pair is the shared vertex */
   {{1, 0}, {1, 5}},
   {{2, 3}, {2, 6}},
   {{1, 0}, {1, 2}},
   {{5, 1}, {5, 6}},
   {{4, 5}, {4, 7}},
   {{0, 4}, {0, 3}},

   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}}
};

/* http://www.mscsoftware.com/training_videos/patran/Reverb_help/index.html#page/Finite%2520Element%2520Modeling/elem_lib_topics.16.6.html
 
  F1    1   2   3
  F2    1   5   4
  F3    2   6   5
  F4    3   4   6
 
  The edges are then defined by:
  E1    1   2
  E2    2   3
  E3    3   1
  E4    1   4
  E5    2   4
  E6    3   4
*/

const npy_uint8 tetra_face_defs[MAX_NUM_FACES][2][2] = {
   {{1, 0}, {1, 2}},
   {{1, 0}, {1, 3}},
   {{2, 1}, {2, 3}},
   {{3, 0}, {3, 2}},

   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}}
};

/* http://www.mscsoftware.com/training_videos/patran/Reverb_help/index.html#page/Finite%2520Element%2520Modeling/elem_lib_topics.16.7.html
  F1    1   2   3   *
  F2    4   5   6   *
  F3    1   8   4   7
  F4    2   9   5   8
  F5    3   7   6   9
 
  The edges are then defined by:
  E1    2   1
  E2    1   3
  E3    3   2
  E4    5   4
  E5    4   6
  E6    6   5
  E7    2   5
  E8    1   4
  E9    3   6
*/

const npy_uint8 wedge_face_defs[MAX_NUM_FACES][2][2] = {
   {{0, 1}, {0, 2}},
   {{3, 4}, {3, 5}},
   {{0, 1}, {0, 3}},
   {{2, 0}, {2, 5}},
   {{1, 2}, {1, 4}},

   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}},
   {{-1, -1}, {-1 -1}}
};
