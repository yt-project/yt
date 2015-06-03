// This array is used to triangulate the hexahedral mesh elements
// Each element has six faces with two triangles each.
// The vertex ordering convention is assumed to follow that used
// here: http://homepages.cae.wisc.edu/~tautges/papers/cnmev3.pdf
// Note that this is the case for Exodus II data.
int triangulate_hex[12][3] = {
  {0, 1, 2}, {0, 2, 3}, // Face is 0 1 2 3 
  {4, 5, 6}, {4, 6, 7}, // Face is 4 5 6 7
  {0, 1, 5}, {0, 5, 4}, // Face is 0 1 5 4
  {1, 2, 6}, {1, 6, 5}, // Face is 1 2 6 5
  {0, 3, 7}, {0, 7, 4}, // Face is 0 3 7 4
  {3, 2, 6}, {3, 6, 7}  // Face is 3 2 6 7
};

// Similarly, this is used to triangulate the tetrahedral cells
int triangulate_tetra[4][3] = {
  {0, 1, 2}, 
  {0, 1, 3},
  {0, 2, 3},
  {1, 2, 3}
};
