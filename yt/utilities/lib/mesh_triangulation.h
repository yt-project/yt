#define MAX_NUM_TRI 12
#define HEX_NV 8
#define HEX_NT 12
#define TETRA_NV 4
#define TETRA_NT 4
#define WEDGE_NV 6
#define WEDGE_NT 8

// This array is used to triangulate the hexahedral mesh elements
// Each element has six faces with two triangles each.
// The vertex ordering convention is assumed to follow that used
// here: http://homepages.cae.wisc.edu/~tautges/papers/cnmev3.pdf
// Note that this is the case for Exodus II data.
int triangulate_hex[MAX_NUM_TRI][3] = {
  {0, 2, 1}, {0, 3, 2}, // Face is 3 2 1 0 
  {4, 5, 6}, {4, 6, 7}, // Face is 4 5 6 7
  {0, 1, 5}, {0, 5, 4}, // Face is 0 1 5 4
  {1, 2, 6}, {1, 6, 5}, // Face is 1 2 6 5
  {0, 7, 3}, {0, 4, 7}, // Face is 3 0 4 7
  {3, 6, 2}, {3, 7, 6}  // Face is 2 3 7 6
};

// Similarly, this is used to triangulate the tetrahedral cells
int triangulate_tetra[MAX_NUM_TRI][3] = {
  {0, 1, 3}, 
  {2, 3, 1},
  {0, 3, 2},
  {0, 2, 1},

  {-1, -1, -1},
  {-1, -1, -1},
  {-1, -1, -1},
  {-1, -1, -1},
  {-1, -1, -1},
  {-1, -1, -1},
  {-1, -1, -1},
  {-1, -1, -1}
};

// Triangulate wedges
int triangulate_wedge[MAX_NUM_TRI][3] = {
  {3, 0, 1},
  {4, 3, 1},
  {2, 5, 4},
  {2, 4, 1},
  {0, 3, 2},
  {2, 3, 5},
  {3, 4, 5},
  {0, 2, 1},

  {-1, -1, -1},
  {-1, -1, -1},
  {-1, -1, -1},
  {-1, -1, -1}
};


// This is used to select faces from a 20-sided hex element
int hex20_faces[6][8] = {
  {0, 1, 5, 4, 12, 8,  13, 16}, 
  {1, 2, 6, 5, 13, 9,  14, 17},
  {3, 2, 6, 7, 15, 10, 14, 18},
  {0, 3, 7, 4, 12, 11, 15, 19},
  {4, 5, 6, 7, 19, 16, 17, 18},
  {0, 1, 2, 3, 11, 8,  9,  10}
};
