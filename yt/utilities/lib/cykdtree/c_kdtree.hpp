#include <vector>
#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include "c_utils.hpp"

#define LEAF_MAX 4294967295

template <typename T>
T deserialize_scalar(std::istream &is) {
  T scalar;
  is.read((char*)&scalar, sizeof(T));
  return scalar;
}

template <typename T>
void serialize_scalar(std::ostream &os, const T &scalar) {
  os.write((char*)&scalar, sizeof(scalar));
}

template <typename T>
T* deserialize_pointer_array(std::istream &is, uint64_t len) {
  T* arr = (T*)malloc(len*sizeof(T));
  is.read((char*)&arr[0], len*sizeof(T));
  return arr;
}

template <typename T>
void serialize_pointer_array(std::ostream &os, const T* array, uint64_t len) {
  os.write((char*)array, len*sizeof(T));
}

class Node {
public:
  bool is_empty;
  bool is_leaf;
  uint32_t leafid;
  uint32_t ndim;
  double *left_edge;
  double *right_edge;
  uint64_t left_idx;
  uint64_t children;
  bool *periodic_left;
  bool *periodic_right;
  std::vector<std::vector<uint32_t> > left_neighbors;
  std::vector<std::vector<uint32_t> > right_neighbors;
  std::vector<uint32_t> all_neighbors;
  std::vector<Node*> left_nodes;
  // innernode parameters
  uint32_t split_dim;
  double split;
  Node *less;
  Node *greater;
  // empty node constructor
  Node() {
    is_empty = true;
    is_leaf = false;
    leafid = LEAF_MAX;
    ndim = 0;
    left_edge = NULL;
    right_edge = NULL;
    periodic_left = NULL;
    periodic_right = NULL;
    less = NULL;
    greater = NULL;
  }
  // emtpy node with some info
  Node(uint32_t ndim0, double *le, double *re, bool *ple, bool *pre) {
    is_empty = true;
    is_leaf = false;
    leafid = 4294967295;
    ndim = ndim0;
    left_edge = (double*)malloc(ndim*sizeof(double));
    right_edge = (double*)malloc(ndim*sizeof(double));
    periodic_left = (bool*)malloc(ndim*sizeof(bool));
    periodic_right = (bool*)malloc(ndim*sizeof(bool));
    memcpy(left_edge, le, ndim*sizeof(double));
    memcpy(right_edge, re, ndim*sizeof(double));
    memcpy(periodic_left, ple, ndim*sizeof(bool));
    memcpy(periodic_right, pre, ndim*sizeof(bool));
    less = NULL;
    greater = NULL;
    for (uint32_t i=0; i<ndim; i++) {
      left_nodes.push_back(NULL);
    }
  }
  // innernode constructor
  Node(uint32_t ndim0, double *le, double *re, bool *ple, bool *pre,
       uint64_t Lidx, uint32_t sdim0, double split0, Node *lnode, Node *gnode,
       std::vector<Node*> left_nodes0) {
    is_empty = false;
    is_leaf = false;
    leafid = 4294967295;
    ndim = ndim0;
    left_idx = Lidx;

    split_dim = sdim0;
    split = split0;
    less = lnode;
    greater = gnode;
    children = lnode->children + gnode->children;

    left_edge = (double*)malloc(ndim*sizeof(double));
    right_edge = (double*)malloc(ndim*sizeof(double));
    periodic_left = (bool*)malloc(ndim*sizeof(bool));
    periodic_right = (bool*)malloc(ndim*sizeof(bool));
    memcpy(left_edge, le, ndim*sizeof(double));
    memcpy(right_edge, re, ndim*sizeof(double));
    memcpy(periodic_left, ple, ndim*sizeof(bool));
    memcpy(periodic_right, pre, ndim*sizeof(bool));
    for (uint32_t d = 0; d < ndim; d++)
      left_nodes.push_back(left_nodes0[d]);

    left_neighbors = std::vector<std::vector<uint32_t> >(ndim);
    right_neighbors = std::vector<std::vector<uint32_t> >(ndim);
  }
  // leafnode constructor
  Node(uint32_t ndim0, double *le, double *re, bool *ple, bool *pre,
       uint64_t Lidx, uint64_t n, int leafid0,
       std::vector<Node*> left_nodes0) {
    is_empty = false;
    is_leaf = true;
    leafid = leafid0;
    ndim = ndim0;
    split = 0.0;
    split_dim = 0;
    left_idx = Lidx;
    less = NULL;
    greater = NULL;

    children = n;

    left_edge = (double*)malloc(ndim*sizeof(double));
    right_edge = (double*)malloc(ndim*sizeof(double));
    periodic_left = (bool*)malloc(ndim*sizeof(bool));
    periodic_right = (bool*)malloc(ndim*sizeof(bool));
    memcpy(left_edge, le, ndim*sizeof(double));
    memcpy(right_edge, re, ndim*sizeof(double));
    memcpy(periodic_left, ple, ndim*sizeof(bool));
    memcpy(periodic_right, pre, ndim*sizeof(bool));
    for (uint32_t d = 0; d < ndim; d++)
      left_nodes.push_back(left_nodes0[d]);

    left_neighbors = std::vector<std::vector<uint32_t> >(ndim);
    right_neighbors = std::vector<std::vector<uint32_t> >(ndim);

    for (uint32_t d = 0; d < ndim; d++) {
      if ((left_nodes[d]) && (!(left_nodes[d]->is_empty)))
    	add_neighbors(left_nodes[d], d);
    }
  }
  Node(std::istream &is) {
    // Note that Node instances intialized via this method do not have
    // any neighbor information. We will build neighbor information later
    // by walking the tree
    bool check_bit = deserialize_scalar<bool>(is);
    if (!check_bit) {
      // something has gone terribly wrong so we crash
      abort();
    }
    is_empty = deserialize_scalar<bool>(is);
    is_leaf = deserialize_scalar<bool>(is);
    leafid = deserialize_scalar<uint32_t>(is);
    ndim = deserialize_scalar<uint32_t>(is);
    left_edge = deserialize_pointer_array<double>(is, ndim);
    right_edge = deserialize_pointer_array<double>(is, ndim);
    left_idx = deserialize_scalar<uint64_t>(is);
    children = deserialize_scalar<uint64_t>(is);
    periodic_left = deserialize_pointer_array<bool>(is, ndim);
    periodic_right = deserialize_pointer_array<bool>(is, ndim);
    split_dim = deserialize_scalar<uint32_t>(is);
    split = deserialize_scalar<double>(is);
    less = NULL;
    greater = NULL;
    left_neighbors = std::vector<std::vector<uint32_t> >(ndim);
    right_neighbors = std::vector<std::vector<uint32_t> >(ndim);
    for (uint32_t i=0; i<ndim; i++) {
      left_nodes.push_back(NULL);
    }
  }
  void serialize(std::ostream &os) {
    // prepend actual data for Node with true so we can indicate
    // NULL nodes in the data stream, checking with istream.peek()
    serialize_scalar<bool>(os, true);
    serialize_scalar<bool>(os, is_empty);
    serialize_scalar<bool>(os, is_leaf);
    serialize_scalar<uint32_t>(os, leafid);
    serialize_scalar<uint32_t>(os, ndim);
    serialize_pointer_array<double>(os, left_edge, ndim);
    serialize_pointer_array<double>(os, right_edge, ndim);
    serialize_scalar<uint64_t>(os, left_idx);
    serialize_scalar<uint64_t>(os, children);
    serialize_pointer_array<bool>(os, periodic_left, ndim);
    serialize_pointer_array<bool>(os, periodic_right, ndim);
    serialize_scalar<uint32_t>(os, split_dim);
    serialize_scalar<double>(os, split);
  }
  ~Node() {
    if (left_edge)
      free(left_edge);
    if (right_edge)
      free(right_edge);
    if (periodic_left)
      free(periodic_left);
    if (periodic_right)
      free(periodic_right);
  }
  friend std::ostream &operator<<(std::ostream &os, const Node &node) {
    // this is available for nicely formatted debugging, use serialize
    // to save data to disk
    os << "is_empty:      " << node.is_empty << std::endl;
    os << "is_leaf:       " << node.is_leaf << std::endl;
    os << "leafid:        " << node.leafid << std::endl;
    os << "ndim:          " << node.ndim << std::endl;
    os << "left_edge:     ";
    for (uint32_t i = 0; i < node.ndim; i++) {
      os << node.left_edge[i] << " ";
    }
    os << std::endl;
    os << "right_edge:    ";
    for (uint32_t i = 0; i < node.ndim; i++) {
      os << node.right_edge[i] << " ";
    }
    os << std::endl;
    os << "left_idx:      " << node.left_idx << std::endl;
    os << "children:      " << node.children << std::endl;
    os << "periodic_left: ";
    for (uint32_t i = 0; i < node.ndim; i++) {
      os << node.periodic_left[i] << " ";
    }
    os << std::endl;
    os << "periodic_right: ";
    for (uint32_t i = 0; i < node.ndim; i++) {
      os << node.periodic_right[i] << " ";
    }
    os << std::endl;
    os << "split_dim:     " << node.split_dim << std::endl;
    os << "split:         " << node.split << std::endl;
    for (uint32_t i=0; i < node.left_nodes.size(); i++) {
      os << node.left_nodes[i] << std::endl;
      if (node.left_nodes[i]) {
        os << node.left_nodes[i]->left_idx << std::endl;
        os << node.left_nodes[i]->children << std::endl;
      }
    }

    return os;
  }

  Node* copy() {
    Node *out;
    if (is_empty) {
      if (left_edge) {
	out = new Node(ndim, left_edge, right_edge,
		       periodic_left, periodic_right);
      } else {
	out = new Node();
      }
    } else if (is_leaf) {
      std::vector<Node*> left_nodes_copy;
      for (uint32_t d = 0; d < ndim; d++)
	left_nodes_copy.push_back(NULL);
      out = new Node(ndim, left_edge, right_edge,
		     periodic_left, periodic_right,
		     left_idx, children, leafid,
		     left_nodes_copy);
    } else {
      Node *lnode = less->copy();
      Node *gnode = greater->copy();
      std::vector<Node*> left_nodes_copy;
      for (uint32_t d = 0; d < ndim; d++)
	left_nodes_copy.push_back(NULL);
      out = new Node(ndim, left_edge, right_edge,
		     periodic_left, periodic_right,
		     left_idx, split_dim, split, lnode, gnode,
		     left_nodes_copy);
      std::vector<uint32_t>::iterator it;
      for (uint32_t d = 0; d < ndim; d++) {
	for (it = left_neighbors[d].begin();
	     it != left_neighbors[d].end(); it++) {
	  out->left_neighbors[d].push_back(*it);
	}
	for (it = right_neighbors[d].begin();
	     it != right_neighbors[d].end(); it++) {
	  out->right_neighbors[d].push_back(*it);
	}
      }
    }

    return out;
  }

  void update_ids(uint32_t add_to) {
    leafid += add_to;
    uint32_t i;
    for (uint32_t d = 0; d < ndim; d++) {
      for (i = 0; i < left_neighbors[d].size(); i++)
	left_neighbors[d][i] += add_to;
      for (i = 0; i < right_neighbors[d].size(); i++)
	right_neighbors[d][i] += add_to;
    }
    for (i = 0; i < all_neighbors.size(); i++)
      all_neighbors[i] += add_to;
  }

  void print_neighbors() {
    uint32_t i, j;
    // Left
    printf("left:  [");
    for (i = 0; i < ndim; i++) {
      printf("[");
      for (j = 0; j < left_neighbors[i].size(); j++)
	printf("%u ", left_neighbors[i][j]);
      printf("] ");
    }
    printf("]\n");
    // Right
    printf("right: [");
    for (i = 0; i < ndim; i++) {
      printf("[");
      for (j = 0; j < right_neighbors[i].size(); j++)
	printf("%u ", right_neighbors[i][j]);
      printf("] ");
    }
    printf("]\n");
  }

  void add_neighbors(Node* curr, uint32_t dim) {
    if (curr->is_leaf) {
      left_neighbors[dim].push_back(curr->leafid);
      curr->right_neighbors[dim].push_back(leafid);
    } else {
      if (curr->split_dim == dim) {
	add_neighbors(curr->greater, dim);
      } else {
	if (curr->split > this->right_edge[curr->split_dim])
	  add_neighbors(curr->less, dim);
	else if (curr->split < this->left_edge[curr->split_dim])
	  add_neighbors(curr->greater, dim);
	else {
	  add_neighbors(curr->less, dim);
	  add_neighbors(curr->greater, dim);
	}
      }
    }
  }

  void clear_neighbors() {
    uint32_t d;
    for (d = 0; d < ndim; d++) {
      left_neighbors[d].clear();
      right_neighbors[d].clear();
    }
  }

  bool is_left_node(Node *lnode, uint32_t ldim) {
    uint32_t d;
    for (d = 0; d < ndim; d++) {
      if (d == ldim)
	continue;
      if (right_edge[d] < lnode->left_edge[d])
	return false;
      if (left_edge[d] > lnode->right_edge[d])
	return false;
    }
    return true;
  }

  void select_unique_neighbors() {
    if (!is_leaf)
      return;

    uint32_t d;
    std::vector<uint32_t>::iterator last;
    for (d = 0; d < ndim; d++) {
      // left
      std::sort(left_neighbors[d].begin(), left_neighbors[d].end());
      last = std::unique(left_neighbors[d].begin(), left_neighbors[d].end());
      left_neighbors[d].erase(last, left_neighbors[d].end());
      // right
      std::sort(right_neighbors[d].begin(), right_neighbors[d].end());
      last = std::unique(right_neighbors[d].begin(), right_neighbors[d].end());
      right_neighbors[d].erase(last, right_neighbors[d].end());
    }
  }

  void join_neighbors() {
    if (!is_leaf)
      return;

    uint32_t d;
    std::vector<uint32_t>::iterator last;
    // Create concatenated vector and remove duplicates
    all_neighbors = left_neighbors[0];
    for (d = 1; d < ndim; d++)
      all_neighbors.insert(all_neighbors.end(), left_neighbors[d].begin(), left_neighbors[d].end());
    for (d = 0; d < ndim; d++)
      all_neighbors.insert(all_neighbors.end(), right_neighbors[d].begin(), right_neighbors[d].end());

    // Get unique
    std::sort(all_neighbors.begin(), all_neighbors.end());
    last = std::unique(all_neighbors.begin(), all_neighbors.end());
    all_neighbors.erase(last, all_neighbors.end());

  }

  bool check_overlap(Node other, uint32_t dim) {
    if (other.right_edge[dim] < left_edge[dim])
      return false;
    else if (other.left_edge[dim] > right_edge[dim])
      return false;
    else
      return true;
  }

};

void write_tree_nodes(std::ostream &os, Node* node) {
  if (node) {
    // depth first search of tree below node, writing each node to os
    // as we go
    node->serialize(os);
    write_tree_nodes(os, node->less);
    write_tree_nodes(os, node->greater);
  }
  else {
    // write null character to indicate empty node
    serialize_scalar<bool>(os, false);
  }
}

Node* read_tree_nodes(std::istream &is,
                      std::vector<Node*> &leaves,
                      std::vector<Node*> &left_nodes) {
  Node* node = new Node(is);
  node->left_nodes = left_nodes;
  bool is_leaf = true;

  if (is.peek()) {
    // read left subtree
    node->less = read_tree_nodes(is, leaves, left_nodes);
    is_leaf = false;
  }
  else {
    // no left children
    is.get();
    node->less = NULL;
  }

  if (is.peek()) {
    // read right subtree
    std::vector<Node*> greater_left_nodes = left_nodes;
    greater_left_nodes[node->split_dim] = node->less;
    node->greater = read_tree_nodes(is, leaves, greater_left_nodes);
    is_leaf = false;
  }
  else {
    // no right children
    is.get();
    node->greater = NULL;
  }

  if (is_leaf) {
    leaves.push_back(node);
    for (uint32_t d = 0; d < node->ndim; d++) {
      if ((node->left_nodes[d]) && (!(node->left_nodes[d]->is_empty))) {
        node->add_neighbors(node->left_nodes[d], d);
      }
    }
  }

  return node;
}

void free_tree_nodes(Node* node) {
  if (node)
    {
      free_tree_nodes(node->less);
      free_tree_nodes(node->greater);
      delete node;
    }
}

class KDTree
{
public:
  bool is_partial;
  bool skip_dealloc_root;
  bool use_sliding_midpoint;
  uint64_t* all_idx;
  uint64_t npts;
  uint32_t ndim;
  uint64_t left_idx;
  int64_t data_version;
  bool *periodic_left;
  bool *periodic_right;
  uint32_t leafsize;
  double* domain_left_edge;
  double* domain_right_edge;
  double* domain_width;
  bool* periodic;
  bool any_periodic;
  double* domain_mins;
  double* domain_maxs;
  uint32_t num_leaves;
  std::vector<Node*> leaves;
  Node* root;

  // KDTree() {}
  KDTree(double *pts, uint64_t *idx, uint64_t n, uint32_t m,
	 uint32_t leafsize0, double *left_edge, double *right_edge,
	 bool *periodic_left0, bool *periodic_right0,
	 double *domain_mins0, double *domain_maxs0, int64_t dversion,
	 bool use_sliding_midpoint0 = false, bool dont_build = false)
  {
    is_partial = true;
    skip_dealloc_root = false;
    use_sliding_midpoint = use_sliding_midpoint0;

    all_idx = idx;
    npts = n;
    ndim = m;
    leafsize = leafsize0;
    domain_left_edge = (double*)malloc(ndim*sizeof(double));
    domain_right_edge = (double*)malloc(ndim*sizeof(double));
    periodic_left = (bool*)malloc(ndim*sizeof(bool));
    periodic_right = (bool*)malloc(ndim*sizeof(bool));
    data_version = dversion;
    periodic = (bool*)malloc(ndim*sizeof(bool));
    domain_mins = NULL;
    domain_maxs = NULL;
    domain_width = (double*)malloc(ndim*sizeof(double));
    num_leaves = 0;

    memcpy(domain_left_edge, left_edge, ndim*sizeof(double));
    memcpy(domain_right_edge, right_edge, ndim*sizeof(double));
    memcpy(periodic_left, periodic_left0, ndim*sizeof(bool));
    memcpy(periodic_right, periodic_right0, ndim*sizeof(bool));

    if (domain_mins0) {
      domain_mins = (double*)malloc(ndim*sizeof(double));
      memcpy(domain_mins, domain_mins0, ndim*sizeof(double));
    } else if (pts) {
      domain_mins = min_pts(pts, n, m);
    }
    if (domain_maxs0) {
      domain_maxs = (double*)malloc(ndim*sizeof(double));
      memcpy(domain_maxs, domain_maxs0, ndim*sizeof(double));
    } else if (pts) {
      domain_maxs = max_pts(pts, n, m);
    }

    any_periodic = false;
    for (uint32_t d = 0; d < ndim; d++) {
      if ((periodic_left[d]) && (periodic_right[d])) {
	periodic[d] = true;
	any_periodic = true;
      } else {
	periodic[d] = false;
      }
    }

    for (uint32_t d = 0; d < ndim; d++)
      domain_width[d] = domain_right_edge[d] - domain_left_edge[d];

    if ((pts) && (!(dont_build)))
      build_tree(pts);

  }
  KDTree(double *pts, uint64_t *idx, uint64_t n, uint32_t m, uint32_t leafsize0,
	 double *left_edge, double *right_edge, bool *periodic0, int64_t dversion,
	 bool use_sliding_midpoint0 = false, bool dont_build = false)
  {
    is_partial = false;
    skip_dealloc_root = false;
    use_sliding_midpoint = use_sliding_midpoint0;
    left_idx = 0;

    all_idx = idx;
    npts = n;
    ndim = m;
    leafsize = leafsize0;
    domain_left_edge = (double*)malloc(ndim*sizeof(double));
    domain_right_edge = (double*)malloc(ndim*sizeof(double));
    data_version = dversion;
    periodic_left = (bool*)malloc(ndim*sizeof(bool));
    periodic_right = (bool*)malloc(ndim*sizeof(bool));
    periodic = (bool*)malloc(ndim*sizeof(bool));
    domain_mins = NULL;
    domain_maxs = NULL;
    domain_width = (double*)malloc(ndim*sizeof(double));
    num_leaves = 0;

    memcpy(domain_left_edge, left_edge, ndim*sizeof(double));
    memcpy(domain_right_edge, right_edge, ndim*sizeof(double));
    memcpy(periodic, periodic0, ndim*sizeof(bool));

    if (pts) {
      domain_mins = min_pts(pts, n, m);
      domain_maxs = max_pts(pts, n, m);
    }

    any_periodic = false;
    for (uint32_t d = 0; d < ndim; d++) {
      if (periodic[d]) {
	periodic_left[d] = true;
	periodic_right[d] = true;
	any_periodic = true;
      } else {
	periodic_left[d] = false;
	periodic_right[d] = false;
      }
    }

    for (uint32_t d = 0; d < ndim; d++)
      domain_width[d] = domain_right_edge[d] - domain_left_edge[d];

    if ((pts) && (!(dont_build)))
      build_tree(pts);

  }
  KDTree(std::istream &is)
  {
    data_version = deserialize_scalar<int64_t>(is);
    is_partial = deserialize_scalar<bool>(is);
    use_sliding_midpoint = deserialize_scalar<bool>(is);
    npts = deserialize_scalar<uint64_t>(is);
    all_idx = deserialize_pointer_array<uint64_t>(is, npts);
    ndim = deserialize_scalar<uint32_t>(is);
    left_idx = deserialize_scalar<uint64_t>(is);
    periodic = deserialize_pointer_array<bool>(is, ndim);
    periodic_left = deserialize_pointer_array<bool>(is, ndim);
    periodic_right = deserialize_pointer_array<bool>(is, ndim);
    any_periodic = deserialize_scalar<bool>(is);
    leafsize = deserialize_scalar<uint32_t>(is);
    domain_left_edge = deserialize_pointer_array<double>(is, ndim);
    domain_right_edge = deserialize_pointer_array<double>(is, ndim);
    domain_width = deserialize_pointer_array<double>(is, ndim);
    domain_mins = deserialize_pointer_array<double>(is, ndim);
    domain_maxs = deserialize_pointer_array<double>(is, ndim);
    num_leaves = deserialize_scalar<uint32_t>(is);
    std::vector<Node*> left_nodes;
    for (uint32_t i=0; i < ndim; i++) {
      left_nodes.push_back(NULL);
    }
    root = read_tree_nodes(is, leaves, left_nodes);
    finalize_neighbors();
  }
  void serialize(std::ostream &os)
  {
    serialize_scalar<int64_t>(os, data_version);
    serialize_scalar<bool>(os, is_partial);
    serialize_scalar<bool>(os, use_sliding_midpoint);
    serialize_scalar<uint64_t>(os, npts);
    serialize_pointer_array<uint64_t>(os, all_idx, npts);
    serialize_scalar<uint32_t>(os, ndim);
    serialize_scalar<uint64_t>(os, left_idx);
    serialize_pointer_array<bool>(os, periodic, ndim);
    serialize_pointer_array<bool>(os, periodic_left, ndim);
    serialize_pointer_array<bool>(os, periodic_right, ndim);
    serialize_scalar<bool>(os, any_periodic);
    serialize_scalar<uint32_t>(os, leafsize);
    serialize_pointer_array<double>(os, domain_left_edge, ndim);
    serialize_pointer_array<double>(os, domain_right_edge, ndim);
    serialize_pointer_array<double>(os, domain_width, ndim);
    serialize_pointer_array<double>(os, domain_mins, ndim);
    serialize_pointer_array<double>(os, domain_maxs, ndim);
    serialize_scalar<uint32_t>(os, num_leaves);
    write_tree_nodes(os, root);
  }
  ~KDTree()
  {
    if (!(skip_dealloc_root))
      free_tree_nodes(root);
    free(domain_left_edge);
    free(domain_right_edge);
    free(domain_width);
    if (domain_mins)
      free(domain_mins);
    if (domain_maxs)
      free(domain_maxs);
    free(periodic);
    free(periodic_left);
    free(periodic_right);
  }

  void consolidate_edges(double *leaves_le, double *leaves_re) {
    for (uint32_t k = 0; k < num_leaves; k++) {
      memcpy(leaves_le+ndim*leaves[k]->leafid,
             leaves[k]->left_edge,
             ndim*sizeof(double));
      memcpy(leaves_re+ndim*leaves[k]->leafid,
             leaves[k]->right_edge,
             ndim*sizeof(double));
    }
  }

  void build_tree(double* all_pts) {
    uint32_t d;
    double *LE = (double*)malloc(ndim*sizeof(double));
    double *RE = (double*)malloc(ndim*sizeof(double));
    bool *PLE = (bool*)malloc(ndim*sizeof(bool));
    bool *PRE = (bool*)malloc(ndim*sizeof(bool));
    double *mins = (double*)malloc(ndim*sizeof(double));
    double *maxs = (double*)malloc(ndim*sizeof(double));
    std::vector<Node*> left_nodes;

    if (!(domain_mins))
      domain_mins = min_pts(all_pts, npts, ndim);
    if (!(domain_maxs))
      domain_maxs = max_pts(all_pts, npts, ndim);

    for (d = 0; d < ndim; d++) {
      LE[d] = domain_left_edge[d];
      RE[d] = domain_right_edge[d];
      PLE[d] = periodic_left[d];
      PRE[d] = periodic_right[d];
      mins[d] = domain_mins[d];
      maxs[d] = domain_maxs[d];
      left_nodes.push_back(NULL);
    }

    root = build(0, npts, LE, RE, PLE, PRE, all_pts,
                 mins, maxs, left_nodes);

    free(LE);
    free(RE);
    free(PLE);
    free(PRE);
    free(mins);
    free(maxs);

    // Finalize neighbors
    finalize_neighbors();

  }

  void finalize_neighbors() {
    uint32_t d;

    // Add periodic neighbors
    if (any_periodic)
      set_neighbors_periodic();

    // Remove duplicate neighbors
    for (d = 0; d < num_leaves; d++) {
      leaves[d]->select_unique_neighbors();
      leaves[d]->join_neighbors();
    }
  }

  void clear_neighbors() {
    std::vector<Node*>::iterator it;
    for (it = leaves.begin(); it != leaves.end(); it++)
      (*it)->clear_neighbors();
  }

  void set_neighbors_periodic()
  {
    uint32_t d0;
    Node* leaf;
    Node *prev;
    uint64_t i, j;

    // Periodic neighbors
    for (i = 0; i < num_leaves; i++) {
      leaf = leaves[i];
      for (d0 = 0; d0 < ndim; d0++) {
	if (!leaf->periodic_left[d0])
	  continue;
	for (j = i; j < num_leaves; j++) {
	  prev = leaves[j];
	  if (!prev->periodic_right[d0])
	    continue;
	  add_neighbors_periodic(leaf, prev, d0);
	}
      }
    }
  }

  void add_neighbors_periodic(Node *leaf, Node *prev, uint32_t d0) {
    uint32_t d, ndim_escape;
    bool match;
    if (!leaf->periodic_left[d0])
      return;
    if (!prev->periodic_right[d0])
      return;
    match = true;
    ndim_escape = 0;
    for (d = 0; d < ndim; d++) {
      if (d == d0)
	continue;
      if (leaf->left_edge[d] >= prev->right_edge[d]) {
	if (!(leaf->periodic_right[d] && prev->periodic_left[d])) {
	  match = false;
	  break;
	} else {
	  ndim_escape++;
	}
      }
      if (leaf->right_edge[d] <= prev->left_edge[d]) {
	if (!(prev->periodic_right[d] && leaf->periodic_left[d])) {
	  match = false;
	  break;
	} else {
	  ndim_escape++;
	}
      }
    }
    if ((match) && (ndim_escape < (ndim-1))) {
      // printf("%d: %d, %d (%d)\n", d0, leaf->leafid, prev->leafid, ndim_escape);
      leaf->left_neighbors[d0].push_back(prev->leafid);
      prev->right_neighbors[d0].push_back(leaf->leafid);
    }
  }

  Node* build(uint64_t Lidx, uint64_t n,
              double *LE, double *RE,
              bool *PLE, bool *PRE,
              double* all_pts,
              double *mins, double *maxes,
              std::vector<Node*> left_nodes)
  {
    // Create leaf
    if (n < leafsize) {
      Node* out = new Node(ndim, LE, RE, PLE, PRE, Lidx, n, num_leaves,
			   left_nodes);
      num_leaves++;
      leaves.push_back(out);
      return out;
    } else {
      // Split
      uint32_t dmax, d;
      int64_t split_idx = 0;
      double split_val = 0.0;
      dmax = split(all_pts, all_idx, Lidx, n, ndim, mins, maxes,
		   split_idx, split_val, use_sliding_midpoint);
      if (maxes[dmax] == mins[dmax]) {
	// all points singular
	Node* out = new Node(ndim, LE, RE, PLE, PRE, Lidx, n, num_leaves,
			     left_nodes);
	num_leaves++;
	leaves.push_back(out);
	return out;
      }

      // Get new boundaries
      uint64_t Nless = split_idx-Lidx+1;
      uint64_t Ngreater = n - Nless;
      double *lessmaxes = (double*)malloc(ndim*sizeof(double));
      double *lessright = (double*)malloc(ndim*sizeof(double));
      bool *lessPRE = (bool*)malloc(ndim*sizeof(bool));
      double *greatermins = (double*)malloc(ndim*sizeof(double));
      double *greaterleft = (double*)malloc(ndim*sizeof(double));
      bool *greaterPLE = (bool*)malloc(ndim*sizeof(bool));
      std::vector<Node*> greater_left_nodes;
      for (d = 0; d < ndim; d++) {
	lessmaxes[d] = maxes[d];
	lessright[d] = RE[d];
	lessPRE[d] = PRE[d];
	greatermins[d] = mins[d];
	greaterleft[d] = LE[d];
	greaterPLE[d] = PLE[d];
	greater_left_nodes.push_back(left_nodes[d]);
      }
      lessmaxes[dmax] = split_val;
      lessright[dmax] = split_val;
      lessPRE[dmax] = false;
      greatermins[dmax] = split_val;
      greaterleft[dmax] = split_val;
      greaterPLE[dmax] = false;

      // Build less and greater nodes
      Node* less = build(Lidx, Nless, LE, lessright, PLE, lessPRE,
                         all_pts, mins, lessmaxes, left_nodes);
      greater_left_nodes[dmax] = less;
      Node* greater = build(Lidx+Nless, Ngreater, greaterleft, RE,
                            greaterPLE, PRE, all_pts,
                            greatermins, maxes, greater_left_nodes);

      // Create innernode referencing child nodes
      Node* out = new Node(ndim, LE, RE, PLE, PRE, Lidx, dmax, split_val,
			   less, greater, left_nodes);

      free(lessright);
      free(greaterleft);
      free(lessPRE);
      free(greaterPLE);
      free(lessmaxes);
      free(greatermins);
      return out;
    }
  }

  double* wrap_pos(double* pos) {
    uint32_t d;
    double* wrapped_pos = (double*)malloc(ndim*sizeof(double));
    for (d = 0; d < ndim; d++) {
      if (periodic[d]) {
	if (pos[d] < domain_left_edge[d]) {
	  wrapped_pos[d] = domain_right_edge[d] - fmod((domain_right_edge[d] - pos[d]),domain_width[d]);
	} else {
	  wrapped_pos[d] = domain_left_edge[d] + fmod((pos[d] - domain_left_edge[d]),domain_width[d]);
	}
      } else {
	wrapped_pos[d] = pos[d];
      }
    }
    return wrapped_pos;
  }

  Node* search(double* pos0, bool dont_wrap = false)
  {
    uint32_t i;
    Node* out = NULL;
    bool valid;
    // Wrap positions
    double* pos;
    if ((!dont_wrap) && (any_periodic))
      pos = wrap_pos(pos0); // allocates new array
    else
      pos = pos0;
    // Ensure that pos is in root, early return NULL if it's not
    valid = true;
    for (i = 0; i < ndim; i++) {
      if (pos[i] < root->left_edge[i]) {
	valid = false;
	break;
      }
      if (pos[i] >= root->right_edge[i]) {
	valid = false;
	break;
      }
    }
    // Traverse tree looking for leaf containing pos
    if (valid) {
      out = root;
      while (!(out->is_leaf)) {
	if (pos[out->split_dim] < out->split)
	  out = out->less;
	else
	  out = out->greater;
      }
    }

    if ((!dont_wrap) && (any_periodic))
      free(pos);
    return out;
  }

  std::vector<uint32_t> get_neighbor_ids(double* pos)
  {
    Node* leaf;
    std::vector<uint32_t> neighbors;
    leaf = search(pos);
    if (leaf)
      neighbors = leaf->all_neighbors;
    return neighbors;
  }

};
