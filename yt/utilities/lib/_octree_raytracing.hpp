#include <vector>
#include <algorithm>
#include <memory>
#include <array>
#include <iostream>
#include <cassert>
#include <math.h>
#include <random>
#include <string>
#include <fstream>


typedef double F;
const int Ndim = 3;

/*  A simple node struct that contains a key and a fixed number of children,
    typically Nchildren = 2**Ndim
 */
template <typename keyType, int Nchildren>
struct GenericNode
{
    // Tree data
    GenericNode<keyType, Nchildren>** children = nullptr;
    GenericNode<keyType, Nchildren>* parent = nullptr;

    // Node data
    keyType key;
    int level = 0;
    bool terminal = false;
    int index = -1;
};

template <typename keyType>
struct RayInfo {
    std::vector<keyType> keys;
    std::vector<F> t;

    RayInfo() {};
    RayInfo(int N) {
        if (N > 0) {
            keys.reserve(N);
            t.reserve(2*N);
        }
    }
};

struct Ray {
    std::array<F, Ndim> o; // Origin
    std::array<F, Ndim> d; // Direction
    F tmin = -1e99;
    F tmax = 1e99;

    Ray(const F* _o, const F* _d, const F _tmin, const F _tmax) : tmin(_tmin), tmax(_tmax) {
        for (auto idim = 0; idim < Ndim; ++idim) {
            o[idim] = _o[idim];
        }
        F dd = 0;
        for (auto idim = 0; idim < Ndim; ++idim) {
            dd += _d[idim] * _d[idim];
        }
        dd = sqrt(dd);
        for (auto idim = 0; idim < Ndim; ++idim) {
            d[idim] = _d[idim] / dd;
        }
    };

    Ray() {};
};

/*  Converts an array of integer position into a flattened index.
    The fast varying index is the last one.
 */
inline unsigned char ijk2iflat(const std::array<bool, Ndim> ijk) {
    unsigned char iflat = 0;
    for (auto i : ijk) {
        iflat += i;
        iflat <<= 1;
    }
    return iflat >> 1;
};

/*  Converts a flattened index into an array of integer position.
    The fast varying index is the last one.
*/
inline std::array<bool, Ndim> iflat2ijk(unsigned char iflat) {
    std::array<bool, Ndim> ijk;
    for (auto idim = Ndim-1; idim >= 0; --idim) {
        ijk[idim] = iflat & 0b1;
        iflat >>= 1;
    }
    return ijk;
};

/*  A class to build an octree and cast rays through it. */
template <class keyType>
class Octree {
    typedef GenericNode<keyType, (1<<Ndim)> Node;
    typedef std::vector<keyType> keyVector;
    typedef std::array<F, Ndim> Pos;
    typedef std::array<int, Ndim> iPos;
    typedef std::array<unsigned char, Ndim> ucPos;

private:
    const unsigned char twotondim;
    const int maxDepth;
    Pos size;
    Pos DLE; // Domain left edge
    Pos DRE; // Domain right edge
    Node* root;
    int global_index = 0;

public:
    Octree(int _maxDepth, F* _DLE, F* _DRE) :
            twotondim (1<<Ndim),
            maxDepth(_maxDepth) {
        for (auto idim = 0; idim < Ndim; ++idim) {
            DRE[idim] = _DRE[idim];
            DLE[idim] = _DLE[idim];
            size[idim] = _DRE[idim] - _DLE[idim];
        }

        root = new Node();
        // Allocate root's children
        root->children = (Node**) malloc(sizeof(Node*)*twotondim);
        for (auto i = 0; i < twotondim; ++i) root->children = nullptr;
    }


    ~Octree() {
        recursive_remove_node(root);
    };

    /*
        Insert a new node in the tree.
    */
    Node* insert_node(const iPos ipos, const int lvl, keyType key) {
        assert(lvl <= maxDepth);

        // std::cerr << "Inserting at level: " << lvl << "/" << maxDepth << std::endl;
        // this is 0b100..., where the 1 is at position maxDepth
        uint64_t mask = 1<<(maxDepth - 1);

        iPos ijk = ipos;
        std::array<bool, Ndim> bitMask;

        Node* node = root;
        Node* child = nullptr;

        // Go down the tree
        for (auto ibit = maxDepth-1; ibit >= maxDepth - lvl; --ibit) {
            // Find children based on bits
            for (auto idim = 0; idim < Ndim; ++idim) {
                bitMask[idim] = ijk[idim] & mask;
            }
            mask >>= 1;
            auto iflat = ijk2iflat(bitMask);

            // Create child if it does not exist yet
            child = create_get_node(node, iflat);
            node = child;
        }

        // Mark last node as terminal
        node->terminal = true;
        node->key = key;

        return node;
    }

    Node* insert_node(const int* ipos, const int lvl, keyType key) {
        std::array<int, Ndim> ipos_as_arr;
        for (auto idim = 0; idim < Ndim; ++idim) ipos_as_arr[idim] = ipos[idim];
        return insert_node(ipos_as_arr, lvl, key);
    }

    void insert_node_no_ret(const int* ipos, const int lvl, keyType key) {
        insert_node(ipos, lvl, key);
    }

    // Perform multiple ray cast
    RayInfo<keyType>** cast_rays(const F *origins, const F *directions, const int Nrays) {
        RayInfo<keyType> **ray_infos = (RayInfo<keyType>**) malloc(sizeof(RayInfo<keyType>*)*Nrays);
        int Nfound = 0;
        #pragma omp parallel for shared(ray_infos, Nfound) schedule(static, 100)
        for (auto i = 0; i < Nrays; ++i) {
            std::vector<F> tList;
            ray_infos[i] = new RayInfo<keyType>(Nfound);
            auto ri = ray_infos[i];
            Ray r(&origins[3*i], &directions[3*i], -1e99, 1e99);
            cast_ray(&r, ri->keys, ri->t);

            // Keep track of the number of cells hit to preallocate the next ray info container
            Nfound = std::max(Nfound, (int) ri->keys.size());
        }
        return ray_infos;
    }

    // Perform single ray tracing
    void cast_ray(Ray *r, keyVector &keyList, std::vector<F> &tList) {
        // Boolean mask for direction
        unsigned char a = 0;
        unsigned char bmask = twotondim >> 1;

        // Put ray in positive direction and store info in bitmask "a"
        for (auto idim = 0; idim < Ndim; ++idim) {
            if (r->d[idim] < 0.0) {
                r->o[idim] = size[idim]-r->o[idim];
                r->d[idim] = -r->d[idim];
                a |= bmask;
            }
            bmask >>= 1;
        }

        // Compute intersection points
        Pos t0, t1;
        for (auto idim = 0; idim < Ndim; ++idim){
            t0[idim] = (DLE[idim] - r->o[idim]) / r->d[idim];
            t1[idim] = (DRE[idim] - r->o[idim]) / r->d[idim];
        }

        // If entry point is smaller than exit point, find path in octree
        if (*std::max_element(t0.begin(), t0.end()) < *std::min_element(t1.begin(), t1.end()))
            proc_subtree(t0[0], t0[1], t0[2],
                         t1[0], t1[1], t1[2],
                         root, a, keyList, tList);
    }

    void cast_ray(double* o, double* d, keyVector &keyList, std::vector<F> &tList) {
        Ray r(o, d, -1e99, 1e99);
        cast_ray(&r, keyList, tList);
    }

private:

    /*
        Upsert a node as a child of another.

        This will create a new node as a child of the current one, or return
        an existing one if it already exists
    */
    Node* create_get_node(Node* parent, const int iflat) {
        // Create children if not already existing
        if (parent->children == nullptr) {
            parent->children = (Node**) malloc(sizeof(Node*)*twotondim);
            for (auto i = 0; i < twotondim; ++i) parent->children[i] = nullptr;
        }

        if (parent->children[iflat] == nullptr) {
            Node* node = new Node();
            node->level = parent->level + 1;
            node->index = global_index;
            node->parent = parent;
            ++global_index;

            parent->children[iflat] = node;
        }
        return parent->children[iflat];
    }

    /*
        Recursively free memory.
    */
    void recursive_remove_node(Node* node) {
        if (node->children) {
            for (auto i = 0; i < twotondim; ++i) {
                auto child = node->children[i];
                if (child) {
                    recursive_remove_node(child);
                }
                delete child;
            }
            free(node->children);
        }
    }

    /*
        Traverse the tree, assuming that the ray intersects
        From http://wscg.zcu.cz/wscg2000/Papers_2000/X31.pdf
    */
    void proc_subtree(const F tx0, const F ty0, const F tz0,
                      const F tx1, const F ty1, const F tz1,
                      const Node *n, const unsigned char a,
                      keyVector &keyList, std::vector<F> &tList, int lvl=0) {
        // Check if exit face is not in our back
        if (tx1 < 0 || ty1 < 0 || tz1 < 0) return;

        // Exit if the node is null (happens if it hasn't been added to the tree)
        if (!n) return;

        // Process leaf node
        if (n->terminal) {
            keyList.push_back(n->key);
            // Push entry & exit t
            tList.push_back(std::max(std::max(tx0, ty0), tz0));
            tList.push_back(std::min(std::min(tx1, ty1), tz1));
            assert(n->children == nullptr);
            return;
        }

        // Early break for leafs without children that aren't terminal
        if (n->children == nullptr) return;

        // Compute middle intersection
        F txm, tym, tzm;
        txm = (tx0 + tx1) * 0.5;
        tym = (ty0 + ty1) * 0.5;
        tzm = (tz0 + tz1) * 0.5;

        unsigned char iNode = first_node(tx0, ty0, tz0, txm, tym, tzm);

        // Iterate over children
        do {

            switch (iNode)
            {
            case 0:
                proc_subtree(tx0, ty0, tz0, txm, tym, tzm, n->children[a], a, keyList, tList, lvl+1);
                iNode = next_node(txm, tym, tzm, 4, 2, 1);
                break;
            case 1:
                proc_subtree(tx0, ty0, tzm, txm, tym, tz1, n->children[1^a], a, keyList, tList, lvl+1);
                iNode = next_node(txm, tym, tz1, 5, 3, 8);
                break;
            case 2:
                proc_subtree(tx0, tym, tz0, txm, ty1, tzm, n->children[2^a], a, keyList, tList, lvl+1);
                iNode = next_node(txm, ty1, tzm, 6, 8, 3);
                break;
            case 3:
                proc_subtree(tx0, tym, tzm, txm, ty1, tz1, n->children[3^a], a, keyList, tList, lvl+1);
                iNode = next_node(txm, ty1, tz1, 7, 8, 8);
                break;
            case 4:
                proc_subtree(txm, ty0, tz0, tx1, tym, tzm, n->children[4^a], a, keyList, tList, lvl+1);
                iNode = next_node(tx1, tym, tzm, 8, 6, 5);
                break;
            case 5:
                proc_subtree(txm, ty0, tzm, tx1, tym, tz1, n->children[5^a], a, keyList, tList, lvl+1);
                iNode = next_node(tx1, tym, tz1, 8, 7, 8);
                break;
            case 6:
                proc_subtree(txm, tym, tz0, tx1, ty1, tzm, n->children[6^a], a, keyList, tList, lvl+1);
                iNode = next_node(tx1, ty1, tzm, 8, 8, 7);
                break;
            case 7:
                proc_subtree(txm, tym, tzm, tx1, ty1, tz1, n->children[7^a], a, keyList, tList, lvl+1);
                iNode = 8;
                break;
            }
        } while (iNode < twotondim);
    }

    // From "An Efficient Parametric Algorithm for Octree Traversal" by Revelles, Urena, & Lastra
    inline unsigned char first_node(const F tx0, const F ty0, const F tz0,
                                    const F txm, const F tym, const F tzm) {
        unsigned char index = 0;
        if (tx0 >= std::max(ty0, tz0)) {         // enters YZ plane
            if (tym < tx0) index |= 0b010;
            if (tzm < tx0) index |= 0b001;
        } else if (ty0 >= std::max(tx0, tz0))  { // enters XZ plane
            if (txm < ty0) index |= 0b100;
            if (tzm < ty0) index |= 0b001;
        } else {                                 // enters XY plane
            if (txm < tz0) index |= 0b100;
            if (tym < tz0) index |= 0b010;
        }
        return index;
    }
    // From "An Efficient Parametric Algorithm for Octree Traversal" by Revelles, Urena, & Lastra
    inline unsigned char next_node(const F tx, const F ty, const F tz,
                                   const uint8_t ix, const uint8_t iy, const uint8_t iz) {
        if(tx < std::min(ty, tz)) {         // YZ plane
            return ix;
        } else if (ty < std::min(tx, tz)) { // XZ plane
            return iy;
        } else {                            // XY plane
            return iz;
        }
    }
};
