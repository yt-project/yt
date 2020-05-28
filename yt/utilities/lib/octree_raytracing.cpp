#include <vector>
#include <algorithm>
#include <memory>
#include <array>
#include <iostream>
#include <cassert>
#include <cmath>
#include <random>
#include <string>
#include <fstream>


typedef double F;

/*  A simple node struct that contains a key and a fixed number of children,
    typically Nchildren = 2**Ndim
 */
template <typename keyType, int Nchildren>
struct GenericNode
{
    using _Node = struct GenericNode<keyType, Nchildren>;

    // Tree data
    _Node** children = nullptr;
    _Node* parent = nullptr;

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

template <int Ndim>
struct Ray {
    std::array<F, Ndim> o; // Origin
    std::array<F, Ndim> d; // Direction
    F tmin = -1e99;
    F tmax = 1e99;

    Ray(const std::array<F, Ndim> _o, const std::array<F, Ndim> _d, const F _tmin, const F _tmax) : o(_o), tmin(_tmin), tmax(_tmax) {
        F dd = 0;
        for (auto idim = 0; idim < Ndim; ++idim) {
            dd += _d[idim] * _d[idim];
        }
        dd = std::sqrt(dd);
        for (auto idim = 0; idim < Ndim; ++idim) {
            d[idim] = _d[idim] / dd;
        }
    };

    Ray(const F* _o, const F* _d, const F _tmin, const F _tmax) : tmin(_tmin), tmax(_tmax) {
        for (auto idim = 0; idim < Ndim; ++idim) {
            o[idim] = _o[idim];
        }
        F dd = 0;
        for (auto idim = 0; idim < Ndim; ++idim) {
            dd += _d[idim] * _d[idim];
        }
        dd = std::sqrt(dd);
        for (auto idim = 0; idim < Ndim; ++idim) {
            d[idim] = _d[idim] / dd;
        }
    };

    Ray() {};
};

/*  Converts an array of integer position into a flattened index.
    The fast varying index is the last one.
 */
template <int Ndim>
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
template <int Ndim>
inline std::array<bool, Ndim> iflat2ijk(unsigned char iflat) {
    std::array<bool, Ndim> ijk;
    for (auto idim = Ndim-1; idim >= 0; --idim) {
        ijk[idim] = iflat & 0b1;
        iflat >>= 1;
    }
    return ijk;
};

/*  A class to build an octree and cast rays through it. */
template <class keyType, int Ndim>
class Octree {
    using Node = struct GenericNode<keyType, (1<<Ndim)>;
    using keyVector = std::vector<keyType>;
    using Pos = std::array<F, Ndim>;
    using iPos = std::array<int, Ndim>;
    using ucPos = std::array<unsigned char, Ndim>;

private:
    const unsigned char twotondim;
    const int maxDepth;
    Pos size;
    Pos DLE; // Domain left edge
    Pos DRE; // Domain right edge
    Node* root;
    int global_index = 0;

public:
    Octree(int _maxDepth, F* _size) :
            twotondim (1<<Ndim),
            maxDepth(_maxDepth) {
        for (auto idim = 0; idim < Ndim; ++idim) size[idim] = _size[idim];

        root = new Node();
        // Allocate root's children
        root->children = (Node**) malloc(sizeof(Node*)*twotondim);
        for (auto i = 0; i < twotondim; ++i) root->children = nullptr;

        DLE.fill(0);
        DRE = size;
    }

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
        uint mask = 1<<(maxDepth - 1);

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
            auto iflat = ijk2iflat<Ndim>(bitMask);

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
        Node* n = insert_node(ipos, lvl, key);
    }

    // Perform multiple ray cast
    RayInfo<keyType>** cast_rays(const F *origins, const F *directions, const int Nrays) {
        RayInfo<keyType> **ray_infos = (RayInfo<keyType>**)malloc(sizeof(RayInfo<keyType>*)*Nrays);
        int Nfound = 0;
        #pragma omp parallel for
        for (auto i = 0; i < Nrays; ++i) {
            std::vector<F> tList;
            ray_infos[i] = new RayInfo<keyType>(Nfound);
            auto ri = ray_infos[i];
            Ray<Ndim> r(&origins[3*i], &directions[3*i], -1e99, 1e99);
            cast_ray(&r, ri->keys, ri->t);
            Nfound = std::max(Nfound, (int) ri->keys.size());
        }
        return ray_infos;
    }

    // Perform single ray tracing
    void cast_ray(Ray<Ndim> *r, keyVector &keyList, std::vector<F> &tList) {
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

private:

    /*
        Upsert a node as a child of another.

        This will create a new node as a child of the current one, or return
        an existing one if it already exists
    */
    Node* create_get_node(Node* parent, int iflat) {
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
                free(child);
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
        
        // Early break for leafs without children
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
        } else {                                       // enters XY plane
            if (txm < tz0) index |= 0b100;
            if (tym < tz0) index |= 0b010;
        }
        return index;
    }
    // From "An Efficient Parametric Algorithm for Octree Traversal" by Revelles, Urena, & Lastra
    inline unsigned char next_node(const F tx, const F ty, const F tz,
                                   const u_char ix, const u_char iy, const u_char iz) {
        if(tx < std::min(ty, tz)) {         // YZ plane
            return ix;
        } else if (ty < std::min(tx, tz)) { // XZ plane
            return iy;
        } else {                            // XY plane
            return iz;
        }
    }
};


// Define some instances for easy use in Python
typedef Ray<3> Ray3D;
typedef RayInfo<int> Ray3DInt;
typedef RayInfo<long> Ray3DLong;
template<typename T>
using Octree3D = Octree<T, 3>;

// Instantiate stuff
template class Octree<int, 3>;

void test1() {

    // std::array<bool, 3> bitMask = {true, false, true};

    for (unsigned char i = 0; i < 8; i++){
        // auto tmp = iflat2ijk<3>(i);
        // std::cout << (int)i << " -> " << tmp[0] << tmp[1] << tmp[2] << " -> " << (int)ijk2iflat<3>(iflat2ijk<3>(i)) << std::endl;
        assert(ijk2iflat<3>(iflat2ijk<3>(i)) == i);
    }

}

void test2() {
    // Shouldnt crash
    int index = 0;
    int N = 4;
    F size[3] = {1, 1, 1};
    Octree<std::array<int, 3>, 3> o(N, size);
    for (auto i = 0; i < 1<<N; ++i) {
        for (auto j = 0; j < 1<<N; ++j) {
            for (auto k = 0; k < 1<<N; ++k) {
                o.insert_node({i, j, k}, N, std::array<int, 3>({i, j, k}));
                ++index;
            }
        }
    }
}

void test3(){
    int N = 4;
    F size[3] = {1, 1, 1};
    Octree<std::array<int, 3>, 3> o(N, size);
    F ox, oy, oz;
    F rx, ry, rz;
    // std::cin >> ox >> oy >> oz;
    // std::cin >> rx >> ry >> rz;

    ox = 0.01;
    oy = 0.84;
    oz = 0.95;

    rx = 1.;
    ry = -1.2;
    rz = -1.5;

    F oo[3];
    oo[0] = ox;
    oo[1] = oy;
    oo[2] = oz;
    F rr[3];
    rr[0] = rx;
    rr[1] = ry;
    rr[2] = rz;
    Ray<3> r(oo, rr, -1e99, 1e99);
    std::cerr<< "Casting ray in direction:\t" << rx << ", " << ry << ", " << rz << "(len=" << std::sqrt(r.d[0]*r.d[0] + r.d[1]*r.d[1] + r.d[2]*r.d[2]) << ")" <<std::endl;
    std::cerr<< "Origin:\t" << ", " << ox << ", " << oy << ", " << oz << std::endl;

    std::vector<std::array<int, 3>> ret;
    std::vector<double> tList;
    // o.cast_ray(&r, ret, tList);
}


void test4() { 
    int N = 6;
    F size[3] = {1, 1, 1};
    Octree<std::array<int, 3>, 3> o(N, size);

    // Filling half of octree at level 3
    for (auto i = 0; i < 1<<N; ++i) {
        for (auto j = 0; j < 1<<N; ++j) {
            for (auto k = 0; k < 1<<N; ++k) {
                if (i < 1<<(N-1)) {
                    o.insert_node({i, j, k}, N, {i, j, k});
                } else if (i % 2) {
                    o.insert_node({i, j, k}, N-1, {i, j, k});
                }
            }
        }
    }

    // Now trying to cast ray
    std::vector<std::array<F, 3>> pos(1024*1024);
    std::vector<std::array<F, 3>> dir(1024*1024);

    std::mt19937 gen(16091992);
    auto dis = std::uniform_real_distribution<> (0., 1.);
    for (auto i = 0; i < (int) pos.size(); ++i) {
        for (auto idim = 0; idim < 3; ++idim) {
            pos[i][idim] = dis(gen);
            dir[i][idim] = dis(gen)*2-1;
        }
    }
    // auto ret = o.cast_rays(pos, dir, (int) pos.size());
    // for (auto k: ret) {
    //     std::cout << k[0] << " " << k[1] << " " << k[2] << std::endl;
    // }
}

void test5() {
    std::ifstream inFile;

    inFile.open("/tmp/ipos.txt");
    int ix, iy, iz;
    std::vector<std::array<int, 3>> ipos;
    while (inFile >> ix >> iy >> iz) {
        ipos.push_back({ix, iy, iz});
    }

    inFile.close();

    inFile.open("/tmp/lvl.txt");
    std::vector<int> ilvl;
    while (inFile >> ix) {
        ilvl.push_back(ix);
    }

    // Create octree
    double size[3] = {1, 1, 1};
    Octree3D<int> oct(16, size);
    for (auto i = 0; i < (int) ipos.size(); ++i) {
        oct.insert_node(&ipos[i][0], ilvl[i], i);
    }

    // Now cast a ray
    double o[3] = {0.5, 0.5, 0.5};
    double d[3] = {1., 2., 3.};
    Ray3D r(o, d, -1e99, 1e99);
    std::vector<int> keyList;
    std::vector<double> tList;
    oct.cast_ray(&r, keyList, tList);
}

int main() {
    // std::cout << "########################## TEST 1 ##########################" << std::endl;
    // test1();
    // std::cout << "########################## TEST 2 ##########################" << std::endl;
    // test2();
    // std::cout << "########################## TEST 3 ##########################" << std::endl;
    // test3();
    // std::cout << "########################## TEST 4 ##########################" << std::endl;
    // test4();
    std::cout << "########################## TEST 5 ##########################" << std::endl;
    test5();
    return 0;
}
