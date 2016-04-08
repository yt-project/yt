"""
Wrapper for EWAH Bool Array: https://github.com/lemire/EWAHBoolArray



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

cimport numpy as np
cimport cython
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string

# Streams req for c++ IO
cdef extern from "<ostream>" namespace "std":
    cdef cppclass ostream[T]:
        pass
cdef extern from "<istream>" namespace "std":
    cdef cppclass istream[T]:
        pass

cdef extern from "<sstream>" namespace "std":
    cdef cppclass stringstream:
        stringstream() except +
        string str()
        ostream write(char *, size_t)
        istream read(char *, size_t)
        bint eof()

cdef extern from "ewah.h":
    cppclass EWAHBoolArraySetBitForwardIterator[uword]:
        # EWAHBoolArraySetBitForwardIterator()
        EWAHBoolArraySetBitForwardIterator(const EWAHBoolArraySetBitForwardIterator &o)
        size_t operator*()
        EWAHBoolArraySetBitForwardIterator &operator++()
        bint operator==(EWAHBoolArraySetBitForwardIterator &x)
        bint operator!=(EWAHBoolArraySetBitForwardIterator &x)
    # ctypedef EWAHBoolArraySetBitForwardIterator[unsigned long long] const_iterator
    cdef cppclass EWAHBoolArray[uword]:
        # We are going to skip the varargs here; it is too tricky to assemble.
        bint get(const size_t pos)
        bint set(size_t i)
        void makeSameSize(EWAHBoolArray &a)
        vector[size_t] toArray()
        void logicaland(EWAHBoolArray &a, EWAHBoolArray &container)
        void logicalor(EWAHBoolArray &a, EWAHBoolArray &container)
        void logicalxor(EWAHBoolArray &a, EWAHBoolArray &container)
        bint intersects(EWAHBoolArray &a)
        void reset()
        size_t sizeInBits()
        size_t sizeInBytes()
        bint operator==(EWAHBoolArray &x)
        bint operator!=(EWAHBoolArray &x)
        void append(EWAHBoolArray &x)
        # Recommended container is "vector[size_t]"
        void appendRowIDs[container](container &out, const size_t offset)
        void appendSetBits[container](container &out, const size_t offset)
        size_t numberOfOnes()
        void logicalnot(EWAHBoolArray &x)
        void inplace_logicalnot()
        void swap(EWAHBoolArray &x)
        void read(stringstream &incoming, bint savesizeinbits)
        void readBuffer(stringstream &incoming, const size_t buffersize)
        void write(stringstream &out, bint savesizeinbits)
        void writeBuffer(stringstream &out)
        vector[uword] &getBuffer()
        # const_iterator begin()
        # const_iterator end()
        EWAHBoolArraySetBitForwardIterator begin()
        EWAHBoolArraySetBitForwardIterator end()

ctypedef EWAHBoolArray[np.uint64_t] ewah_bool_array
#ctypedef EWAHBoolArray[np.uint64_t].EWAHBoolArraySetBitForwardIterator ewah_bool_iterator
ctypedef EWAHBoolArraySetBitForwardIterator[np.uint64_t] ewah_bool_iterator
ctypedef vector[size_t] bitset_array
ctypedef map[np.uint64_t, ewah_bool_array] ewah_map
ctypedef stringstream sstream