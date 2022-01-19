"""
Wrapper for EWAH Bool Array: https://github.com/lemire/EWAHBoolArray



"""


from libc.stdint cimport uint32_t, uint64_t
from libcpp cimport bool
from libcpp.map cimport map as cmap
from libcpp.string cimport string
from libcpp.vector cimport vector


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

cdef extern from "ewah.h" namespace "ewah":
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
        size_t addWord(uword newdata)
        vector[uword] &getBuffer()
        # const_iterator begin()
        # const_iterator end()
        EWAHBoolArraySetBitForwardIterator begin()
        EWAHBoolArraySetBitForwardIterator end()

cdef extern from "boolarray.h" namespace "ewah":
    cppclass BoolArray[uword]:
        void setSizeInBits(size_t sizeib)
        void set(size_t pos)
        void unset(size_t pos)
        bool get(size_t pos)
        void reset()
        size_t sizeInBits()
        size_t sizeInBytes()
        size_t numberOfOnes()
        void inplace_logicalxor(BoolArray &other)
        void inplace_logicalnot()
        size_t padWithZeroes(size_t totalbits)
        uword getWord(size_t pos)
        size_t wordinbits

cimport cython
cimport numpy as np

IF UNAME_SYSNAME == "Windows":
    ctypedef uint32_t ewah_word_type
ELSE:
    ctypedef np.uint32_t ewah_word_type
ctypedef EWAHBoolArray[ewah_word_type] ewah_bool_array
ctypedef EWAHBoolArraySetBitForwardIterator[ewah_word_type] ewah_bool_iterator
ctypedef vector[size_t] bitset_array
ctypedef cmap[np.uint64_t, ewah_bool_array] ewah_map
ctypedef stringstream sstream
ctypedef BoolArray[ewah_word_type] bool_array
