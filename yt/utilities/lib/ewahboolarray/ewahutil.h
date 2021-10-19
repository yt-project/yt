/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 *
 * Some code from the public domain tuklib.
 */

#ifndef EWAHUTIL_H
#define EWAHUTIL_H

#include <iso646.h> // mostly for Microsoft compilers
#include <limits.h>
#include <stdint.h> // part of Visual Studio 2010 and better
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef _MSC_VER
#include <intrin.h>
#endif


#if ((ULONG_MAX) == (UINT_MAX))
#define UWORD uint32_t
#else
#define UWORD uint64_t
#endif

namespace ewah {

static inline uint32_t ctz64(uint64_t n) {
#if defined(__GNUC__) && UINT_MAX >= UINT32_MAX && ULLONG_MAX >= UINT64_MAX
  return static_cast<uint32_t>(__builtin_ctzll(n));
#elif defined(_WIN64) && defined(_MSC_VER) && _MSC_VER >= 1400 &&              \
    ULONG_MAX >= UINT64_MAX
  uint32_t i;
  _BitScanForward64((unsigned long *)&i, n);
  return i;
#else
  uint32_t i = 1;
  if ((n & static_cast<uint64_t>(4294967295)) == 0) {
    n >>= 32;
    i += 32;
  }
  if ((n & static_cast<uint64_t>(0x0000FFFFUL)) == 0) {
    n >>= 16;
    i += 16;
  }

  if ((n & static_cast<uint64_t>(0x000000FFUL)) == 0) {
    n >>= 8;
    i += 8;
  }

  if ((n & static_cast<uint64_t>(0x0000000FUL)) == 0) {
    n >>= 4;
    i += 4;
  }

  if ((n & static_cast<uint64_t>(0x00000003UL)) == 0) {
    n >>= 2;
    i += 2;
  }
  i -= (n & 0x1);
  return i;
#endif
}

static inline uint32_t ctz32(uint32_t n) {
#if defined(__GNUC__) && UINT_MAX >= UINT32_MAX
  return static_cast<uint32_t>(__builtin_ctz(n));

#elif defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
  uint32_t i;
  __asm__("bsfl %1, %0" : "=r"(i) : "rm"(n));
  return i;

#elif defined(_MSC_VER) && _MSC_VER >= 1400
  uint32_t i;
  _BitScanForward((unsigned long *)&i, n);
  return i;

#else
  uint32_t i = 1;

  if ((n & static_cast<uint32_t>(0x0000FFFF)) == 0) {
    n >>= 16;
    i += 16;
  }

  if ((n & static_cast<uint32_t>(0x000000FF)) == 0) {
    n >>= 8;
    i += 8;
  }

  if ((n & static_cast<uint32_t>(0x0000000F)) == 0) {
    n >>= 4;
    i += 4;
  }

  if ((n & static_cast<uint32_t>(0x00000003)) == 0) {
    n >>= 2;
    i += 2;
  }

  i -= (n & 1);

  return i;
#endif
}

static inline uint32_t ctz16(uint16_t n) {
#if defined(__GNUC__) && UINT_MAX >= UINT32_MAX
  return static_cast<uint32_t>(__builtin_ctz(n));

#elif defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
  uint32_t i;
  __asm__("bsfl %1, %0" : "=r"(i) : "rm"(n));
  return i;

#elif defined(_MSC_VER) && _MSC_VER >= 1400
  uint32_t i;
  _BitScanForward((unsigned long *)&i, n);
  return i;

#else
  uint32_t i = 1;

  if ((n & static_cast<uint16_t>(0x000000FF)) == 0) {
    n >>= 8;
    i += 8;
  }

  if ((n & static_cast<uint16_t>(0x0000000F)) == 0) {
    n >>= 4;
    i += 4;
  }

  if ((n & static_cast<uint16_t>(0x00000003)) == 0) {
    n >>= 2;
    i += 2;
  }
  i -= (n & 1);

  return i;
#endif
}

#ifdef __GNUC__
/**
 * count the number of bits set to one (32 bit version)
 */
inline uint32_t countOnes(uint32_t x) {
  return static_cast<uint32_t>(__builtin_popcount(x));
}
#elif defined(_MSC_VER) && _MSC_VER >= 1400 && !defined(_M_ARM)&& !defined(_M_ARM64)
inline uint32_t countOnes(uint32_t x) { return __popcnt(x); }
#else
inline uint32_t countOnes(uint32_t v) {
  v = v - ((v >> 1) & 0x55555555);
  v = (v & 0x33333333) + ((v >> 2) & 0x33333333);
  return static_cast<uint32_t>((((v + (v >> 4)) & 0x0F0F0F0F) * 0x01010101) >>
                               24);
}
#endif

#ifdef __GNUC__
/**
 * count the number of bits set to one (64 bit version)
 */
inline uint32_t countOnes(uint64_t x) {
  return static_cast<uint32_t>(__builtin_popcountll(x));
}
#elif defined(_WIN64) && defined(_MSC_VER) && _MSC_VER >= 1400 && !defined(_M_ARM64)
inline uint32_t countOnes(uint64_t x) {
  return static_cast<uint32_t>(__popcnt64(static_cast<__int64>(x)));
}
#else
inline uint32_t countOnes(uint64_t v) {
  v = v - ((v >> 1) & 0x5555555555555555);
  v = (v & 0x3333333333333333) + ((v >> 2) & 0x3333333333333333);
  v = ((v + (v >> 4)) & 0x0F0F0F0F0F0F0F0F);
  return static_cast<uint32_t>((v * (0x0101010101010101)) >> 56);
}
#endif

inline uint32_t countOnes(uint16_t v) {
  return countOnes(static_cast<uint32_t>(v));
}

inline uint32_t numberOfTrailingZeros(uint32_t x) {
  if (x == 0)
    return 32;
  return ctz32(x);
}

inline uint32_t numberOfTrailingZeros(uint64_t x) {
  if (x == 0)
    return 64;
  return ctz64(x);
}

inline uint32_t numberOfTrailingZeros(uint16_t x) {
  if (x == 0)
    return 16;
  return ctz16(x);
}

/**
 * Returns the binary representation of a binary word.
 */
template <class uword> std::string toBinaryString(const uword w) {
  std::ostringstream convert;
  for (uint32_t k = 0; k < sizeof(uword) * 8; ++k) {
    if (w & (static_cast<uword>(1) << k))
      convert << "1";
    else
      convert << "0";
  }
  return convert.str();
}
} // namespace ewah
#endif
