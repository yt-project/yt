/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 *     with contributions from Zarian Waheed and others.
 */

#ifndef EWAH_H
#define EWAH_H

#include <algorithm>
#include <queue>
#include <vector>

#include "boolarray.h"
#include "ewahutil.h"

#include "runninglengthword.h"

namespace ewah {

template <class uword> class EWAHBoolArrayIterator;

template <class uword> class EWAHBoolArraySetBitForwardIterator;

class BitmapStatistics;

template <class uword> class EWAHBoolArrayRawIterator;

/**
 * This class is a compressed bitmap.
 * This is where compression
 * happens.
 * The underlying data structure is an STL vector.
 */
template <class uword = uint32_t> class EWAHBoolArray {
public:
  EWAHBoolArray() : buffer(1, 0), sizeinbits(0), lastRLW(0) {}

  static EWAHBoolArray bitmapOf(size_t n, ...) {
    EWAHBoolArray ans;
    va_list vl;
    va_start(vl, n);
    for (size_t i = 0; i < n; i++) {
      ans.set(static_cast<size_t>(va_arg(vl, int)));
    }
    va_end(vl);
    return ans;
  }

  /**
   * Recover wasted memory usage. Fit buffers to the actual data.
   */
  void trim() { buffer.shrink_to_fit(); }

  /**
   * Query the value of bit i. This runs in time proportional to
   * the size of the bitmap. This is not meant to be use in
   * a performance-sensitive context.
   *
   *  (This implementation is based on zhenjl's Go version of JavaEWAH.)
   *
   */
  bool get(const size_t pos) const {
    if (pos >= static_cast<size_t>(sizeinbits))
      return false;
    const size_t wordpos = pos / wordinbits;
    size_t WordChecked = 0;
    EWAHBoolArrayRawIterator<uword> j = raw_iterator();
    while (j.hasNext()) {
      BufferedRunningLengthWord<uword> &rle = j.next();
      WordChecked += static_cast<size_t>(rle.getRunningLength());
      if (wordpos < WordChecked)
        return rle.getRunningBit();
      if (wordpos < WordChecked + rle.getNumberOfLiteralWords()) {
        const uword w = j.dirtyWords()[wordpos - WordChecked];
        return (w & (static_cast<uword>(1) << (pos % wordinbits))) != 0;
      }
      WordChecked += static_cast<size_t>(rle.getNumberOfLiteralWords());
    }
    return false;
  }

  /**
   * Returns true if no bit is set.
   */
  bool empty() const {
    size_t pointer(0);
    while (pointer < buffer.size()) {
      ConstRunningLengthWord<uword> rlw(buffer[pointer]);
      if (rlw.getRunningBit()) {
        if (rlw.getRunningLength() > 0)
          return false;
      }
      ++pointer;
      for (size_t k = 0; k < rlw.getNumberOfLiteralWords(); ++k) {
        if (buffer[pointer] != 0)
          return false;
        ++pointer;
      }
    }
    return true;
  }

  /**
   * Set the ith bit to true (starting at zero).
   * Auto-expands the bitmap. It has constant running time complexity.
   * Note that you must set the bits in increasing order:
   * set(1), set(2) is ok; set(2), set(1) is not ok.
   * set(100), set(100) is also not ok.
   *
   * Note: by design EWAH is not an updatable data structure in
   * the sense that once bit 1000 is set, you cannot change the value
   * of bits 0 to 1000.
   *
   * Returns true if the value of the bit was changed, and false otherwise.
   * (In practice, if you set the bits in strictly increasing order, it
   * should always return true.)
   */
  bool set(size_t i);

  /**
   * Transform into a string that presents a list of set bits.
   * The running time is linear in the compressed size of the bitmap.
   */
  operator std::string() const {
    std::stringstream ss;
    ss << *this;
    return ss.str();
  }
  friend std::ostream &operator<<(std::ostream &out, const EWAHBoolArray &a) {

    out << "{";
    for (EWAHBoolArray::const_iterator i = a.begin(); i != a.end();) {
      out << *i;
      ++i;
      if (i != a.end())
        out << ",";
    }
    out << "}";

    return out;
  }
  /**
   * Make sure the two bitmaps have the same size (padding with zeroes
   * if necessary). It has constant running time complexity.
   *
   * This is useful when calling "logicalnot" functions.
   *
   * This can an adverse effect of performance, especially when computing
   * intersections.
   */
  void makeSameSize(EWAHBoolArray &a) {
    if (a.sizeinbits < sizeinbits)
      a.padWithZeroes(sizeinbits);
    else if (sizeinbits < a.sizeinbits)
      padWithZeroes(a.sizeinbits);
  }

  enum { RESERVEMEMORY = true }; // for speed

  typedef EWAHBoolArraySetBitForwardIterator<uword> const_iterator;

  /**
   * Returns an iterator that can be used to access the position of the
   * set bits. The running time complexity of a full scan is proportional to the
   * number
   * of set bits: be aware that if you have long strings of 1s, this can be
   * very inefficient.
   *
   * It can be much faster to use the toArray method if you want to
   * retrieve the set bits.
   */
  const_iterator begin() const {
    return EWAHBoolArraySetBitForwardIterator<uword>(&buffer);
  }

  /**
   * Basically a bogus iterator that can be used together with begin()
   * for constructions such as for(EWAHBoolArray<uword>::iterator i = b.begin();
   * i!=b.end(); ++i) {}
   */
  const_iterator &end() const {
    return EWAHBoolArraySetBitForwardIterator<uword>::end();
  }

  /**
   * Retrieve the set bits. Can be much faster than iterating through
   * the set bits with an iterator.
   */
  std::vector<size_t> toArray() const;

  /**
   * computes the logical and with another compressed bitmap
   * answer goes into container
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   *
   * The sizeInBits() of the result is equal to the maximum that of the current
   * bitmap's sizeInBits() and that of a.sizeInBits().
   */
  void logicaland(const EWAHBoolArray &a, EWAHBoolArray &container) const;

  /**
   * computes the logical and with another compressed bitmap
   * Return the answer
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   *
   * The sizeInBits() of the result is equal to the maximum that of the current
   * bitmap's sizeInBits() and that of a.sizeInBits().
   */
  EWAHBoolArray logicaland(const EWAHBoolArray &a) const {
    EWAHBoolArray answer;
    logicaland(a, answer);
    return answer;
  }

  /**
   * calls logicaland
   */
  EWAHBoolArray operator&(const EWAHBoolArray &a) const {
    return logicaland(a);
  }

  /**
   * computes the logical and with another compressed bitmap
   * answer goes into container
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   *
   * The sizeInBits() of the result should be equal to that of the current
   * bitmap irrespective of a.sizeInBits().
   *
   */
  void logicalandnot(const EWAHBoolArray &a, EWAHBoolArray &container) const;

  /**
   * calls logicalandnot
   */
  EWAHBoolArray operator-(const EWAHBoolArray &a) const {
    return logicalandnot(a);
  }

  /**
   * computes the logical and not with another compressed bitmap
   * Return the answer
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   *
   * The sizeInBits() of the result should be equal to that of the current
   * bitmap irrespective of a.sizeInBits().
   *
   */
  EWAHBoolArray logicalandnot(const EWAHBoolArray &a) const {
    EWAHBoolArray answer;
    logicalandnot(a, answer);
    return answer;
  }

  /**
   * tests whether the bitmaps "intersect" (have at least one 1-bit at the same
   * position). This function does not modify the existing bitmaps.
   * It is faster than calling logicaland.
   */
  bool intersects(const EWAHBoolArray &a) const;

  /**
   * computes the logical or with another compressed bitmap
   * answer goes into container
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   *
   * If you have many bitmaps, see fast_logicalor_tocontainer.
   *
   * The sizeInBits() of the result is equal to the maximum that of the current
   * bitmap's sizeInBits() and that of a.sizeInBits().
   */
  void logicalor(const EWAHBoolArray &a, EWAHBoolArray &container) const;

  /**
   * computes the size (in number of set bits) of the logical or with another
   * compressed bitmap
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   */
  size_t logicalorcount(const EWAHBoolArray &a) const;

  /**
   * computes the size (in number of set bits) of the logical and with another
   * compressed bitmap
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   */
  size_t logicalandcount(const EWAHBoolArray &a) const;

  /**
   * computes the size (in number of set bits) of the logical and not with
   * another compressed bitmap
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   */
  size_t logicalandnotcount(const EWAHBoolArray &a) const;

  /**
   * computes the size (in number of set bits) of the logical xor with another
   * compressed bitmap
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   */
  size_t logicalxorcount(const EWAHBoolArray &a) const;

  /**
   * computes the logical or with another compressed bitmap
   * Return the answer
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   *
   * If you have many bitmaps, see fast_logicalor.
   *
   * The sizeInBits() of the result is equal to the maximum that of the current
   * bitmap's sizeInBits() and that of a.sizeInBits().
   */
  EWAHBoolArray logicalor(const EWAHBoolArray &a) const {
    EWAHBoolArray answer;
    logicalor(a, answer);
    return answer;
  }

  /**
   * calls logicalor
   */
  EWAHBoolArray operator|(const EWAHBoolArray &a) const { return logicalor(a); }

  /**
   * computes the logical xor with another compressed bitmap
   * answer goes into container
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   *
   * The sizeInBits() of the result is equal to the maximum that of the current
   * bitmap's sizeInBits() and that of a.sizeInBits().
   */
  void logicalxor(const EWAHBoolArray &a, EWAHBoolArray &container) const;

  /**
   * computes the logical xor with another compressed bitmap
   * Return the answer
   * Running time complexity is proportional to the sum of the compressed
   * bitmap sizes.
   *
   * The sizeInBits() of the result is equal to the maximum that of the current
   * bitmap's sizeInBits() and that of a.sizeInBits().
   */
  EWAHBoolArray logicalxor(const EWAHBoolArray &a) const {
    EWAHBoolArray answer;
    logicalxor(a, answer);
    return answer;
  }

  /**
   * calls logicalxor
   */
  EWAHBoolArray operator^(const EWAHBoolArray &a) const {
    return logicalxor(a);
  }
  /**
   * clear the content of the bitmap. It does not
   * release the memory.
   */
  void reset() {
    buffer.clear();
    buffer.push_back(0);
    sizeinbits = 0;
    lastRLW = 0;
  }

  /**
   * convenience method.
   *
   * returns the number of words added (storage cost increase)
   */
  inline size_t addWord(const uword newdata,
                        const uint32_t bitsthatmatter = 8 * sizeof(uword));

  inline void printout(std::ostream &o = std::cout) {
    toBoolArray().printout(o);
  }

  /**
   * Prints a verbose description of the content of the compressed bitmap.
   */
  void debugprintout() const;

  /**
   * Return the size in bits of this bitmap (this refers
   * to the uncompressed size in bits).
   *
   * You can increase it with padWithZeroes()
   */
  inline size_t sizeInBits() const { return sizeinbits; }

  /**
   * Return the size of the buffer in bytes. This
   * is equivalent to the storage cost, minus some overhead.
   * See sizeOnDisk to get the actual storage cost with overhead.
   */
  inline size_t sizeInBytes() const { return buffer.size() * sizeof(uword); }

  /**
   * same as addEmptyWord, but you can do several in one shot!
   * returns the number of words added (storage cost increase)
   */
  size_t addStreamOfEmptyWords(const bool v, size_t number);

  /**
   * add a stream of dirty words, returns the number of words added
   * (storage cost increase)
   */
  size_t addStreamOfDirtyWords(const uword *v, const size_t number);

  /**
   * add a stream of dirty words, each one negated, returns the number of words
   * added
   * (storage cost increase)
   */
  size_t addStreamOfNegatedDirtyWords(const uword *v, const size_t number);

  /**
   * make sure the size of the array is totalbits bits by padding with zeroes.
   * returns the number of words added (storage cost increase).
   *
   * This is useful when calling "logicalnot" functions.
   *
   * This can an adverse effect of performance, especially when computing
   * intersections.
   *
   */
  size_t padWithZeroes(const size_t totalbits);

  /**
   * Compute the size on disk assuming that it was saved using
   * the method "write".
   */
  size_t sizeOnDisk(const bool savesizeinbits = true) const;

  /**
   * Save this bitmap to a stream. The file format is
   * | sizeinbits | buffer length | buffer content|
   * the sizeinbits part can be omitted if "savesizeinbits=false".
   * Both sizeinbits and buffer length are saved using the uint64_t data
   * type.
   * Returns how many bytes were handed out to the stream.
   */
  size_t write(std::ostream &out, const bool savesizeinbits = true) const;

  /**
   * same as write(std::ostream...), except that you provide a char pointer
   * and a "capacity" (in bytes). The function never writes at or beyond
   * "out+capacity". If the storage needed exceeds the given capacity, the value
   * zero is returned: it should be considered an error. Otherwise, the number
   * of bytes copied is returned.
   */
  size_t write(char *out, size_t capacity,
               const bool savesizeinbits = true) const;

  /**
   * This only writes the content of the buffer (see write()) method.
   * It is for advanced users.
   */
  void writeBuffer(std::ostream &out) const;

  /**
   * size (in words) of the underlying STL vector.
   */
  size_t bufferSize() const { return buffer.size(); }

  /**
   * this is the counterpart to the write method.
   * if you set savesizeinbits=false, then you are responsible
   * for setting the value of the attribute sizeinbits (see method
   * setSizeInBits).
   *
   * Returns how many bytes were queried from the stream.
   */
  size_t read(std::istream &in, const bool savesizeinbits = true);

  /**
   * same as read(std::istream...), except that you provide a char pointer
   * and a "capacity" (in bytes). The function never reads at or beyond
   * "in+capacity". If the detected storage exceeds the  given capacity, the
   * value zero is returned: it should be considered an error. Otherwise, the
   * number of bytes read is returned.
   */
  size_t read(const char *in, size_t capacity,
              const bool savesizeinbits = true);

  /**
   * read the buffer from a stream, see method writeBuffer.
   * this is for advanced users.
   */
  void readBuffer(std::istream &in, const size_t buffersize);

  /**
   * We define two EWAHBoolArray as being equal if they have the same set bits.
   * Alternatively, B1==B2 if and only if cardinality(B1 XOR B2) ==0.
   */
  bool operator==(const EWAHBoolArray &x) const;

  /**
   * We define two EWAHBoolArray as being different if they do not have the same
   * set bits.
   * Alternatively, B1!=B2 if and only if cardinality(B1 XOR B2) >0.
   */
  bool operator!=(const EWAHBoolArray &x) const;

  bool operator==(const BoolArray<uword> &x) const;

  bool operator!=(const BoolArray<uword> &x) const;

  /**
   * Iterate over the uncompressed words.
   * Can be considerably faster than begin()/end().
   * Running time complexity of a full scan is proportional to the
   * uncompressed size of the bitmap.
   */
  EWAHBoolArrayIterator<uword> uncompress() const;

  /**
   * To iterate over the compressed data.
   * Can be faster than any other iterator.
   * Running time complexity of a full scan is proportional to the
   * compressed size of the bitmap.
   */
  EWAHBoolArrayRawIterator<uword> raw_iterator() const;

  /**
   * Appends the content of some other compressed bitmap
   * at the end of the current bitmap.
   */
  void append(const EWAHBoolArray &x);

  /**
   * For research purposes. This computes the number of
   * dirty words and the number of compressed words.
   */
  BitmapStatistics computeStatistics() const;

  /**
   * For convenience, this fully uncompresses the bitmap.
   * Not fast!
   */
  BoolArray<uword> toBoolArray() const;

  /**
   * Convert to a list of positions of "set" bits.
   * The recommended container is vector<size_t>.
   *
   * See also toArray().
   */
  template <class container>
  void appendSetBits(container &out, const size_t offset = 0) const;

  /**
   * Returns a vector containing the position of the set
   * bits in increasing order. This just calls "toArray".
   */
  std::vector<size_t> toVector() const { return toArray(); }

  /**
   * Returns the number of bits set to the value 1.
   * The running time complexity is proportional to the
   * compressed size of the bitmap.
   *
   * This is sometimes called the cardinality.
   */
  size_t numberOfOnes() const;

  /**
   * Swap the content of this bitmap with another bitmap.
   * No copying is done. (Running time complexity is constant.)
   */
  void swap(EWAHBoolArray &x);

  const std::vector<uword> &getBuffer() const { return buffer; }

  enum { wordinbits = sizeof(uword) * 8 };

  /**
   * Please don't copy your bitmaps! The running time
   * complexity of a copy is the size of the compressed bitmap.
   **/
  EWAHBoolArray(const EWAHBoolArray &other)
      : buffer(other.buffer), sizeinbits(other.sizeinbits),
        lastRLW(other.lastRLW) {}

  /**
   * Copies the content of one bitmap onto another. Running time complexity
   * is proportional to the size of the compressed bitmap.
   * please, never hard-copy this object. Use the swap method if you must.
   */
  EWAHBoolArray &operator=(const EWAHBoolArray &x) {
    buffer = x.buffer;
    sizeinbits = x.sizeinbits;
    lastRLW = x.lastRLW;
    return *this;
  }

  /**
   * Move constructor.
   */
  EWAHBoolArray(EWAHBoolArray &&other)
      : buffer(std::move(other.buffer)), sizeinbits(other.sizeinbits),
        lastRLW(other.lastRLW) {}

  /**
   * Move assignment operator.
   */
  EWAHBoolArray &operator=(EWAHBoolArray &&x) {
    buffer = std::move(x.buffer);
    sizeinbits = x.sizeinbits;
    lastRLW = x.lastRLW;
    return *this;
  }

  /**
   * This is equivalent to the operator =. It is used
   * to keep in mind that assignment can be expensive.
   *
   *if you don't care to copy the bitmap (performance-wise), use this!
   */
  void expensive_copy(const EWAHBoolArray &x) {
    buffer = x.buffer;
    sizeinbits = x.sizeinbits;
    lastRLW = x.lastRLW;
  }

  /**
   * Write the logical not of this bitmap in the provided container.
   *
   * This function takes into account the sizeInBits value.
   * You may need to call "padWithZeroes" to adjust the sizeInBits.
   */
  void logicalnot(EWAHBoolArray &x) const;

  /**
   * Write the logical not of this bitmap in the provided container.
   *
   * This function takes into account the sizeInBits value.
   * You may need to call "padWithZeroes" to adjust the sizeInBits.
   */
  EWAHBoolArray<uword> logicalnot() const {
    EWAHBoolArray answer;
    logicalnot(answer);
    return answer;
  }

  /**
   * Apply the logical not operation on this bitmap.
   * Running time complexity is proportional to the compressed size of the
   *bitmap.
   * The current bitmap is not modified.
   *
   * This function takes into account the sizeInBits value.
   * You may need to call "padWithZeroes" to adjust the sizeInBits.
   **/
  void inplace_logicalnot();

  /**
   * set size in bits. This does not affect the compressed size. It
   * runs in constant time. This should not normally be used, except
   * as part of a deserialization process.
   */
  inline void setSizeInBits(const size_t size) { sizeinbits = size; }

  /**
   * Like addStreamOfEmptyWords but
   * addStreamOfEmptyWords but does not return the cost increase,
   * does not update sizeinbits
   */
  inline void fastaddStreamOfEmptyWords(const bool v, size_t number);
  /**
   * LikeaddStreamOfDirtyWords but does not return the cost increase,
   * does not update sizeinbits.
   */
  inline void fastaddStreamOfDirtyWords(const uword *v, const size_t number);

private:
  void assertWordCount(std::string message) const;
  void correctWordCount();
  size_t numberOfWords() const;
  // private because does not increment the size in bits
  // returns the number of words added (storage cost increase)
  inline size_t addLiteralWord(const uword newdata);

  // private because does not increment the size in bits
  // returns the number of words added (storage cost increase)
  size_t addEmptyWord(const bool v);
  // this second version "might" be faster if you hate OOP.
  // in my tests, it turned out to be slower!
  // private because does not increment the size in bits
  // inline void addEmptyWordStaticCalls(bool v);

  std::vector<uword> buffer;
  size_t sizeinbits;
  size_t lastRLW;
};
} // namespace ewah
#include "ewah-inl.h"

#endif
