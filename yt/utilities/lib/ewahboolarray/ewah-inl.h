#ifndef EWAH_INL_H
#define EWAH_INL_H

#include "ewah.h"

namespace ewah {

/**
 * computes the logical or (union) between "n" bitmaps (referenced by a
 * pointer).
 * The answer gets written out in container. This might be faster than calling
 * logicalor n-1 times.
 */
template <class uword>
void fast_logicalor_tocontainer(size_t n, const EWAHBoolArray<uword> **inputs,
                                EWAHBoolArray<uword> &container);

/**
 * computes the logical or (union) between "n" bitmaps (referenced by a
 * pointer).
 * Returns the answer. This might be faster than calling
 * logicalor n-1 times.
 */
template <class uword>
EWAHBoolArray<uword> fast_logicalor(size_t n,
                                    const EWAHBoolArray<uword> **inputs) {
  EWAHBoolArray<uword> answer;
  fast_logicalor_tocontainer(n, inputs, answer);
  return answer;
}

/**
 * Iterate over words of bits from a compressed bitmap.
 */
template <class uword> class EWAHBoolArrayIterator {
public:
  /**
   * is there a new word?
   */
  bool hasNext() const { return pointer < myparent.size(); }

  /**
   * return next word.
   */
  uword next() {
    uword returnvalue;
    if (compressedwords < rl) {
      ++compressedwords;
      if (b)
        returnvalue = notzero;
      else
        returnvalue = zero;
    } else {
      ++literalwords;
      ++pointer;
      returnvalue = myparent[pointer];
    }
    if ((compressedwords == rl) && (literalwords == lw)) {
      ++pointer;
      if (pointer < myparent.size())
        readNewRunningLengthWord();
    }
    return returnvalue;
  }

  EWAHBoolArrayIterator(const EWAHBoolArrayIterator<uword> &other)
      : pointer(other.pointer), myparent(other.myparent),
        compressedwords(other.compressedwords),
        literalwords(other.literalwords), rl(other.rl), lw(other.lw),
        b(other.b) {}

  static const uword zero = 0;
  static const uword notzero = static_cast<uword>(~zero);

private:
  EWAHBoolArrayIterator(const std::vector<uword> &parent);
  void readNewRunningLengthWord();
  friend class EWAHBoolArray<uword>;
  size_t pointer;
  const std::vector<uword> &myparent;
  uword compressedwords;
  uword literalwords;
  uword rl, lw;
  bool b;
};

/**
 * Used to go through the set bits. Not optimally fast, but convenient.
 */
template <class uword> class EWAHBoolArraySetBitForwardIterator {
public:
  typedef std::forward_iterator_tag iterator_category;
  typedef size_t *pointer;
  typedef size_t &reference_type;
  typedef size_t value_type;
  typedef ptrdiff_t difference_type;
  typedef EWAHBoolArraySetBitForwardIterator<uword> type_of_iterator;
  /**
   * Provides the location of the set bit.
   */
  inline size_t operator*() const { return answer; }

  bool operator<(const type_of_iterator &o) const {
    if (!o.hasValue)
      return true;
    if (!hasValue)
      return false;
    return answer < o.answer;
  }

  bool operator<=(const type_of_iterator &o) const {
    if (!o.hasValue)
      return true;
    if (!hasValue)
      return false;
    return answer <= o.answer;
  }

  bool operator>(const type_of_iterator &o) const { return !((*this) <= o); }

  bool operator>=(const type_of_iterator &o) const { return !((*this) < o); }

  EWAHBoolArraySetBitForwardIterator &operator++() { //++i
    if (hasNext)
      next();
    else
      hasValue = false;
    return *this;
  }

  EWAHBoolArraySetBitForwardIterator operator++(int) { // i++
    EWAHBoolArraySetBitForwardIterator old(*this);
    if (hasNext)
      next();
    else
      hasValue = false;
    return old;
  }

  bool operator==(const EWAHBoolArraySetBitForwardIterator<uword> &o) const {
    if ((!hasValue) && (!o.hasValue))
      return true;
    return (hasValue == o.hasValue) && (answer == o.answer);
  }

  bool operator!=(const EWAHBoolArraySetBitForwardIterator<uword> &o) const {
    return !(*this == o);
  }

  static EWAHBoolArraySetBitForwardIterator<uword> &end() {
    static EWAHBoolArraySetBitForwardIterator<uword> e;
    return e;
  }

  EWAHBoolArraySetBitForwardIterator(const std::vector<uword> *parent,
                                     size_t startpointer = 0)
      : word(0), position(0), runningLength(0), literalPosition(0),
        wordPosition(startpointer), wordLength(0), buffer(parent),
        hasNext(false), hasValue(false), answer(0) {
    if (wordPosition < buffer->size()) {
      setRunningLengthWord();
      hasNext = moveToNext();
      if (hasNext) {
        next();
        hasValue = true;
      }
    }
  }

  EWAHBoolArraySetBitForwardIterator()
      : word(0), position(0), runningLength(0), literalPosition(0),
        wordPosition(0), wordLength(0), buffer(NULL), hasNext(false),
        hasValue(false), answer(0) {}

  inline bool runningHasNext() const { return position < runningLength; }

  inline bool literalHasNext() {
    while (word == 0 && wordPosition < wordLength) {
      word = (*buffer)[wordPosition++];
      literalPosition = position;
      position += WORD_IN_BITS;
    }
    return word != 0;
  }

  inline void setRunningLengthWord() {
    uword rlw = (*buffer)[wordPosition];
    runningLength =
        (size_t)WORD_IN_BITS * RunningLengthWord<uword>::getRunningLength(rlw) +
        position;
    if (!RunningLengthWord<uword>::getRunningBit(rlw)) {
      position = runningLength;
    }
    wordPosition++; // point to first literal word
    wordLength =
        static_cast<uword>(wordPosition + RunningLengthWord<uword>::getNumberOfLiteralWords(rlw));
  }

  inline bool moveToNext() {
    while (!runningHasNext() && !literalHasNext()) {
      if (wordPosition >= buffer->size()) {
        return false;
      }
      setRunningLengthWord();
    }
    return true;
  }

  void next() { // update answer
    if (runningHasNext()) {
      answer = position++;
      if (runningHasNext())
        return;
    } else {
      uword t = static_cast<uword>(word & (~word + 1));
      answer = literalPosition + countOnes((UWORD)(t - 1));
      word ^= t;
    }
    hasNext = moveToNext();
  }

  enum { WORD_IN_BITS = sizeof(uword) * 8 };
  uword word; // lit word
  size_t position;
  size_t runningLength;
  size_t literalPosition;
  size_t wordPosition; // points to word in buffer
  uword wordLength;
  const std::vector<uword> *buffer;
  bool hasNext;
  bool hasValue;
  size_t answer;
};

/**
 * This object is returned by the compressed bitmap as a
 * statistical descriptor.
 */
class BitmapStatistics {
public:
  BitmapStatistics()
      : totalliteral(0), totalcompressed(0), runningwordmarker(0),
        maximumofrunningcounterreached(0) {}
  size_t getCompressedSize() const { return totalliteral + runningwordmarker; }
  size_t getUncompressedSize() const { return totalliteral + totalcompressed; }
  size_t getNumberOfDirtyWords() const { return totalliteral; }
  size_t getNumberOfCleanWords() const { return totalcompressed; }
  size_t getNumberOfMarkers() const { return runningwordmarker; }
  size_t getOverRuns() const { return maximumofrunningcounterreached; }
  size_t totalliteral;
  size_t totalcompressed;
  size_t runningwordmarker;
  size_t maximumofrunningcounterreached;
};

template <class uword> bool EWAHBoolArray<uword>::set(size_t i) {
  if (i < sizeinbits)
    return false;
  const size_t dist = (i + wordinbits) / wordinbits -
                      (sizeinbits + wordinbits - 1) / wordinbits;
  sizeinbits = i + 1;
  if (dist > 0) { // easy
    if (dist > 1) {
      fastaddStreamOfEmptyWords(false, dist - 1);
    }
    addLiteralWord(
        static_cast<uword>(static_cast<uword>(1) << (i % wordinbits)));
    return true;
  }
  RunningLengthWord<uword> lastRunningLengthWord(buffer[lastRLW]);
  if (lastRunningLengthWord.getNumberOfLiteralWords() == 0) {
    lastRunningLengthWord.setRunningLength(
        static_cast<uword>(lastRunningLengthWord.getRunningLength() - 1));
    addLiteralWord(
        static_cast<uword>(static_cast<uword>(1) << (i % wordinbits)));
    return true;
  }
  buffer[buffer.size() - 1] |=
      static_cast<uword>(static_cast<uword>(1) << (i % wordinbits));
  // check if we just completed a stream of 1s
  if (buffer[buffer.size() - 1] == static_cast<uword>(~0)) {
    // we remove the last dirty word
    buffer[buffer.size() - 1] = 0;
    buffer.resize(buffer.size() - 1);
    lastRunningLengthWord.setNumberOfLiteralWords(static_cast<uword>(
        lastRunningLengthWord.getNumberOfLiteralWords() - 1));
    // next we add one clean word
    addEmptyWord(true);
  }
  return true;
}

template <class uword> void EWAHBoolArray<uword>::inplace_logicalnot() {
  size_t pointer(0), lastrlw(0);
  while (pointer < buffer.size()) {
    RunningLengthWord<uword> rlw(buffer[pointer]);
    lastrlw = pointer; // we save this up
    if (rlw.getRunningBit())
      rlw.setRunningBit(false);
    else
      rlw.setRunningBit(true);
    ++pointer;
    for (size_t k = 0; k < rlw.getNumberOfLiteralWords(); ++k) {
      buffer[pointer] = static_cast<uword>(~buffer[pointer]);
      ++pointer;
    }
  }
  if (sizeinbits % wordinbits != 0) {
    RunningLengthWord<uword> rlw(buffer[lastrlw]);
    const uword maskbogus =
        static_cast<uword>((static_cast<uword>(1) << (sizeinbits % wordinbits)) - 1);
    if (rlw.getNumberOfLiteralWords() > 0) { // easy case
      buffer[lastrlw + 1 + rlw.getNumberOfLiteralWords() - 1] &= maskbogus;
    } else {
      rlw.setRunningLength(rlw.getRunningLength() - 1);
      addLiteralWord(maskbogus);
    }
  }
}

template <class uword> size_t EWAHBoolArray<uword>::numberOfWords() const {
  size_t tot(0);
  size_t pointer(0);
  while (pointer < buffer.size()) {
    ConstRunningLengthWord<uword> rlw(buffer[pointer]);
    tot += rlw.size();
    pointer += 1 + rlw.getNumberOfLiteralWords();
  }
  return tot;
}

template <class uword>
void EWAHBoolArray<uword>::assertWordCount(std::string message) const {
#ifdef EWAHASSERT
  size_t tot = numberOfWords();
  size_t expected = (sizeinbits + wordinbits - 1) / wordinbits;
  if (expected != tot) {
    std::cerr << "[assertWordCount] wordinbits " << wordinbits << std::endl;
    std::cerr << "[assertWordCount] sizeinbits " << sizeinbits << std::endl;
    std::cerr << "[assertWordCount] " << message << std::endl;
    std::cerr << "[assertWordCount] number of words " << tot << std::endl;
    std::cerr << "[assertWordCount] expected number of words " << expected
              << std::endl;
    debugprintout();
    throw std::runtime_error("bug");
  }
#endif
}

template <class uword> void EWAHBoolArray<uword>::correctWordCount() {
  size_t tot = numberOfWords();
  size_t expected = (sizeinbits + wordinbits - 1) / wordinbits;
  if (expected != tot) {
    if (tot < expected) {
      fastaddStreamOfEmptyWords(false, expected - tot);
    } else {
      RunningLengthWord<uword> lastRunningLengthWord(buffer[lastRLW]);
      lastRunningLengthWord.setRunningLength(static_cast<uword>(
          lastRunningLengthWord.getRunningLength() + expected - tot));
    }
  }
}

template <class uword> size_t EWAHBoolArray<uword>::numberOfOnes() const {
  size_t tot(0);
  size_t pointer(0);
  while (pointer < buffer.size()) {
    ConstRunningLengthWord<uword> rlw(buffer[pointer]);
    if (rlw.getRunningBit()) {
      tot += static_cast<size_t>(rlw.getRunningLength() * wordinbits);
    }
    ++pointer;
    for (size_t k = 0; k < rlw.getNumberOfLiteralWords(); ++k) {
      tot += countOnes((UWORD)buffer[pointer]);
      ++pointer;
    }
  }
  return tot;
}

template <class uword>
std::vector<size_t> EWAHBoolArray<uword>::toArray() const {
  std::vector<size_t> ans;
  size_t pos(0);
  size_t pointer(0);
  const size_t buffersize = buffer.size();
  while (pointer < buffersize) {
    ConstRunningLengthWord<uword> rlw(buffer[pointer]);
    const size_t productofrl =
        static_cast<size_t>(rlw.getRunningLength() * wordinbits);
    if (rlw.getRunningBit()) {
      size_t upper_limit = pos + productofrl;
      for (; pos < upper_limit; ++pos) {
        ans.push_back(pos);
      }
    } else {
      pos += productofrl;
    }
    ++pointer;
    const size_t rlwlw = rlw.getNumberOfLiteralWords();
    for (size_t k = 0; k < rlwlw; ++k) {
      uword myword = buffer[pointer];
      while (myword != 0) {
        uint64_t t = myword & (~myword + 1);
        uint32_t r = numberOfTrailingZeros(t);
        ans.push_back(pos + r);
        myword ^= t;
      }
      pos += wordinbits;
      ++pointer;
    }
  }
  return ans;
}

template <class uword>
void EWAHBoolArray<uword>::logicalnot(EWAHBoolArray &x) const {
  x.reset();
  x.buffer.reserve(buffer.size());
  EWAHBoolArrayRawIterator<uword> i = this->raw_iterator();
  if (!i.hasNext())
    return; // nothing to do
  while (true) {
    BufferedRunningLengthWord<uword> &rlw = i.next();
    if (i.hasNext()) {
      if (rlw.getRunningLength() > 0)
        x.fastaddStreamOfEmptyWords(!rlw.getRunningBit(),
                                    rlw.getRunningLength());
      if (rlw.getNumberOfLiteralWords() > 0) {
        const uword *dw = i.dirtyWords();
        for (size_t k = 0; k < rlw.getNumberOfLiteralWords(); ++k) {
          x.addLiteralWord(~dw[k]);
        }
      }
    } else {
      if (rlw.getNumberOfLiteralWords() == 0) {
        if ((this->sizeinbits % wordinbits != 0) && !rlw.getRunningBit()) {
          if (rlw.getRunningLength() > 1)
            x.fastaddStreamOfEmptyWords(!rlw.getRunningBit(),
                                        rlw.getRunningLength() - 1);
          const uword maskbogus =
              static_cast<uword>((static_cast<uword>(1) << (this->sizeinbits % wordinbits)) - 1);
          x.addLiteralWord(maskbogus);
          break;
        } else {
          if (rlw.getRunningLength() > 0)
            x.fastaddStreamOfEmptyWords(!rlw.getRunningBit(),
                                        rlw.getRunningLength());
          break;
        }
      }
      if (rlw.getRunningLength() > 0)
        x.fastaddStreamOfEmptyWords(!rlw.getRunningBit(),
                                    rlw.getRunningLength());
      const uword *dw = i.dirtyWords();
      for (size_t k = 0; k + 1 < rlw.getNumberOfLiteralWords(); ++k) {
        x.addLiteralWord(~dw[k]);
      }
      const uword maskbogus =
          (this->sizeinbits % wordinbits != 0)
              ? static_cast<uword>((static_cast<uword>(1) << (this->sizeinbits % wordinbits)) - 1)
              : ~static_cast<uword>(0);
      x.addLiteralWord(static_cast<uword>((~dw[rlw.getNumberOfLiteralWords() - 1]) & maskbogus));
      break;
    }
  }
  x.sizeinbits = this->sizeinbits;
}

template <class uword>
size_t EWAHBoolArray<uword>::addWord(const uword newdata,
                                     const uint32_t bitsthatmatter) {
  sizeinbits += bitsthatmatter;
  if (newdata == 0) {
    return addEmptyWord(0);
  } else if (newdata == static_cast<uword>(~0)) {
    return addEmptyWord(1);
  } else {
    return addLiteralWord(newdata);
  }
}

template <class uword>
inline void EWAHBoolArray<uword>::writeBuffer(std::ostream &out) const {
  if (!buffer.empty())
    out.write(reinterpret_cast<const char *>(&buffer[0]),
              sizeof(uword) * buffer.size());
}

template <class uword>
inline void EWAHBoolArray<uword>::readBuffer(std::istream &in,
                                             const size_t buffersize) {
  buffer.resize(buffersize);
  if (buffersize > 0)
    in.read(reinterpret_cast<char *>(&buffer[0]), sizeof(uword) * buffersize);
}

template <class uword>
size_t EWAHBoolArray<uword>::write(std::ostream &out,
                                   const bool savesizeinbits) const {
  size_t written = 0;
  if (savesizeinbits) {
    uint64_t sb = static_cast<uint64_t>(sizeinbits);
    out.write(reinterpret_cast<const char *>(&sb), sizeof(sb));
    written += sizeof(uint64_t);
  }
  const size_t buffersize = buffer.size();
  uint64_t bs = static_cast<uint64_t>(buffersize);
  out.write(reinterpret_cast<const char *>(&bs), sizeof(bs));
  written += sizeof(uint64_t);

  if (buffersize > 0) {
    out.write(reinterpret_cast<const char *>(&buffer[0]),
              static_cast<std::streamsize>(sizeof(uword) * buffersize));
    written += sizeof(uword) * buffersize;
  }
  return written;
}

template <class uword>
size_t EWAHBoolArray<uword>::write(char *out, size_t capacity,
                                   const bool savesizeinbits) const {
  size_t written = 0;
  if (savesizeinbits) {
    uint64_t sb = static_cast<uint64_t>(sizeinbits);
    if (capacity < sizeof(sb))
      return 0;
    capacity -= sizeof(sb);
    memcpy(out, &sb, sizeof(sb));
    out += sizeof(sb);
    written += sizeof(uint64_t);
  }
  const size_t buffersize = buffer.size();
  uint64_t bs = static_cast<uint64_t>(buffersize);
  if (capacity < sizeof(bs))
    return 0;
  capacity -= sizeof(bs);
  memcpy(out, &buffersize, sizeof(bs));
  out += sizeof(bs);
  written += sizeof(uint64_t);

  if (buffersize > 0) {
    if (capacity < sizeof(uword) * buffersize)
      return 0;
    memcpy(out, &buffer[0], sizeof(uword) * buffersize);
    written += sizeof(uword) * buffersize;
  }
  return written;
}

template <class uword>
size_t EWAHBoolArray<uword>::read(std::istream &in, const bool savesizeinbits) {
  size_t read = 0;
  if (savesizeinbits) {
    uint64_t tmp;
    in.read(reinterpret_cast<char *>(&tmp), sizeof(tmp));
    read += sizeof(tmp);
    sizeinbits = static_cast<size_t>(tmp);
  } else {
    sizeinbits = 0;
  }
  size_t buffersize(0);
  uint64_t tmp;
  in.read(reinterpret_cast<char *>(&tmp), sizeof(tmp));
  read += sizeof(tmp);
  buffersize = static_cast<size_t>(tmp);
  buffer.resize(buffersize);
  if (buffersize > 0) {
    in.read(reinterpret_cast<char *>(&buffer[0]),
            static_cast<std::streamsize>(sizeof(uword) * buffersize));
    read += sizeof(uword) * buffersize;
  }
  return read;
}

template <class uword>
size_t EWAHBoolArray<uword>::read(const char *in, size_t capacity,
                                  const bool savesizeinbits) {
  size_t read = 0;
  if (savesizeinbits) {
    uint64_t tmp;
    if (capacity < sizeof(tmp))
      return 0;
    capacity -= sizeof(tmp);
    memcpy(reinterpret_cast<char *>(&tmp), in, sizeof(tmp));
    read += sizeof(tmp);
    in += sizeof(tmp);
    sizeinbits = static_cast<size_t>(tmp);
  } else {
    sizeinbits = 0;
  }
  size_t buffersize(0);
  uint64_t tmp;
  if (capacity < sizeof(uint64_t))
    return 0;
  capacity -= sizeof(uint64_t);
  memcpy(reinterpret_cast<char *>(&tmp), in, sizeof(uint64_t));
  in += sizeof(uint64_t);
  read += sizeof(uint64_t);
  buffersize = static_cast<size_t>(tmp);
  buffer.resize(buffersize);
  if (buffersize > 0) {
    if (capacity < sizeof(uword) * buffersize)
      return 0;
    memcpy(&buffer[0], in, sizeof(uword) * buffersize);
    read += sizeof(uword) * buffersize;
  }
  return read;
}

template <class uword>
size_t EWAHBoolArray<uword>::addLiteralWord(const uword newdata) {
  RunningLengthWord<uword> lastRunningLengthWord(buffer[lastRLW]);
  uword numbersofar = lastRunningLengthWord.getNumberOfLiteralWords();
  if (numbersofar >=
      RunningLengthWord<uword>::largestliteralcount) { // 0x7FFF) {
    buffer.push_back(0);
    lastRLW = buffer.size() - 1;
    RunningLengthWord<uword> lastRunningLengthWord2(buffer[lastRLW]);
    lastRunningLengthWord2.setNumberOfLiteralWords(1);
    buffer.push_back(newdata);
    return 2;
  }
  lastRunningLengthWord.setNumberOfLiteralWords(
      static_cast<uword>(numbersofar + 1));
  buffer.push_back(newdata);
  return 1;
}

template <class uword>
size_t EWAHBoolArray<uword>::padWithZeroes(const size_t totalbits) {
  size_t wordsadded = 0;
  if (totalbits <= sizeinbits)
    return wordsadded;

  size_t missingbits = totalbits - sizeinbits;

  RunningLengthWord<uword> rlw(buffer[lastRLW]);
  if (rlw.getNumberOfLiteralWords() > 0) {
    // Consume trailing zeroes of trailing literal word (past sizeinbits)
    size_t remain = sizeinbits % wordinbits;
    if (remain > 0) // Is last word partial?
    {
      size_t avail = wordinbits - remain;
      if (avail > 0) {
        if (missingbits > avail) {
          missingbits -= avail;
        } else {
          missingbits = 0;
        }
        sizeinbits += avail;
      }
    }
  }

  if (missingbits > 0) {
    size_t wordstoadd = missingbits / wordinbits;
    if ((missingbits % wordinbits) != 0)
      ++wordstoadd;

    wordsadded = addStreamOfEmptyWords(false, wordstoadd);
  }
  sizeinbits = totalbits;
  return wordsadded;
}

/**
 * This is a low-level iterator.
 */

template <class uword = uint32_t> class EWAHBoolArrayRawIterator {
public:
  EWAHBoolArrayRawIterator(const EWAHBoolArray<uword> &p)
      : pointer(0), myparent(&p.getBuffer()), rlw((*myparent)[pointer], this) {}
  EWAHBoolArrayRawIterator(const EWAHBoolArrayRawIterator &o)
      : pointer(o.pointer), myparent(o.myparent), rlw(o.rlw) {}

  bool hasNext() const { return pointer < myparent->size(); }

  BufferedRunningLengthWord<uword> &next() {
    rlw.read((*myparent)[pointer]);
    pointer = static_cast<size_t>(pointer + rlw.getNumberOfLiteralWords() + 1);
    return rlw;
  }

  const uword *dirtyWords() const {
    return myparent->data() +
           static_cast<size_t>(pointer - rlw.getNumberOfLiteralWords());
  }

  EWAHBoolArrayRawIterator &operator=(const EWAHBoolArrayRawIterator &other) {
    pointer = other.pointer;
    myparent = other.myparent;
    rlw = other.rlw;
    return *this;
  }

  size_t pointer;
  const std::vector<uword> *myparent;
  BufferedRunningLengthWord<uword> rlw;

  EWAHBoolArrayRawIterator();
};

template <class uword>
EWAHBoolArrayIterator<uword> EWAHBoolArray<uword>::uncompress() const {
  return EWAHBoolArrayIterator<uword>(buffer);
}

template <class uword>
EWAHBoolArrayRawIterator<uword> EWAHBoolArray<uword>::raw_iterator() const {
  return EWAHBoolArrayRawIterator<uword>(*this);
}

template <class uword>
bool EWAHBoolArray<uword>::operator==(const EWAHBoolArray &x) const {
  EWAHBoolArrayRawIterator<uword> i = x.raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = raw_iterator();
  if (!(i.hasNext() and j.hasNext())) { // hopefully this never happens...
    return (i.hasNext() == false) && (j.hasNext() == false);
  }
  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();

  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey = i_is_prey ? rlwi : rlwj;
      BufferedRunningLengthWord<uword> &predator = i_is_prey ? rlwj : rlwi;
      size_t index = 0;
      const bool nonzero =
          ((!predator.getRunningBit())
               ? prey.nonzero_discharge(predator.getRunningLength(), index)
               : prey.nonzero_dischargeNegated(predator.getRunningLength(),
                                               index));
      if (nonzero) {
        return false;
      }
      if (predator.getRunningLength() - index > 0) {
        if (predator.getRunningBit()) {
          return false;
        }
      }
      predator.discardRunningWordsWithReload();
    }
    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k)
        if ((rlwi.getLiteralWordAt(k) ^ rlwj.getLiteralWordAt(k)) != 0)
          return false;
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  const bool i_remains = rlwi.size() > 0;
  BufferedRunningLengthWord<uword> &remaining = i_remains ? rlwi : rlwj;
  return !remaining.nonzero_discharge();
}

template <class uword> void EWAHBoolArray<uword>::swap(EWAHBoolArray &x) {
  buffer.swap(x.buffer);
  size_t tmp = x.sizeinbits;
  x.sizeinbits = sizeinbits;
  sizeinbits = tmp;
  tmp = x.lastRLW;
  x.lastRLW = lastRLW;
  lastRLW = tmp;
}

template <class uword>
void EWAHBoolArray<uword>::append(const EWAHBoolArray &x) {
  if (sizeinbits % wordinbits == 0) {
    // hoping for the best?
    sizeinbits += x.sizeinbits;
    ConstRunningLengthWord<uword> lRLW(buffer[lastRLW]);
    if ((lRLW.getRunningLength() == 0) &&
        (lRLW.getNumberOfLiteralWords() == 0)) {
      // it could be that the running length word is empty, in such a case,
      // we want to get rid of it!
      lastRLW = x.lastRLW + buffer.size() - 1;
      buffer.resize(buffer.size() - 1);
      buffer.insert(buffer.end(), x.buffer.begin(), x.buffer.end());
    } else {
      lastRLW = x.lastRLW + buffer.size();
      buffer.insert(buffer.end(), x.buffer.begin(), x.buffer.end());
    }
  } else {
    std::stringstream ss;
    ss << "This should really not happen! You are trying to append to a bitmap "
          "having a fractional number of words, that is,  "
       << static_cast<int>(sizeinbits) << " bits with a word size in bits of "
       << static_cast<int>(wordinbits) << ". ";
    ss << "Size of the bitmap being appended: " << x.sizeinbits << " bits."
       << std::endl;
    throw std::invalid_argument(ss.str());
  }
}

template <class uword>
EWAHBoolArrayIterator<uword>::EWAHBoolArrayIterator(
    const std::vector<uword> &parent)
    : pointer(0), myparent(parent), compressedwords(0), literalwords(0), rl(0),
      lw(0), b(0) {
  if (pointer < myparent.size())
    readNewRunningLengthWord();
}

template <class uword>
void EWAHBoolArrayIterator<uword>::readNewRunningLengthWord() {
  literalwords = 0;
  compressedwords = 0;
  ConstRunningLengthWord<uword> rlw(myparent[pointer]);
  rl = rlw.getRunningLength();
  lw = rlw.getNumberOfLiteralWords();
  b = rlw.getRunningBit();
  if ((rl == 0) && (lw == 0)) {
    if (pointer < myparent.size() - 1) {
      ++pointer;
      readNewRunningLengthWord();
    } else {
      pointer = myparent.size();
    }
  }
}

template <class uword>
BoolArray<uword> EWAHBoolArray<uword>::toBoolArray() const {
  BoolArray<uword> ans(sizeinbits);
  EWAHBoolArrayIterator<uword> i = uncompress();
  size_t counter = 0;
  while (i.hasNext()) {
    ans.setWord(counter++, i.next());
  }
  return ans;
}

template <class uword>
template <class container>
void EWAHBoolArray<uword>::appendSetBits(container &out,
                                         const size_t offset) const {
  size_t pointer(0);
  size_t currentoffset(offset);
  if (RESERVEMEMORY)
    out.reserve(buffer.size() + 64); // trading memory for speed.
  const size_t buffersize = buffer.size();
  while (pointer < buffersize) {
    ConstRunningLengthWord<uword> rlw(buffer[pointer]);
    const size_t productofrl =
        static_cast<size_t>(rlw.getRunningLength() * wordinbits);
    if (rlw.getRunningBit()) {
      const size_t upper_limit = currentoffset + productofrl;
      for (; currentoffset < upper_limit; ++currentoffset) {
        out.push_back(currentoffset);
      }
    } else {
      currentoffset += productofrl;
    }
    ++pointer;
    const size_t rlwlw = rlw.getNumberOfLiteralWords();
    for (uword k = 0; k < rlwlw; ++k) {
      uword currentword = buffer[pointer];
      while (currentword != 0) {
        uword t = static_cast<uword>(currentword & (~currentword+1));
        uint32_t r = numberOfTrailingZeros(t);
        out.push_back(currentoffset + r);
        currentword ^= t;
      }
      currentoffset += wordinbits;
      ++pointer;
    }
  }
}

template <class uword>
bool EWAHBoolArray<uword>::operator!=(const EWAHBoolArray<uword> &x) const {
  return !(*this == x);
}

template <class uword>
bool EWAHBoolArray<uword>::operator==(const BoolArray<uword> &x) const {
  // could be more efficient
  return (this->toBoolArray() == x);
}

template <class uword>
bool EWAHBoolArray<uword>::operator!=(const BoolArray<uword> &x) const {
  // could be more efficient
  return (this->toBoolArray() != x);
}

template <class uword>
size_t EWAHBoolArray<uword>::addStreamOfEmptyWords(const bool v,
                                                   size_t number) {
  if (number == 0)
    return 0;
  sizeinbits += number * wordinbits;
  size_t wordsadded = 0;
  if ((RunningLengthWord<uword>::getRunningBit(buffer[lastRLW]) != v) &&
      (RunningLengthWord<uword>::size(buffer[lastRLW]) == 0)) {
    RunningLengthWord<uword>::setRunningBit(buffer[lastRLW], v);
  } else if ((RunningLengthWord<uword>::getNumberOfLiteralWords(
                  buffer[lastRLW]) != 0) ||
             (RunningLengthWord<uword>::getRunningBit(buffer[lastRLW]) != v)) {
    buffer.push_back(0);
    ++wordsadded;
    lastRLW = buffer.size() - 1;
    if (v)
      RunningLengthWord<uword>::setRunningBit(buffer[lastRLW], v);
  }
  const uword runlen =
      RunningLengthWord<uword>::getRunningLength(buffer[lastRLW]);

  const uword whatwecanadd =
      number < static_cast<size_t>(
                   RunningLengthWord<uword>::largestrunninglengthcount - runlen)
          ? static_cast<uword>(number)
          : static_cast<uword>(
                RunningLengthWord<uword>::largestrunninglengthcount - runlen);
  RunningLengthWord<uword>::setRunningLength(
      buffer[lastRLW], static_cast<uword>(runlen + whatwecanadd));

  number -= static_cast<size_t>(whatwecanadd);
  while (number >= RunningLengthWord<uword>::largestrunninglengthcount) {
    buffer.push_back(0);
    ++wordsadded;
    lastRLW = buffer.size() - 1;
    if (v)
      RunningLengthWord<uword>::setRunningBit(buffer[lastRLW], v);
    RunningLengthWord<uword>::setRunningLength(
        buffer[lastRLW], RunningLengthWord<uword>::largestrunninglengthcount);
    number -= static_cast<size_t>(
        RunningLengthWord<uword>::largestrunninglengthcount);
  }
  if (number > 0) {
    buffer.push_back(0);
    ++wordsadded;
    lastRLW = buffer.size() - 1;
    if (v)
      RunningLengthWord<uword>::setRunningBit(buffer[lastRLW], v);
    RunningLengthWord<uword>::setRunningLength(buffer[lastRLW],
                                               static_cast<uword>(number));
  }
  return wordsadded;
}

template <class uword>
void EWAHBoolArray<uword>::fastaddStreamOfEmptyWords(const bool v,
                                                     size_t number) {
  if (number == 0)
    return;
  if ((RunningLengthWord<uword>::getRunningBit(buffer[lastRLW]) != v) &&
      (RunningLengthWord<uword>::size(buffer[lastRLW]) == 0)) {
    RunningLengthWord<uword>::setRunningBit(buffer[lastRLW], v);
  } else if ((RunningLengthWord<uword>::getNumberOfLiteralWords(
                  buffer[lastRLW]) != 0) ||
             (RunningLengthWord<uword>::getRunningBit(buffer[lastRLW]) != v)) {
    buffer.push_back(0);
    lastRLW = buffer.size() - 1;
    if (v)
      RunningLengthWord<uword>::setRunningBit(buffer[lastRLW], v);
  }
  const uword runlen =
      RunningLengthWord<uword>::getRunningLength(buffer[lastRLW]);

  const uword whatwecanadd =
      number < static_cast<size_t>(
                   RunningLengthWord<uword>::largestrunninglengthcount - runlen)
          ? static_cast<uword>(number)
          : static_cast<uword>(
                RunningLengthWord<uword>::largestrunninglengthcount - runlen);
  RunningLengthWord<uword>::setRunningLength(
      buffer[lastRLW], static_cast<uword>(runlen + whatwecanadd));

  number -= static_cast<size_t>(whatwecanadd);
  while (number >= RunningLengthWord<uword>::largestrunninglengthcount) {
    buffer.push_back(0);
    lastRLW = buffer.size() - 1;
    if (v)
      RunningLengthWord<uword>::setRunningBit(buffer[lastRLW], v);
    RunningLengthWord<uword>::setRunningLength(
        buffer[lastRLW], RunningLengthWord<uword>::largestrunninglengthcount);
    number -= static_cast<size_t>(
        RunningLengthWord<uword>::largestrunninglengthcount);
  }
  if (number > 0) {
    buffer.push_back(0);
    lastRLW = buffer.size() - 1;
    if (v)
      RunningLengthWord<uword>::setRunningBit(buffer[lastRLW], v);
    RunningLengthWord<uword>::setRunningLength(buffer[lastRLW],
                                               static_cast<uword>(number));
  }
}

template <class uword>
size_t EWAHBoolArray<uword>::addStreamOfDirtyWords(const uword *v,
                                                   const size_t number) {
  if (number == 0)
    return 0;
  uword rlw = buffer[lastRLW];
  size_t NumberOfLiteralWords =
      RunningLengthWord<uword>::getNumberOfLiteralWords(rlw);
  if (NumberOfLiteralWords + number <=
      RunningLengthWord<uword>::largestliteralcount) {
    RunningLengthWord<uword>::setNumberOfLiteralWords(
        rlw, static_cast<uword>(NumberOfLiteralWords + number));
    buffer[lastRLW] = rlw;
    sizeinbits += number * wordinbits;
    buffer.insert(buffer.end(), v, v + number);
    return number;
  }
  // we proceed the long way
  size_t howmanywecanadd =
      RunningLengthWord<uword>::largestliteralcount - NumberOfLiteralWords;
  RunningLengthWord<uword>::setNumberOfLiteralWords(
      rlw, RunningLengthWord<uword>::largestliteralcount);
  buffer[lastRLW] = rlw;
  buffer.insert(buffer.end(), v, v + howmanywecanadd);
  size_t wordadded = howmanywecanadd;
  sizeinbits += howmanywecanadd * wordinbits;
  buffer.push_back(0);
  lastRLW = buffer.size() - 1;
  ++wordadded;
  wordadded +=
      addStreamOfDirtyWords(v + howmanywecanadd, number - howmanywecanadd);
  return wordadded;
}

template <class uword>
void EWAHBoolArray<uword>::fastaddStreamOfDirtyWords(const uword *v,
                                                     const size_t number) {
  if (number == 0)
    return;
  uword rlw = buffer[lastRLW];
  size_t NumberOfLiteralWords =
      RunningLengthWord<uword>::getNumberOfLiteralWords(rlw);
  if (NumberOfLiteralWords + number <=
      RunningLengthWord<uword>::largestliteralcount) {
    RunningLengthWord<uword>::setNumberOfLiteralWords(
        rlw, static_cast<uword>(NumberOfLiteralWords + number));
    buffer[lastRLW] = rlw;
    for (size_t i = 0; i < number; ++i)
      buffer.push_back(v[i]);
    // buffer.insert(buffer.end(), v, v+number); // seems slower than push_back?
    return;
  }
  // we proceed the long way
  size_t howmanywecanadd =
      RunningLengthWord<uword>::largestliteralcount - NumberOfLiteralWords;
  RunningLengthWord<uword>::setNumberOfLiteralWords(
      rlw, RunningLengthWord<uword>::largestliteralcount);
  buffer[lastRLW] = rlw;
  for (size_t i = 0; i < howmanywecanadd; ++i)
    buffer.push_back(v[i]);
  // buffer.insert(buffer.end(), v, v+howmanywecanadd);// seems slower than
  // push_back?
  buffer.push_back(0);
  lastRLW = buffer.size() - 1;
  fastaddStreamOfDirtyWords(v + howmanywecanadd, number - howmanywecanadd);
}

template <class uword>
size_t EWAHBoolArray<uword>::addStreamOfNegatedDirtyWords(const uword *v,
                                                          const size_t number) {
  if (number == 0)
    return 0;
  uword rlw = buffer[lastRLW];
  size_t NumberOfLiteralWords =
      RunningLengthWord<uword>::getNumberOfLiteralWords(rlw);
  if (NumberOfLiteralWords + number <=
      RunningLengthWord<uword>::largestliteralcount) {
    RunningLengthWord<uword>::setNumberOfLiteralWords(
        rlw, static_cast<uword>(NumberOfLiteralWords + number));
    buffer[lastRLW] = rlw;
    sizeinbits += number * wordinbits;
    for (size_t k = 0; k < number; ++k)
      buffer.push_back(~v[k]);
    return number;
  }
  // we proceed the long way
  size_t howmanywecanadd =
      RunningLengthWord<uword>::largestliteralcount - NumberOfLiteralWords;
  RunningLengthWord<uword>::setNumberOfLiteralWords(
      rlw, RunningLengthWord<uword>::largestliteralcount);
  buffer[lastRLW] = rlw;
  for (size_t k = 0; k < howmanywecanadd; ++k)
    buffer.push_back(~v[k]);
  size_t wordadded = howmanywecanadd;
  sizeinbits += howmanywecanadd * wordinbits;
  buffer.push_back(0);
  lastRLW = buffer.size() - 1;
  ++wordadded;
  wordadded +=
      addStreamOfDirtyWords(v + howmanywecanadd, number - howmanywecanadd);
  return wordadded;
}

template <class uword> size_t EWAHBoolArray<uword>::addEmptyWord(const bool v) {
  RunningLengthWord<uword> lastRunningLengthWord(buffer[lastRLW]);
  const bool noliteralword =
      (lastRunningLengthWord.getNumberOfLiteralWords() == 0);
  // first, if the last running length word is empty, we align it
  // this
  uword runlen = lastRunningLengthWord.getRunningLength();
  if ((noliteralword) && (runlen == 0)) {
    lastRunningLengthWord.setRunningBit(v);
  }
  if ((noliteralword) && (lastRunningLengthWord.getRunningBit() == v) &&
      (runlen < RunningLengthWord<uword>::largestrunninglengthcount)) {
    lastRunningLengthWord.setRunningLength(static_cast<uword>(runlen + 1));
    return 0;
  } else {
    // we have to start anew
    buffer.push_back(0);
    lastRLW = buffer.size() - 1;
    RunningLengthWord<uword> lastRunningLengthWord2(buffer[lastRLW]);
    lastRunningLengthWord2.setRunningBit(v);
    lastRunningLengthWord2.setRunningLength(1);
    return 1;
  }
}

template <class uword>
void fast_logicalor_tocontainer(size_t n, const EWAHBoolArray<uword> **inputs,
                                EWAHBoolArray<uword> &container) {
  class EWAHBoolArrayPtr {

  public:
    EWAHBoolArrayPtr(const EWAHBoolArray<uword> *p, bool o) : ptr(p), own(o) {}
    const EWAHBoolArray<uword> *ptr;
    bool own; // whether to clean

    bool operator<(const EWAHBoolArrayPtr &o) const {
      return o.ptr->sizeInBytes() < ptr->sizeInBytes(); // backward on purpose
    }
  };

  if (n == 0) {
    container.reset();
    return;
  }
  if (n == 1) {
    container = *inputs[0];
    return;
  }
  std::priority_queue<EWAHBoolArrayPtr> pq;
  for (size_t i = 0; i < n; i++) {
    // could use emplace
    pq.push(EWAHBoolArrayPtr(inputs[i], false));
  }
  while (pq.size() > 2) {

    EWAHBoolArrayPtr x1 = pq.top();
    pq.pop();

    EWAHBoolArrayPtr x2 = pq.top();
    pq.pop();

    EWAHBoolArray<uword> *buffer = new EWAHBoolArray<uword>();
    x1.ptr->logicalor(*x2.ptr, *buffer);

    if (x1.own) {
      delete x1.ptr;
    }
    if (x2.own) {
      delete x2.ptr;
    }
    pq.push(EWAHBoolArrayPtr(buffer, true));
  }
  EWAHBoolArrayPtr x1 = pq.top();
  pq.pop();

  EWAHBoolArrayPtr x2 = pq.top();
  pq.pop();

  x1.ptr->logicalor(*x2.ptr, container);

  if (x1.own) {
    delete x1.ptr;
  }
  if (x2.own) {
    delete x2.ptr;
  }
}

template <class uword>
void EWAHBoolArray<uword>::logicalor(const EWAHBoolArray &a,
                                     EWAHBoolArray &container) const {
  container.reset();
  if (RESERVEMEMORY)
    container.buffer.reserve(buffer.size() + a.buffer.size());
  EWAHBoolArrayRawIterator<uword> i = a.raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = raw_iterator();
  if (!(i.hasNext() and j.hasNext())) { // hopefully this never happens...
    container.setSizeInBits(sizeInBits());
    return;
  }
  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();

  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey = i_is_prey ? rlwi : rlwj;
      BufferedRunningLengthWord<uword> &predator = i_is_prey ? rlwj : rlwi;
      if (predator.getRunningBit()) {
        container.fastaddStreamOfEmptyWords(true, predator.getRunningLength());
        prey.discardFirstWordsWithReload(predator.getRunningLength());
      } else {
        const size_t index =
            prey.discharge(container, predator.getRunningLength());
        container.fastaddStreamOfEmptyWords(false, predator.getRunningLength() -
                                                       index);
      }
      predator.discardRunningWordsWithReload();
    }

    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k) {
        container.addWord(rlwi.getLiteralWordAt(k) | rlwj.getLiteralWordAt(k));
      }
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  const bool i_remains = rlwi.size() > 0;
  BufferedRunningLengthWord<uword> &remaining = i_remains ? rlwi : rlwj;
  remaining.discharge(container);
  container.setSizeInBits(sizeInBits() > a.sizeInBits() ? sizeInBits()
                                                        : a.sizeInBits());
}

template <class uword>
size_t EWAHBoolArray<uword>::logicalorcount(const EWAHBoolArray &a) const {
  size_t answer = 0;
  EWAHBoolArrayRawIterator<uword> i = a.raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = raw_iterator();
  if (!(i.hasNext() and j.hasNext())) { // hopefully this never happens...
    return 0;
  }
  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();

  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey = i_is_prey ? rlwi : rlwj;
      BufferedRunningLengthWord<uword> &predator = i_is_prey ? rlwj : rlwi;
      if (predator.getRunningBit()) {
        answer += predator.getRunningLength() * wordinbits;
        prey.discardFirstWordsWithReload(predator.getRunningLength());

      } else {
        // const size_t index =
        prey.dischargeCount(predator.getRunningLength(), &answer);
      }
      predator.discardRunningWordsWithReload();
    }

    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k) {
        answer += countOnes(
            (uword)(rlwi.getLiteralWordAt(k) | rlwj.getLiteralWordAt(k)));
      }
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  const bool i_remains = rlwi.size() > 0;
  BufferedRunningLengthWord<uword> &remaining = i_remains ? rlwi : rlwj;
  answer += remaining.dischargeCount();
  return answer;
}

template <class uword>
void EWAHBoolArray<uword>::logicalxor(const EWAHBoolArray &a,
                                      EWAHBoolArray &container) const {
  container.reset();
  if (RESERVEMEMORY)
    container.buffer.reserve(buffer.size() + a.buffer.size());
  EWAHBoolArrayRawIterator<uword> i = a.raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = raw_iterator();
  if (!(i.hasNext() and j.hasNext())) { // hopefully this never happens...
    container.setSizeInBits(sizeInBits());
    return;
  }
  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();
  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey = i_is_prey ? rlwi : rlwj;
      BufferedRunningLengthWord<uword> &predator = i_is_prey ? rlwj : rlwi;
      const size_t index =
          (!predator.getRunningBit())
              ? prey.discharge(container, predator.getRunningLength())
              : prey.dischargeNegated(container, predator.getRunningLength());
      container.fastaddStreamOfEmptyWords(predator.getRunningBit(),
                                          predator.getRunningLength() - index);
      predator.discardRunningWordsWithReload();
    }
    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k)
        container.addWord(rlwi.getLiteralWordAt(k) ^ rlwj.getLiteralWordAt(k));
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  const bool i_remains = rlwi.size() > 0;
  BufferedRunningLengthWord<uword> &remaining = i_remains ? rlwi : rlwj;
  remaining.discharge(container);
  container.setSizeInBits(sizeInBits() > a.sizeInBits() ? sizeInBits()
                                                        : a.sizeInBits());
}

template <class uword>
size_t EWAHBoolArray<uword>::logicalxorcount(const EWAHBoolArray &a) const {
  EWAHBoolArrayRawIterator<uword> i = a.raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = raw_iterator();
  if (!i.hasNext())
    return a.numberOfOnes();
  if (!j.hasNext())
    return this->numberOfOnes();

  size_t answer = 0;

  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();
  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey = i_is_prey ? rlwi : rlwj;
      BufferedRunningLengthWord<uword> &predator = i_is_prey ? rlwj : rlwi;
      size_t index;

      if (predator.getRunningBit()) {
        index =
            prey.dischargeCountNegated(predator.getRunningLength(), &answer);
      } else {
        index = prey.dischargeCount(predator.getRunningLength(), &answer);
      }
      if (predator.getRunningBit())
        answer += (predator.getRunningLength() - index) * wordinbits;

      predator.discardRunningWordsWithReload();
    }
    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k) {
        answer += countOnes(
            (uword)(rlwi.getLiteralWordAt(k) ^ rlwj.getLiteralWordAt(k)));
      }
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  const bool i_remains = rlwi.size() > 0;
  BufferedRunningLengthWord<uword> &remaining = i_remains ? rlwi : rlwj;
  answer += remaining.dischargeCount();
  return answer;
}

template <class uword>
void EWAHBoolArray<uword>::logicaland(const EWAHBoolArray &a,
                                      EWAHBoolArray &container) const {
  container.reset();
  if (RESERVEMEMORY)
    container.buffer.reserve(buffer.size() > a.buffer.size() ? buffer.size()
                                                             : a.buffer.size());
  EWAHBoolArrayRawIterator<uword> i = a.raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = raw_iterator();
  if (!(i.hasNext() and j.hasNext())) { // hopefully this never happens...
    container.setSizeInBits(sizeInBits());
    return;
  }
  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();

  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey(i_is_prey ? rlwi : rlwj);
      BufferedRunningLengthWord<uword> &predator(i_is_prey ? rlwj : rlwi);
      if (!predator.getRunningBit()) {
        container.fastaddStreamOfEmptyWords(false, predator.getRunningLength());
        prey.discardFirstWordsWithReload(predator.getRunningLength());
      } else {
        const size_t index =
            prey.discharge(container, predator.getRunningLength());
        container.fastaddStreamOfEmptyWords(false, predator.getRunningLength() -
                                                       index);
      }
      predator.discardRunningWordsWithReload();
    }
    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k) {
        container.addWord(rlwi.getLiteralWordAt(k) & rlwj.getLiteralWordAt(k));
      }
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  BufferedRunningLengthWord<uword> &remain = rlwj.size() > 0 ? rlwj : rlwi;
  while(remain.size() > 0) {
    container.addStreamOfEmptyWords(false, remain.size());
    if (!remain.next()) { break; }
  }
  container.setSizeInBits(sizeInBits() > a.sizeInBits() ? sizeInBits()
                                                        : a.sizeInBits());
  container.assertWordCount("logicaland");
}

template <class uword>
void EWAHBoolArray<uword>::logicalandnot(const EWAHBoolArray &a,
                                         EWAHBoolArray &container) const {
  container.reset();
  if (RESERVEMEMORY)
    container.buffer.reserve(buffer.size() > a.buffer.size() ? buffer.size()
                                                             : a.buffer.size());
  EWAHBoolArrayRawIterator<uword> i = raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = a.raw_iterator();
  if (!j.hasNext()) {  // the other fellow is empty
    container = *this; // just copy, stupidly, the data
    return;
  }
  if (!(i.hasNext())) { // hopefully this never happens...
    container.setSizeInBits(sizeInBits());
    return;
  }
  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();

  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey(i_is_prey ? rlwi : rlwj);
      BufferedRunningLengthWord<uword> &predator(i_is_prey ? rlwj : rlwi);
      if (((predator.getRunningBit()) && (i_is_prey)) ||
          ((!predator.getRunningBit()) && (!i_is_prey))) {
        container.fastaddStreamOfEmptyWords(false, predator.getRunningLength());
        prey.discardFirstWordsWithReload(predator.getRunningLength());
      } else if (i_is_prey) {
        const size_t index =
            prey.discharge(container, predator.getRunningLength());
        container.fastaddStreamOfEmptyWords(false, predator.getRunningLength() -
                                                       index);
      } else {
        const size_t index =
            prey.dischargeNegated(container, predator.getRunningLength());
        container.fastaddStreamOfEmptyWords(true, predator.getRunningLength() -
                                                      index);
      }
      predator.discardRunningWordsWithReload();
    }
    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k) {
        container.addWord(static_cast<uword>(rlwi.getLiteralWordAt(k) & ~rlwj.getLiteralWordAt(k)));
      }
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  if(rlwi.size() > 0) {
    rlwi.discharge(container);
    container.setSizeInBits(sizeInBits());
  } else {
    while(rlwj.size() > 0) {
      container.addStreamOfEmptyWords(false, rlwj.size());
      if (!rlwj.next()) { break; }
    }
    container.setSizeInBits(a.sizeInBits());
  }
  container.assertWordCount("logicalandnot");
}

template <class uword>
size_t EWAHBoolArray<uword>::logicalandnotcount(const EWAHBoolArray &a) const {
  EWAHBoolArrayRawIterator<uword> i = raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = a.raw_iterator();
  if (!j.hasNext()) { // the other fellow is empty
    return this->numberOfOnes();
  }
  if (!(i.hasNext())) { // hopefully this never happens...
    return 0;
  }
  size_t answer = 0;
  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();

  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey(i_is_prey ? rlwi : rlwj);
      BufferedRunningLengthWord<uword> &predator(i_is_prey ? rlwj : rlwi);
      if (((predator.getRunningBit()) && (i_is_prey)) ||
          ((!predator.getRunningBit()) && (!i_is_prey))) {
        prey.discardFirstWordsWithReload(predator.getRunningLength());
      } else if (i_is_prey) {
        prey.dischargeCount(predator.getRunningLength(), &answer);
      } else {
        const size_t index =
            prey.dischargeCountNegated(predator.getRunningLength(), &answer);
        answer += (predator.getRunningLength() - index) * wordinbits;
      }
      predator.discardRunningWordsWithReload();
    }
    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k) {
        answer += countOnes(
            (uword)(rlwi.getLiteralWordAt(k) & (~rlwj.getLiteralWordAt(k))));
      }
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  const bool i_remains = rlwi.size() > 0;
  if (i_remains) {
    answer += rlwi.dischargeCount();
  }
  return answer;
}

template <class uword>
size_t EWAHBoolArray<uword>::logicalandcount(const EWAHBoolArray &a) const {
  EWAHBoolArrayRawIterator<uword> i = a.raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = raw_iterator();
  if (!(i.hasNext() and j.hasNext())) { // hopefully this never happens...
    return 0;
  }
  size_t answer = 0;
  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();

  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey(i_is_prey ? rlwi : rlwj);
      BufferedRunningLengthWord<uword> &predator(i_is_prey ? rlwj : rlwi);
      if (!predator.getRunningBit()) {
        prey.discardFirstWordsWithReload(predator.getRunningLength());
      } else {
        // const size_t index =
        prey.dischargeCount(predator.getRunningLength(), &answer);
      }
      predator.discardRunningWordsWithReload();
    }
    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k) {
        answer += countOnes(
            (uword)(rlwi.getLiteralWordAt(k) & rlwj.getLiteralWordAt(k)));
      }
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  return answer;
}

template <class uword>
bool EWAHBoolArray<uword>::intersects(const EWAHBoolArray &a) const {
  EWAHBoolArrayRawIterator<uword> i = a.raw_iterator();
  EWAHBoolArrayRawIterator<uword> j = raw_iterator();
  if (!(i.hasNext() and j.hasNext())) { // hopefully this never happens...
    return false;
  }
  // at this point, this should be safe:
  BufferedRunningLengthWord<uword> &rlwi = i.next();
  BufferedRunningLengthWord<uword> &rlwj = j.next();

  while ((rlwi.size() > 0) && (rlwj.size() > 0)) {
    while ((rlwi.getRunningLength() > 0) || (rlwj.getRunningLength() > 0)) {
      const bool i_is_prey = rlwi.getRunningLength() < rlwj.getRunningLength();
      BufferedRunningLengthWord<uword> &prey(i_is_prey ? rlwi : rlwj);
      BufferedRunningLengthWord<uword> &predator(i_is_prey ? rlwj : rlwi);
      if (!predator.getRunningBit()) {
        prey.discardFirstWordsWithReload(predator.getRunningLength());
      } else {
        size_t index = 0;
        bool isnonzero =
            prey.nonzero_discharge(predator.getRunningLength(), index);
        if (isnonzero)
          return true;
      }
      predator.discardRunningWordsWithReload();
    }
    const uword nbre_literal = std::min(rlwi.getNumberOfLiteralWords(),
                                         rlwj.getNumberOfLiteralWords());
    if (nbre_literal > 0) {
      for (size_t k = 0; k < nbre_literal; ++k) {
        if ((rlwi.getLiteralWordAt(k) & rlwj.getLiteralWordAt(k)) != 0)
          return true;
      }
      rlwi.discardLiteralWordsWithReload(nbre_literal);
      rlwj.discardLiteralWordsWithReload(nbre_literal);
    }
  }
  return false;
}

template <class uword>
BitmapStatistics EWAHBoolArray<uword>::computeStatistics() const {
  BitmapStatistics bs;
  EWAHBoolArrayRawIterator<uword> i = raw_iterator();
  while (i.hasNext()) {
    BufferedRunningLengthWord<uword> &brlw(i.next());
    ++bs.runningwordmarker;
    bs.totalliteral += brlw.getNumberOfLiteralWords();
    bs.totalcompressed += brlw.getRunningLength();
    if (brlw.getRunningLength() ==
        RunningLengthWord<uword>::largestrunninglengthcount) {
      ++bs.maximumofrunningcounterreached;
    }
  }
  return bs;
}

template <class uword> void EWAHBoolArray<uword>::debugprintout() const {
  std::cout << "==printing out EWAHBoolArray==" << std::endl;
  std::cout << "Number of compressed words: " << buffer.size() << std::endl;
  std::cout << "Size in bits: " << sizeinbits << std::endl;

  size_t pointer = 0;
  while (pointer < buffer.size()) {
    ConstRunningLengthWord<uword> rlw(buffer[pointer]);
    bool b = rlw.getRunningBit();
    const uword rl = rlw.getRunningLength();
    const uword lw = rlw.getNumberOfLiteralWords();
    std::cout << "pointer = " << pointer << " running bit=" << b
              << " running length=" << rl << " lit. words=" << lw << std::endl;
    for (uword j = 0; j < lw; ++j) {
      const uword &w = buffer[pointer + j + 1];
      std::cout << toBinaryString(w) << std::endl;
    }
    pointer += lw + 1;
  }
  std::cout << "==END==" << std::endl;
}

template <class uword>
size_t EWAHBoolArray<uword>::sizeOnDisk(const bool savesizeinbits) const {
  return (savesizeinbits ? sizeof(uint64_t) : 0) + sizeof(uint64_t) +
         sizeof(uword) * buffer.size();
}
} // namespace ewah
#endif
