/**
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 *
 * (c) Daniel Lemire, http://lemire.me/en/
 */
#ifndef RUNNINGLENGTHWORD_H_
#define RUNNINGLENGTHWORD_H_
#include <iostream>
namespace ewah {

/**
 * For expert users.
 * This class is used to represent a special type of word storing
 * a run length. It is defined by the Enhanced Word Aligned  Hybrid (EWAH)
 * format. You don't normally need to access this class.
 */
template <class uword> class RunningLengthWord {
public:
  RunningLengthWord(uword &data) : mydata(data) {}

  RunningLengthWord(const RunningLengthWord &rlw) : mydata(rlw.mydata) {}

  RunningLengthWord &operator=(const RunningLengthWord &rlw) {
    mydata = rlw.mydata;
    return *this;
  }

  /**
   * Which bit is being repeated?
   */
  bool getRunningBit() const { return mydata & static_cast<uword>(1); }

  /**
   * how many words should be filled by the running bit
   */
  static inline bool getRunningBit(uword data) {
    return data & static_cast<uword>(1);
  }

  /**
   * how many words should be filled by the running bit
   */
  uword getRunningLength() const {
    return static_cast<uword>((mydata >> 1) & largestrunninglengthcount);
  }

  /**
   * followed by how many literal words?
   */
  static inline uword getRunningLength(uword data) {
    return static_cast<uword>((data >> 1) & largestrunninglengthcount);
  }

  /**
   * followed by how many literal words?
   */
  uword getNumberOfLiteralWords() const {
    return static_cast<uword>(mydata >> (1 + runninglengthbits));
  }

  /**
   * Total of getRunningLength() and getNumberOfLiteralWords()
   */
  uword size() const {
    return static_cast<uword>(getRunningLength() + getNumberOfLiteralWords());
  }

  /**
   * Total of getRunningLength() and getNumberOfLiteralWords()
   */
  static inline uword size(uword data) {
    return static_cast<uword>(getRunningLength(data) +
                              getNumberOfLiteralWords(data));
  }

  /**
   * followed by how many literal words?
   */
  static inline uword getNumberOfLiteralWords(uword data) {
    return static_cast<uword>(data >> (1 + runninglengthbits));
  }

  /**
   * running length of which type of bits
   */
  void setRunningBit(bool b) {
    if (b)
      mydata |= static_cast<uword>(1);
    else
      mydata &= static_cast<uword>(~1);
  }

  void discardFirstWords(uword x) {
    const uword rl(getRunningLength());
    if (rl >= x) {
      setRunningLength(rl - x);
      return;
    }
    x -= rl;
    setRunningLength(0);
    setNumberOfLiteralWords(getNumberOfLiteralWords() - x);
  }

  /**
   * running length of which type of bits
   */
  static inline void setRunningBit(uword &data, bool b) {
    if (b)
      data |= static_cast<uword>(1);
    else
      data &= static_cast<uword>(~1);
  }

  void setRunningLength(uword l) {
    mydata |= shiftedlargestrunninglengthcount;
    mydata &=
        static_cast<uword>((l << 1) | notshiftedlargestrunninglengthcount);
  }

  // static call for people who hate objects
  static inline void setRunningLength(uword &data, uword l) {
    data |= shiftedlargestrunninglengthcount;
    data &= static_cast<uword>((l << 1) | notshiftedlargestrunninglengthcount);
  }

  void setNumberOfLiteralWords(uword l) {
    mydata |= notrunninglengthplusrunningbit;
    mydata &= static_cast<uword>((l << (runninglengthbits + 1)) |
                                 runninglengthplusrunningbit);
  }
  // static call for people who hate objects
  static inline void setNumberOfLiteralWords(uword &data, uword l) {
    data |= notrunninglengthplusrunningbit;
    data &= static_cast<uword>(l << (runninglengthbits + 1)) |
            runninglengthplusrunningbit;
  }

  static const uint32_t runninglengthbits = sizeof(uword) * 4;
  static const uint32_t literalbits = sizeof(uword) * 8 - 1 - runninglengthbits;
  static const uword largestliteralcount =
      (static_cast<uword>(1) << literalbits) - 1;
  static const uword largestrunninglengthcount =
      (static_cast<uword>(1) << runninglengthbits) - 1;
  static const uword shiftedlargestrunninglengthcount =
      largestrunninglengthcount << 1;
  static const uword notshiftedlargestrunninglengthcount =
      static_cast<uword>(~shiftedlargestrunninglengthcount);
  static const uword runninglengthplusrunningbit =
      (static_cast<uword>(1) << (runninglengthbits + 1)) - 1;
  static const uword notrunninglengthplusrunningbit =
      static_cast<uword>(~runninglengthplusrunningbit);
  static const uword notlargestrunninglengthcount =
      static_cast<uword>(~largestrunninglengthcount);

  uword &mydata;
};

/**
 * Same as RunningLengthWord, except that the values cannot be modified.
 */
template <class uword = uint32_t> class ConstRunningLengthWord {
public:
  ConstRunningLengthWord() : mydata(0) {}

  ConstRunningLengthWord(const uword data) : mydata(data) {}

  ConstRunningLengthWord(const ConstRunningLengthWord &rlw)
      : mydata(rlw.mydata) {}

  /**
   * Which bit is being repeated?
   */
  bool getRunningBit() const { return mydata & static_cast<uword>(1); }

  /**
   * how many words should be filled by the running bit
   */
  uword getRunningLength() const {
    return static_cast<uword>(
        (mydata >> 1) & RunningLengthWord<uword>::largestrunninglengthcount);
  }

  /**
   * followed by how many literal words?
   */
  uword getNumberOfLiteralWords() const {
    return static_cast<uword>(
        mydata >> (1 + RunningLengthWord<uword>::runninglengthbits));
  }

  /**
   * Total of getRunningLength() and getNumberOfLiteralWords()
   */
  uword size() const { return getRunningLength() + getNumberOfLiteralWords(); }

  uword mydata;
};

template <class uword> class EWAHBoolArray;

template <class uword> class EWAHBoolArrayRawIterator;

/**
 * Same as RunningLengthWord, except that the values are buffered for quick
 * access.
 */
template <class uword = uint32_t> class BufferedRunningLengthWord {
public:
  enum { wordinbits = sizeof(uword) * 8 };

  BufferedRunningLengthWord(const uword &data,
                            EWAHBoolArrayRawIterator<uword> *p)
      : RunningBit(data & static_cast<uword>(1)),
        RunningLength(static_cast<uword>(
            (data >> 1) & RunningLengthWord<uword>::largestrunninglengthcount)),
        NumberOfLiteralWords(static_cast<uword>(
            data >> (1 + RunningLengthWord<uword>::runninglengthbits))),
        parent(p) {}
  BufferedRunningLengthWord(const RunningLengthWord<uword> &p)
      : RunningBit(p.mydata & static_cast<uword>(1)),
        RunningLength((p.mydata >> 1) &
                      RunningLengthWord<uword>::largestrunninglengthcount),
        NumberOfLiteralWords(p.mydata >>
                             (1 + RunningLengthWord<uword>::runninglengthbits)),
        parent(p.parent) {}

  void discharge(EWAHBoolArray<uword> &container) {
    while (size() > 0) {
      // first run
      size_t pl = getRunningLength();
      container.fastaddStreamOfEmptyWords(getRunningBit(), pl);
      size_t pd = getNumberOfLiteralWords();
      writeLiteralWords(pd, container);
      if (!next())
        break;
    }
  }

  size_t dischargeCount() {
    size_t answer = 0;
    while (size() > 0) {
      // first run
      if (getRunningBit()) {
        answer += wordinbits * getRunningLength();
      }
      size_t pd = getNumberOfLiteralWords();
      for (size_t i = 0; i < pd; ++i)
        answer += countOnes((uword)getLiteralWordAt(i));
      if (!next())
        break;
    }
    return answer;
  }

  size_t dischargeCountNegated() {
    size_t answer = 0;
    while (size() > 0) {
      // first run
      if (!getRunningBit()) {
        answer += wordinbits * getRunningLength();
      }
      size_t pd = getNumberOfLiteralWords();
      for (size_t i = 0; i < pd; ++i)
        answer += countOnes((uword)(~getLiteralWordAt(i)));
      if (!next())
        break;
    }
    return answer;
  }

  // Symbolically write out up to max words, returns how many were written,
  // write to count the number bits written (we assume that count was initially
  // zero)
  size_t dischargeCount(size_t max, size_t *count) {
    size_t index = 0;
    while (true) {
      if (index + RunningLength > max) {
        const size_t offset = max - index;
        if (getRunningBit())
          *count += offset * wordinbits;
        RunningLength -= offset;
        return max;
      }
      if (getRunningBit())
        *count += RunningLength * wordinbits;
      index += RunningLength;
      if (NumberOfLiteralWords + index > max) {
        const size_t offset = max - index;
        for (size_t i = 0; i < offset; ++i)
          *count += countOnes((uword)getLiteralWordAt(i));
        RunningLength = 0;
        NumberOfLiteralWords -= offset;
        return max;
      }
      for (size_t i = 0; i < NumberOfLiteralWords; ++i)
        *count += countOnes((uword)getLiteralWordAt(i));
      index += NumberOfLiteralWords;
      if (!next())
        break;
    }
    return index;
  }

  size_t dischargeCountNegated(size_t max, size_t *count) {
    size_t index = 0;
    while (true) {
      if (index + RunningLength > max) {
        const size_t offset = max - index;
        if (!getRunningBit())
          *count += offset * wordinbits;
        RunningLength -= offset;
        return max;
      }
      if (!getRunningBit())
        *count += RunningLength * wordinbits;
      index += RunningLength;
      if (NumberOfLiteralWords + index > max) {
        const size_t offset = max - index;
        for (size_t i = 0; i < offset; ++i)
          *count += countOnes((uword)(~getLiteralWordAt(i)));
        RunningLength = 0;
        NumberOfLiteralWords -= offset;
        return max;
      }
      for (size_t i = 0; i < NumberOfLiteralWords; ++i)
        *count += countOnes((uword)(~getLiteralWordAt(i)));
      index += NumberOfLiteralWords;
      if (!next())
        break;
    }
    return index;
  }
  bool nonzero_discharge() {
    while (size() > 0) {
      // first run
      size_t pl = getRunningLength();
      if ((pl > 0) && (getRunningBit()))
        return true;
      size_t pd = getNumberOfLiteralWords();
      if (pd > 0)
        return true;
      discardFirstWordsWithReload(static_cast<uword>(pl + pd));
    }
    return false;
  }

  // Write out up to max words, returns how many were written
  size_t discharge(EWAHBoolArray<uword> &container, size_t max) {
    size_t index = 0;
    while (true) {
      if (index + RunningLength > max) {
        const size_t offset = max - index;
        container.fastaddStreamOfEmptyWords(getRunningBit(), offset);
        RunningLength = static_cast<uword>(RunningLength - offset);
        return max;
      }
      container.fastaddStreamOfEmptyWords(getRunningBit(), RunningLength);
      index += RunningLength;
      if (NumberOfLiteralWords + index > max) {
        const size_t offset = max - index;
        writeLiteralWords(offset, container);
        RunningLength = 0;
        NumberOfLiteralWords = static_cast<uword>(NumberOfLiteralWords - offset);
        return max;
      }
      writeLiteralWords(NumberOfLiteralWords, container);
      index += NumberOfLiteralWords;
      if (!next())
        break;
    }
    return index;
  }

  bool nonzero_discharge(size_t max, size_t &index) {
    index = 0;
    while ((index < max) && (size() > 0)) {
      // first run
      size_t pl = getRunningLength();
      if (index + pl > max) {
        pl = max - index;
      }
      if ((getRunningBit()) && (pl > 0))
        return true;
      index += pl;
      size_t pd = getNumberOfLiteralWords();
      if (pd + index > max) {
        pd = max - index;
      }
      if (pd > 0)
        return true;
      discardFirstWordsWithReload(static_cast<uword>(pl + pd));
    }
    return false;
  }

  // Write out up to max words, returns how many were written
  size_t dischargeNegated(EWAHBoolArray<uword> &container, size_t max) {
    // todo: could be optimized further
    size_t index = 0;
    while ((index < max) && (size() > 0)) {
      // first run
      size_t pl = getRunningLength();
      if (index + pl > max) {
        pl = max - index;
      }
      container.fastaddStreamOfEmptyWords(!getRunningBit(), pl);
      index += pl;
      size_t pd = getNumberOfLiteralWords();
      if (pd + index > max) {
        pd = max - index;
      }
      writeNegatedLiteralWords(pd, container);
      discardFirstWordsWithReload(static_cast<uword>(pl + pd));
      index += pd;
    }
    return index;
  }
  bool nonzero_dischargeNegated(size_t max, size_t &index) {
    while ((index < max) && (size() > 0)) {
      // first run
      size_t pl = getRunningLength();
      if (index + pl > max) {
        pl = max - index;
      }
      if ((!getRunningBit()) && (pl > 0))
        return true;
      index += pl;
      size_t pd = getNumberOfLiteralWords();
      if (pd + index > max) {
        pd = max - index;
      }
      if (pd > 0)
        return true;
      discardFirstWordsWithReload(static_cast<uword>(pl + pd));
      index += pd;
    }
    return false;
  }

  uword getLiteralWordAt(size_t index) { return parent->dirtyWords()[index]; }

  void writeLiteralWords(size_t numWords, EWAHBoolArray<uword> &container) {
    container.fastaddStreamOfDirtyWords(parent->dirtyWords(), numWords);
  }

  void writeNegatedLiteralWords(size_t numWords,
                                EWAHBoolArray<uword> &container) {
    container.addStreamOfNegatedDirtyWords(parent->dirtyWords(), numWords);
  }

  void discardRunningWords() { RunningLength = 0; }

  void discardRunningWordsWithReload() {
    RunningLength = 0;
    if (NumberOfLiteralWords == 0)
      next();
  }

  bool next() {
    if (!parent->hasNext()) {
      NumberOfLiteralWords = 0;
      RunningLength = 0;
      return false;
    }
    parent->next();
    return true;
  }

  void read(const uword &data) {
    RunningBit = data & static_cast<uword>(1);
    RunningLength = static_cast<uword>(
        (data >> 1) & RunningLengthWord<uword>::largestrunninglengthcount);
    NumberOfLiteralWords = static_cast<uword>(
        data >> (1 + RunningLengthWord<uword>::runninglengthbits));
  }

  /**
   * Which bit is being repeated?
   */
  bool getRunningBit() const { return RunningBit; }

  void discardFirstWords(uword x) {
    if (RunningLength >= x) {
      RunningLength = static_cast<uword>(RunningLength - x);
      return;
    }
    x = static_cast<uword>(x - RunningLength);
    RunningLength = 0;
    NumberOfLiteralWords = static_cast<uword>(NumberOfLiteralWords - x);
  }

  /**
   * how many words should be filled by the running bit (see previous method)
   */
  uword getRunningLength() const { return RunningLength; }

  /**
   * followed by how many literal words?
   */
  uword getNumberOfLiteralWords() const { return NumberOfLiteralWords; }

  /**
   * Total of getRunningLength() and getNumberOfLiteralWords()
   */
  uword size() const {
    return static_cast<uword>(RunningLength + NumberOfLiteralWords);
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  const BufferedRunningLengthWord &a) {
    out << "{RunningBit:" << a.RunningBit
        << ",RunningLength:" << a.RunningLength
        << ",NumberOfLiteralWords:" << a.NumberOfLiteralWords << "}";
    return out;
  }
  void discardLiteralWordsWithReload(uword x) {
    assert(NumberOfLiteralWords >= x);
    NumberOfLiteralWords -= x;
    if (NumberOfLiteralWords == 0)
      next();
  }

  void discardFirstWordsWithReload(uword x) {
    while (x > 0) {
      if (RunningLength > x) {
        RunningLength = static_cast<uword>(RunningLength - x);
        return;
      }
      x = static_cast<uword>(x - RunningLength);
      RunningLength = 0;
      size_t toDiscard = x > NumberOfLiteralWords ? NumberOfLiteralWords : x;
      NumberOfLiteralWords =
          static_cast<uword>(NumberOfLiteralWords - toDiscard);
      x = static_cast<uword>(x - toDiscard);
      if ((x > 0) || (size() == 0)) {
        if (!next())
          break;
      }
    }
  }

private:
  bool RunningBit;
  uword RunningLength;
  uword NumberOfLiteralWords;
  EWAHBoolArrayRawIterator<uword> *parent;
};
} // namespace ewah

#endif /* RUNNINGLENGTHWORD_H_ */
