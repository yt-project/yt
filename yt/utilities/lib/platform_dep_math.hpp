/*
This file provides a compatibility layout between MSVC, and different version of GCC.

MSVC does not define isnormal in the std:: namespace, so we cannot import it from <cmath>, but from <math.h> instead.
However with GCC-5, there is a clash between the definition of isnormal in <math.h> and using C++14, so we need to import from cmath instead.
*/

#if _MSC_VER
#include <math.h>
inline bool __isnormal(double x) {
    return isnormal(x);
}
#else
#include <cmath>
inline bool __isnormal(double x) {
    return std::isnormal(x);
}
#endif
