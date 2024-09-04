#include <math.h>
#ifdef MS_WIN32
#include "malloc.h"
// note: the following implicitly requires _MSC_VER >= 1928 (VS 2019, 16.8)
#include <float.h>
#include <stdint.h>
#elif defined(__FreeBSD__)
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#else
#include <stdint.h>
#include "alloca.h"
#include <math.h>
#endif
