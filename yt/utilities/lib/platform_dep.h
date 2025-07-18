#include <math.h>
#ifdef MS_WIN32
#include "malloc.h"
/*
note: the following implicitly sets a mininum VS version: conservative
minimum is _MSC_VER >= 1928 (VS 2019, 16.8), but may work for VS 2015
but that has not been tested. see https://github.com/yt-project/yt/pull/4980
and https://learn.microsoft.com/en-us/cpp/overview/visual-cpp-language-conformance
*/
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
