#if defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#include "malloc.h"
#else
#include "alloca.h"
#endif