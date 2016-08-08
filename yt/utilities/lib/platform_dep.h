#include <math.h>
#ifdef MS_WIN32
#include "malloc.h"
typedef int int32_t;
typedef long long int64_t;
/* Taken from http://siliconandlithium.blogspot.com/2014/05/msvc-c99-mathh-header.html */
#define isnormal(x) ((_fpclass(x) == _FPCLASS_NN) || (_fpclass(x) == _FPCLASS_PN))
static __inline double rint(double x){
    const double two_to_52 = 4.5035996273704960e+15;
    double fa = fabs(x);
    if(fa >= two_to_52){
        return x;
    } else{
        return copysign(two_to_52 + fa - two_to_52, x);
    }
}
static __inline long int lrint(double x){
    return (long)rint(x);
}
static __inline double fmax(double x, double y){
    return (x > y) ? x : y;
}
static __inline double fmin(double x, double y){
    return (x < y) ? x : y;
}
static __inline double log2(double x) {
    return log(x) * M_LOG2E;
}

/* adapted from http://www.johndcook.com/blog/cpp_erf/
   code is under public domain license */

double erf(double x)
{
    /* constants */
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    double t;
    double y;

    /* Save the sign of x */
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x);

    /* A&S formula 7.1.26 */
    t = 1.0/(1.0 + p*x);
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return sign*y;
}

#else
#include <stdint.h>
#include "alloca.h"
#include <math.h>
#endif


