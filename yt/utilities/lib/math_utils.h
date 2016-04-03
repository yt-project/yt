/* Taken from http://siliconandlithium.blogspot.com/2014/05/msvc-c99-mathh-header.html */
#include <math.h>
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


