#ifndef __COSMOLOGY_H__
#define __COSMOLOGY_H__

typedef struct CosmologyParametersStruct CosmologyParameters;

#define COSMOLOGY_DECLARE_PRIMARY_PARAMETER(name) \
void cosmology_set_##name(CosmologyParameters *c, double value)

#define cosmology_set(c,name,value)	\
cosmology_set_##name(c,value)

COSMOLOGY_DECLARE_PRIMARY_PARAMETER(OmegaM);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(OmegaB);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(OmegaL);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(h);
COSMOLOGY_DECLARE_PRIMARY_PARAMETER(DeltaDC);

#undef COSMOLOGY_DECLARE_PRIMARY_PARAMETER

CosmologyParameters *cosmology_allocate();
void cosmology_free(CosmologyParameters *c);

/*
//  Check that all required cosmological parameters have been set.
//  The minimum set is OmegaM, OmegaB, and h. By default, zero OmegaL,
//  OmegaK, and the DC mode are assumed.
*/
int cosmology_is_set(CosmologyParameters *c);


/*
//  Freeze the cosmology and forbid any further changes to it.
//  In codes that include user-customizable segments (like plugins),
//  this function van be used for insuring that a user does not
//  change the cosmology in mid-run.
*/
void cosmology_set_fixed(CosmologyParameters *c);


/*
//  Manual initialization. This does not need to be called, 
//  the initialization is done automatically on the first call
//  to a relevant function.
*/
void cosmology_init(CosmologyParameters *c);


/*
//  Set the range of global scale factors for thread-safe
//  calls to direct functions until the argument leaves the range.
*/
void cosmology_set_thread_safe_range(CosmologyParameters *c, double amin, double amax);

/*
//  Direct functions take the global cosmological scale factor as the argument.
//  These functionsare are thread-safe if called with the argument in the
//  range set by a prior call to cosmology_set_thread_safe_range(...).
//  Calling them with the argument outside that range is ok, but breaks
//  thread-safety assurance.
*/

#define DEFINE_FUN(name)         \
double name(CosmologyParameters *c, double a); \
double inv_##name(CosmologyParameters *c, double v);

DEFINE_FUN(aBox);
DEFINE_FUN(tCode);
DEFINE_FUN(tPhys);
DEFINE_FUN(dPlus);
DEFINE_FUN(qPlus); /* Q+ = a^2 dD+/(H0 dt) */

#undef DEFINE_FUN

/*
//  Conversion macros
*/
#define  abox_from_auni(c,a)   aBox(c,a)
#define tcode_from_auni(c,a)  tCode(c,a)
#define tphys_from_auni(c,a)  tPhys(c,a)
#define dplus_from_auni(c,a)  dPlus(c,a)

#define auni_from_abox(c,v)   inv_aBox(c,v)
#define auni_from_tcode(c,v)  inv_tCode(c,v)
#define auni_from_tphys(c,v)  inv_tPhys(c,v)
#define auni_from_dplus(c,v)  inv_dPlus(c,v)

#define abox_from_tcode(c,tcode)   aBox(c,inv_tCode(c,tcode))
#define tcode_from_abox(c,abox)    tCode(c,inv_aBox(c,abox))

#define tphys_from_abox(c,abox)    tPhys(c,inv_aBox(c,abox))
#define tphys_from_tcode(c,tcode)  tPhys(c,inv_tCode(c,tcode))
#define dplus_from_tcode(c,tcode)  dPlus(c,inv_tCode(c,tcode))

/*
//  Hubble parameter in km/s/Mpc; defined as macro so that it can be
//  undefined if needed to avoid the name clash.
*/
double cosmology_mu(CosmologyParameters *c, double a);
#define Hubble(c,a) (100*c->h*cosmology_mu(c,a)/(a*a))

#endif /* __COSMOLOGY_H__ */
