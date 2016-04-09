#include <math.h>
#include <stdio.h>
#include <string.h>


#ifndef ERROR
#include <stdio.h>
#define ERROR(msg) { fprintf(stderr,"%s\n",msg); exit(1); }
#endif

#ifndef ASSERT
#include <stdio.h>
#define ASSERT(exp) { if(!(exp)) { fprintf(stderr,"Failed assertion %s, line: %d\n",#exp,__LINE__); } }
#endif

#ifndef HEAPALLOC
#include <stdlib.h>
#define HEAPALLOC(type,size)	(type *)malloc((size)*sizeof(type))
#endif

#ifndef NEWARR
#include <stdlib.h>
#define NEWARR(size)   HEAPALLOC(double,size)
#endif

#ifndef DELETE
#include <stdlib.h>
#define DELETE(ptr)    free(ptr)
#endif


#include "cosmology.h"


struct CosmologyParametersStruct
{
  int set;
  int ndex;
  int size;
  double *la;
  double *aUni;
  double *aBox;
  double *tCode;
  double *tPhys;
  double *dPlus;
  double *qPlus;
  double aLow;
  double tCodeOffset;

  double OmegaM;
  double OmegaD;
  double OmegaB;
  double OmegaL;
  double OmegaK;
  double OmegaR;
  double h;
  double DeltaDC;
  int flat;
  double Omh2;
  double Obh2;
};

void cosmology_clear_table(CosmologyParameters *c);
void cosmology_fill_table(CosmologyParameters *c, double amin, double amax);
void cosmology_fill_table_abox(CosmologyParameters *c, int istart, int n);

CosmologyParameters *cosmology_allocate() {
	CosmologyParameters *c = HEAPALLOC(CosmologyParameters,1);
	if ( c != NULL ) {
		memset(c, 0, sizeof(CosmologyParameters));

		c->ndex = 200;
		c->aLow = 1.0e-2;
	}
	return c;
}

void cosmology_free(CosmologyParameters *c) {
	cosmology_clear_table(c);
	DELETE(c);
}

int cosmology_is_set(CosmologyParameters *c)
{
  return (c->OmegaM>0.0 && c->OmegaB>0.0 && c->h>0.0);
}


void cosmology_fail_on_reset(const char *name, double old_value, double new_value)
{
  char str[150];
  sprintf(str,"Trying to change %s from %lg to %lg...\nCosmology has been fixed and cannot be changed.\n",name,old_value,new_value);
  ERROR(str);
}


void cosmology_set_OmegaM(CosmologyParameters *c, double v)
{
  if(v < 1.0e-3) v = 1.0e-3;
  if(fabs(c->OmegaM-v) > 1.0e-5)
    {
      if(c->set) cosmology_fail_on_reset("OmegaM",c->OmegaM,v);
      c->OmegaM = v;
      c->flat = (fabs(c->OmegaM+c->OmegaL-1.0) > 1.0e-5) ? 0 : 1;
      cosmology_clear_table(c);
    }
}


void cosmology_set_OmegaL(CosmologyParameters *c, double v)
{
  if(fabs(c->OmegaL-v) > 1.0e-5)
    {
      if(c->set) cosmology_fail_on_reset("OmegaL",c->OmegaL,v);
      c->OmegaL = v;
      c->flat = (fabs(c->OmegaM+c->OmegaL-1.0) > 1.0e-5) ? 0 : 1;
      cosmology_clear_table(c);
    }
}


void cosmology_set_OmegaB(CosmologyParameters *c, double v)
{
  if(v < 0.0) v = 0.0;
  if(fabs(c->OmegaB-v) > 1.0e-5)
    {
      if(c->set) cosmology_fail_on_reset("OmegaB",c->OmegaB,v);
      c->OmegaB = v;
      cosmology_clear_table(c);
    }
}


void cosmology_set_h(CosmologyParameters *c, double v)
{
  if(fabs(c->h-v) > 1.0e-5)
    {
      if(c->set) cosmology_fail_on_reset("h",c->h,v);
      c->h = v;
      cosmology_clear_table(c);
    }
}


void cosmology_set_DeltaDC(CosmologyParameters *c, double v)
{
  if(fabs(c->DeltaDC-v) > 1.0e-3)
    {
      if(c->set) cosmology_fail_on_reset("DeltaDC",c->DeltaDC,v);
      c->DeltaDC = v;
      cosmology_clear_table(c);
    }
}


void cosmology_init(CosmologyParameters *c)
{
  if(c->size == 0) /* reset only if the state is dirty */
    {
      if(!cosmology_is_set(c)) ERROR("Not all of the required cosmological parameters have been set; the minimum required set is (OmegaM,OmegaB,h).");

      if(c->OmegaB > c->OmegaM) c->OmegaB = c->OmegaM;
      c->OmegaD = c->OmegaM - c->OmegaB;
      if(c->flat)
	{
	  c->OmegaK = 0.0;
	  c->OmegaL = 1.0 - c->OmegaM;
	}
      else
	{
	  c->OmegaK = 1.0 - (c->OmegaM+c->OmegaL);
	}
      c->OmegaR = 4.166e-5/(c->h*c->h);

      c->Omh2 = c->OmegaM*c->h*c->h;
      c->Obh2 = c->OmegaB*c->h*c->h;

      cosmology_fill_table(c,c->aLow,1.0);

      c->tCodeOffset = 0.0;  /*  Do need to set it to zero first */
#ifndef NATIVE_TCODE_NORMALIZATION
      c->tCodeOffset = 0.0 - tCode(c,inv_aBox(c,1.0));
#endif
    }      
}


void cosmology_set_fixed(CosmologyParameters *c)
{
  cosmology_init(c);
  c->set = 1;
}


double cosmology_mu(CosmologyParameters *c, double a)
{
  return sqrt(((a*a*c->OmegaL+c->OmegaK)*a+c->OmegaM)*a+c->OmegaR);
}


double cosmology_dc_factor(CosmologyParameters *c, double dPlus)
{
  double dc = 1.0 + dPlus*c->DeltaDC;
  return 1.0/pow((dc>0.001)?dc:0.001,1.0/3.0);
}


void cosmology_fill_table_integrate(CosmologyParameters *c, double a, double y[], double f[])
{
  double mu = cosmology_mu(c, a);
  double abox = a*cosmology_dc_factor(c, y[2]);
  
  f[0] = a/(abox*abox*mu);
  f[1] = a/mu;
  f[2] = y[3]/(a*mu);
  f[3] = 1.5*c->OmegaM*y[2]/mu;
}

#ifdef _WIN32
double asinh(double x){
    return log(x + sqrt((x * x) + 1.0));
}
#endif

void cosmology_fill_table_piece(CosmologyParameters *c, int istart, int n)
{
  int i, j;
  double tPhysUnit = (3.0856775813e17/(365.25*86400))/c->h;  /* 1/H0 in Julian years */

  double x, aeq = c->OmegaR/c->OmegaM;
  double tCodeFac = 1.0/sqrt(aeq);
  double tPhysFac = tPhysUnit*aeq*sqrt(aeq)/sqrt(c->OmegaM);

  double da, a0, y0[4], y1[4];
  double f1[4], f2[4], f3[4], f4[4];

  for(i=istart; i<n; i++)
    {
      c->aUni[i] = pow(10.0,c->la[i]);
    }

  /*
  //  Small a regime, use analytical formulae for matter + radiation model
  */  
  for(i=istart; c->aUni[i]<(c->aLow+1.0e-9) && i<n; i++)
    {
      x = c->aUni[i]/aeq;

      c->tPhys[i] = tPhysFac*2*x*x*(2+sqrt(x+1))/(3*pow(1+sqrt(x+1),2.0));
      c->dPlus[i] = aeq*(x + 2.0/3.0 + (6*sqrt(1+x)+(2+3*x)*log(x)-2*(2+3*x)*log(1+sqrt(1+x)))/(log(64.0)-9));  /* long last term is the decaying mode generated after euality; it is very small for x > 10, I keep ot just for completeness; */
      c->qPlus[i] = c->aUni[i]*cosmology_mu(c,c->aUni[i])*(1 + ((2+6*x)/(x*sqrt(1+x))+3*log(x)-6*log(1+sqrt(1+x)))/(log(64)-9)); /* this is a^2*dDPlus/dt/H0 */

      c->aBox[i] = c->aUni[i]*cosmology_dc_factor(c,c->dPlus[i]);
      c->tCode[i] = 1.0 - tCodeFac*asinh(sqrt(aeq/c->aBox[i]));
    }
  
  /*
  //  Large a regime, solve ODEs
  */
  ASSERT(i > 0);

  tCodeFac = 0.5*sqrt(c->OmegaM);
  tPhysFac = tPhysUnit;

  y1[0] = c->tCode[i-1]/tCodeFac;
  y1[1] = c->tPhys[i-1]/tPhysFac;
  y1[2] = c->dPlus[i-1];
  y1[3] = c->qPlus[i-1];

  for(; i<n; i++)
    {
      a0 = c->aUni[i-1];
      da = c->aUni[i] - a0;

      /*  RK4 integration */
      for(j=0; j<4; j++) y0[j] = y1[j];
      cosmology_fill_table_integrate(c, a0,y1,f1);

      for(j=0; j<4; j++) y1[j] = y0[j] + 0.5*da*f1[j];
      cosmology_fill_table_integrate(c, a0+0.5*da,y1,f2);

      for(j=0; j<4; j++) y1[j] = y0[j] + 0.5*da*f2[j];
      cosmology_fill_table_integrate(c, a0+0.5*da,y1,f3);

      for(j=0; j<4; j++) y1[j] = y0[j] + da*f3[j];
      cosmology_fill_table_integrate(c, a0+da,y1,f4);

      for(j=0; j<4; j++) y1[j] = y0[j] + da*(f1[j]+2*f2[j]+2*f3[j]+f4[j])/6.0;

      c->tCode[i] = tCodeFac*y1[0];
      c->tPhys[i] = tPhysFac*y1[1];
      c->dPlus[i] = y1[2];
      c->qPlus[i] = y1[3];

      c->aBox[i] = c->aUni[i]*cosmology_dc_factor(c,c->dPlus[i]);
    }
} 


void cosmology_fill_table(CosmologyParameters *c, double amin, double amax)
{
  int i, imin, imax, iold;
  double dla = 1.0/c->ndex;
  double lamin, lamax;
  double *old_la = c->la;
  double *old_aUni = c->aUni;
  double *old_aBox = c->aBox;
  double *old_tCode = c->tCode;
  double *old_tPhys = c->tPhys;
  double *old_dPlus = c->dPlus;
  double *old_qPlus = c->qPlus;
  int old_size = c->size;

  if(amin > c->aLow) amin = c->aLow;
  lamin = dla*floor(c->ndex*log10(amin));
  lamax = dla*ceil(c->ndex*log10(amax));

  c->size = 1 + (int)(0.5+c->ndex*(lamax-lamin)); 
  ASSERT(fabs(lamax-lamin-dla*(c->size-1)) < 1.0e-14);

  c->la = NEWARR(c->size);     ASSERT(c->la != NULL);
  c->aUni = NEWARR(c->size);   ASSERT(c->aUni != NULL);
  c->aBox = NEWARR(c->size);   ASSERT(c->aBox != NULL);
  c->tCode = NEWARR(c->size);  ASSERT(c->tCode != NULL);
  c->tPhys = NEWARR(c->size);  ASSERT(c->tPhys != NULL);
  c->dPlus = NEWARR(c->size);  ASSERT(c->dPlus != NULL);
  c->qPlus = NEWARR(c->size);  ASSERT(c->qPlus != NULL);

  /*
  //  New log10(aUni) table
  */
  for(i=0; i<c->size; i++)
    {
      c->la[i] = lamin + dla*i;
    }

  if(old_size == 0)
    {
      /*
      //  Filling the table for the first time
      */
      cosmology_fill_table_piece(c,0,c->size);
    }
  else
    {
      /*
      //  Find if we need to expand the lower end
      */
      if(lamin < old_la[0])
	{
	  imin = (int)(0.5+c->ndex*(old_la[0]-lamin));
	  ASSERT(fabs(old_la[0]-lamin-dla*imin) < 1.0e-14);
	}
      else imin = 0;

      /*
      //  Find if we need to expand the upper end
      */
      if(lamax > old_la[old_size-1])
	{
	  imax = (int)(0.5+c->ndex*(old_la[old_size-1]-lamin));
	  ASSERT(fabs(old_la[old_size-1]-lamin-dla*imax) < 1.0e-14);
	}
      else imax = c->size - 1;
  
      /*
      //  Re-use the rest
      */
      if(lamin > old_la[0])
	{
	  iold = (int)(0.5+c->ndex*(lamin-old_la[0]));
	  ASSERT(fabs(lamin-old_la[0]-dla*iold) < 1.0e-14);
	}
      else iold = 0;

      memcpy(c->aUni+imin,old_aUni+iold,sizeof(double)*(imax-imin+1));
      memcpy(c->aBox+imin,old_aBox+iold,sizeof(double)*(imax-imin+1));
      memcpy(c->tCode+imin,old_tCode+iold,sizeof(double)*(imax-imin+1));
      memcpy(c->tPhys+imin,old_tPhys+iold,sizeof(double)*(imax-imin+1));
      memcpy(c->dPlus+imin,old_dPlus+iold,sizeof(double)*(imax-imin+1));
      memcpy(c->qPlus+imin,old_qPlus+iold,sizeof(double)*(imax-imin+1));

      DELETE(old_la);
      DELETE(old_aUni);
      DELETE(old_aBox);
      DELETE(old_tCode);
      DELETE(old_tPhys);
      DELETE(old_dPlus);
      DELETE(old_qPlus);

      /*
      //  Fill in additional pieces
      */
      if(imin > 0) cosmology_fill_table_piece(c,0,imin);
      if(imax < c->size-1) cosmology_fill_table_piece(c,imax,c->size);
    }
}


void cosmology_clear_table(CosmologyParameters *c)
{
  if(c->size > 0)
    {
      DELETE(c->la);
      DELETE(c->aUni);
      DELETE(c->aBox);
      DELETE(c->tCode);
      DELETE(c->tPhys);
      DELETE(c->dPlus);
      DELETE(c->qPlus);

      c->size = 0;
      c->la = NULL;
      c->aUni = NULL;
      c->aBox = NULL;
      c->tCode = NULL;
      c->tPhys = NULL;
      c->dPlus = NULL;
      c->qPlus = NULL;
    }
}


void cosmology_check_range(CosmologyParameters *c, double a)
{
  ASSERT((a > 1.0e-9) && (a < 1.0e9));

  if(c->size == 0) cosmology_init(c);

  if(a < c->aUni[0])
    {
      cosmology_fill_table(c,a,c->aUni[c->size-1]);
    }

  if(a > c->aUni[c->size-1])
    {
      cosmology_fill_table(c,c->aUni[0],a);
    }
}


void cosmology_set_thread_safe_range(CosmologyParameters *c, double amin, double amax)
{
  cosmology_check_range(c, amin);
  cosmology_check_range(c, amax);
}


double cosmology_get_value_from_table(CosmologyParameters *c, double a, double table[])
{
  // This is special case code for boundary conditions
  int idx;
  double la = log10(a);
  if (fabs(la - c->la[c->size-1]) < 1.0e-14) {
    return table[c->size-1];
  } else if (fabs(la - c->la[0]) < 1.0e-14) {
    return table[0];
  }

  idx = (int)(c->ndex*(la-c->la[0]));

  // Note that because we do idx+1 below, we need -1 here.
  ASSERT(idx>=0 && (idx<c->size-1));

  /*
  //  Do it as a function of aUni rather than la to ensure exact inversion
  */
  return table[idx] + (table[idx+1]-table[idx])/(c->aUni[idx+1]-c->aUni[idx])*(a-c->aUni[idx]);
}


int cosmology_find_index(CosmologyParameters *c, double v, double table[])
{
  int ic, il = 0;
  int ih = c->size - 1;

  if(v < table[0])
    {
      return -1;
    }
  if(v > table[c->size-1])
    {
      return c->size + 1;
    }

  while((ih-il) > 1)
    {
      ic = (il+ih)/2;
      if(v > table[ic]) /* special, not fully optimal form to avoid checking that il < c->size-1 */
	il = ic;
      else
	ih = ic;
    }

  ASSERT(il+1 < c->size);

  return il;
}


/*
//  Direct and inverse functions
*/
#define DEFINE_FUN(name,offset)			\
double name(CosmologyParameters *c, double a) \
{ \
  cosmology_check_range(c,a); \
  return cosmology_get_value_from_table(c,a,c->name) + offset; \
} \
double inv_##name(CosmologyParameters *c, double v) \
{ \
  int idx; \
  double *table; \
  if(c->size == 0) cosmology_init(c); \
  v -= offset; \
  table = c->name; \
  idx = cosmology_find_index(c,v,table); \
  while(idx < 0) \
    { \
      cosmology_check_range(c,0.5*c->aUni[0]); \
      table = c->name; \
      idx = cosmology_find_index(c,v,table); \
    } \
  while(idx > c->size) \
    { \
      cosmology_check_range(c,2.0*c->aUni[c->size-1]); \
      table = c->name; \
      idx = cosmology_find_index(c,v,table); \
    } \
  return c->aUni[idx] + (c->aUni[idx+1]-c->aUni[idx])/(table[idx+1]-table[idx])*(v-table[idx]); \
}

DEFINE_FUN(aBox,0.0);
DEFINE_FUN(tCode,c->tCodeOffset);
DEFINE_FUN(tPhys,0.0);
DEFINE_FUN(dPlus,0.0);
DEFINE_FUN(qPlus,0.0);

#undef DEFINE_FUN
