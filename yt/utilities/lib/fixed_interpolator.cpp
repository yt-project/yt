/*******************************************************************************
*******************************************************************************/

//
// A small, tiny, itty bitty module for computation-intensive interpolation
// that I can't seem to make fast in Cython
//

#include "fixed_interpolator.hpp"

#define VINDEX(A,B,C) data[((((A)+ci[0])*(ds[1]+1)+((B)+ci[1]))*(ds[2]+1)+ci[2]+(C))]
//  (((C*ds[1])+B)*ds[0]+A)
#define OINDEX(A,B,C) data[(A)*(ds[1]+1)*(ds[2]+1)+(B)*ds[2]+(B)+(C)]

npy_float64 fast_interpolate(int ds[3], int ci[3], npy_float64 dp[3],
                             npy_float64 *data)
{
    int i;
    npy_float64 dv, dm[3];
    for(i=0;i<3;i++)dm[i] = (1.0 - dp[i]);
    dv  = 0.0;
    dv += VINDEX(0,0,0) * (dm[0]*dm[1]*dm[2]);
    dv += VINDEX(0,0,1) * (dm[0]*dm[1]*dp[2]);
    dv += VINDEX(0,1,0) * (dm[0]*dp[1]*dm[2]);
    dv += VINDEX(0,1,1) * (dm[0]*dp[1]*dp[2]);
    dv += VINDEX(1,0,0) * (dp[0]*dm[1]*dm[2]);
    dv += VINDEX(1,0,1) * (dp[0]*dm[1]*dp[2]);
    dv += VINDEX(1,1,0) * (dp[0]*dp[1]*dm[2]);
    dv += VINDEX(1,1,1) * (dp[0]*dp[1]*dp[2]);
    /*assert(dv < -20);*/
    return dv;
}

npy_float64 offset_interpolate(int ds[3], npy_float64 dp[3], npy_float64 *data)
{
    int i;
    npy_float64 dv, vz[4];

    dv = 1.0 - dp[2];
    vz[0] = dv*OINDEX(0,0,0) + dp[2]*OINDEX(0,0,1);
    vz[1] = dv*OINDEX(0,1,0) + dp[2]*OINDEX(0,1,1);
    vz[2] = dv*OINDEX(1,0,0) + dp[2]*OINDEX(1,0,1);
    vz[3] = dv*OINDEX(1,1,0) + dp[2]*OINDEX(1,1,1);

    dv = 1.0 - dp[1];
    vz[0] = dv*vz[0] + dp[1]*vz[1];
    vz[1] = dv*vz[2] + dp[1]*vz[3];

    dv = 1.0 - dp[0];
    vz[0] = dv*vz[0] + dp[0]*vz[1];

    return vz[0];
}

void offset_fill(int ds[3], npy_float64 *data, npy_float64 gridval[8])
{
    gridval[0] = OINDEX(0,0,0);
    gridval[1] = OINDEX(1,0,0);
    gridval[2] = OINDEX(1,1,0);
    gridval[3] = OINDEX(0,1,0);
    gridval[4] = OINDEX(0,0,1);
    gridval[5] = OINDEX(1,0,1);
    gridval[6] = OINDEX(1,1,1);
    gridval[7] = OINDEX(0,1,1);
}

void vertex_interp(npy_float64 v1, npy_float64 v2, npy_float64 isovalue,
                   npy_float64 vl[3], npy_float64 dds[3],
                   npy_float64 x, npy_float64 y, npy_float64 z,
                   int vind1, int vind2)
{
    /*if (fabs(isovalue - v1) < 0.000001) return 0.0;
    if (fabs(isovalue - v2) < 0.000001) return 1.0;
    if (fabs(v1 - v2) < 0.000001) return 0.0;*/
    int i;
    static npy_float64 cverts[8][3] =
        {{0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
         {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1}};

    npy_float64 mu = ((isovalue - v1) / (v2 - v1));

    if (fabs(1.0 - isovalue/v1) < 0.000001) mu = 0.0;
    if (fabs(1.0 - isovalue/v2) < 0.000001) mu = 1.0;
    if (fabs(v1/v2) < 0.000001) mu = 0.0;

    vl[0] = x; vl[1] = y; vl[2] = z;
    for (i=0;i<3;i++)
        vl[i] += dds[i] * cverts[vind1][i]
               + dds[i] * mu*(cverts[vind2][i] - cverts[vind1][i]);
}

npy_float64 trilinear_interpolate(int ds[3], int ci[3], npy_float64 dp[3],
				  npy_float64 *data)
{
    /* dims is one less than the dimensions of the array */
    int i;
    npy_float64 dm[3], vz[4];
  //dp is the distance to the plane.  dm is val, dp = 1-val
    for(i=0;i<3;i++)dm[i] = (1.0 - dp[i]);

  //First interpolate in z
    vz[0] = dm[2]*VINDEX(0,0,0) + dp[2]*VINDEX(0,0,1);
    vz[1] = dm[2]*VINDEX(0,1,0) + dp[2]*VINDEX(0,1,1);
    vz[2] = dm[2]*VINDEX(1,0,0) + dp[2]*VINDEX(1,0,1);
    vz[3] = dm[2]*VINDEX(1,1,0) + dp[2]*VINDEX(1,1,1);

  //Then in y
    vz[0] = dm[1]*vz[0] + dp[1]*vz[1];
    vz[1] = dm[1]*vz[2] + dp[1]*vz[3];

  //Then in x
    vz[0] = dm[0]*vz[0] + dp[0]*vz[1];
    /*assert(dv < -20);*/
    return vz[0];
}

void eval_gradient(int ds[3], npy_float64 dp[3],
				  npy_float64 *data, npy_float64 *grad)
{
    // We just take some small value

    int i;
    npy_float64 denom, plus, minus, backup, normval;

    normval = 0.0;
    for (i = 0; i < 3; i++) {
      backup = dp[i];
      grad[i] = 0.0;
      if (dp[i] >= 0.95) {plus = dp[i]; minus = dp[i] - 0.05;}
      else if (dp[i] <= 0.05) {plus = dp[i] + 0.05; minus = 0.0;}
      else {plus = dp[i] + 0.05; minus = dp[i] - 0.05;}
      //fprintf(stderr, "DIM: %d %0.3lf %0.3lf\n", i, plus, minus);
      denom = plus - minus;
      dp[i] = plus;
      grad[i] += offset_interpolate(ds, dp, data) / denom;
      dp[i] = minus;
      grad[i] -= offset_interpolate(ds, dp, data) / denom;
      dp[i] = backup;
      normval += grad[i]*grad[i];
    }
    if (normval != 0.0){
      normval = sqrt(normval);
      for (i = 0; i < 3; i++) grad[i] /= -normval;
      //fprintf(stderr, "Normval: %0.3lf %0.3lf %0.3lf %0.3lf\n",
      //        normval, grad[0], grad[1], grad[2]);
    }else{
      grad[0]=grad[1]=grad[2]=0.0;
    }
}
