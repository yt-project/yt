/* -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2010 Krzysztof M. Gorski, Eric Hivon,
 *                          Benjamin D. Wandelt, Anthony J. Banday, 
 *                          Matthias Bartelmann, 
 *                          Reza Ansari & Kenneth M. Ganga 
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.jpl.nasa.gov
 *
 *----------------------------------------------------------------------------- */
/* ang2pix_nest.c */

/* Standard Includes */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Local Includes */
#include "healpix_vectors.h"

void ang2pix_nest( const long nside, double theta, double phi, long *ipix) {

  /* =======================================================================
   * subroutine ang2pix_nest(nside, theta, phi, ipix)
   * =======================================================================
   * gives the pixel number ipix (NESTED) corresponding to angles theta and phi
   *
   * the computation is made to the highest resolution available (nside=8192)
   * and then degraded to that required (by integer division)
   * this doesn't cost more, and it makes sure that the treatement of round-off 
   * will be consistent for every resolution
   * =======================================================================
   */
  
  double z, za, z0, tt, tp, tmp;
  int    face_num,jp,jm;
  long   ifp, ifm;
  int    ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;
  double piover2 = 0.5*M_PI, pi = M_PI, twopi = 2.0*M_PI;
  int    ns_max = 8192;
  static int x2pix[128], y2pix[128];
  static char setup_done = 0;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  if( theta<0. || theta>pi ) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }
  if( !setup_done ) {
    mk_xy2pix(x2pix,y2pix);
    setup_done = 1;
  }
  
  z  = cos(theta);
  za = fabs(z);
  z0 = 2./3.;
  if( phi>=twopi ) phi = phi - twopi;
  if( phi<0. )    phi = phi + twopi;
  tt = phi / piover2; /* in [0,4[ */
  
  if( za<=z0 ) { /* equatorial region */
    
    /* (the index of edge lines increase when the longitude=phi goes up) */
    jp = (int)floor(ns_max*(0.5 + tt - z*0.75)); /* ascending edge line index */
    jm = (int)floor(ns_max*(0.5 + tt + z*0.75)); /* descending edge line index */
    
    /* finds the face */
    ifp = jp / ns_max; /* in {0,4} */
    ifm = jm / ns_max;
    
    if( ifp==ifm ) face_num = (int)fmod(ifp,4) + 4; /* faces 4 to 7 */
    else if( ifp<ifm ) face_num = (int)fmod(ifp,4); /* (half-)faces 0 to 3 */
    else face_num = (int)fmod(ifm,4) + 8;           /* (half-)faces 8 to 11 */
    
    ix = (int)fmod(jm, ns_max);
    iy = ns_max - (int)fmod(jp, ns_max) - 1;
  }
  else { /* polar region, za > 2/3 */
    
    ntt = (int)floor(tt);
    if( ntt>=4 ) ntt = 3;
    tp = tt - ntt;
    tmp = sqrt( 3.*(1. - za) ); /* in ]0,1] */
    
    /* (the index of edge lines increase when distance from the closest pole
     * goes up)
     */
    /* line going toward the pole as phi increases */
    jp = (int)floor( ns_max * tp          * tmp ); 

    /* that one goes away of the closest pole */
    jm = (int)floor( ns_max * (1. - tp) * tmp );
    jp = (int)(jp < ns_max-1 ? jp : ns_max-1);
    jm = (int)(jm < ns_max-1 ? jm : ns_max-1);
    
    /* finds the face and pixel's (x,y) */
    if( z>=0 ) {
      face_num = ntt; /* in {0,3} */
      ix = ns_max - jm - 1;
      iy = ns_max - jp - 1;
    }
    else {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
    }
  }
  
  ix_low = (int)fmod(ix,128);
  ix_hi  =     ix/128;
  iy_low = (int)fmod(iy,128);
  iy_hi  =     iy/128;

  ipf = (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)+ (x2pix[ix_low]+y2pix[iy_low]);
  ipf = (long)(ipf / pow(ns_max/nside,2));     /* in {0, nside**2 - 1} */
  *ipix =(long)( ipf + face_num*pow(nside,2)); /* in {0, 12*nside**2 - 1} */
}
