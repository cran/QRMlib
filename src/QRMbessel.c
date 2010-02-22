/* QRMlib: version 1.4 */
/* this file is a component of QRMlib */
/* Copyright (C) 2005-06 Alexander McNeil */
/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA. */

/* Contact: Alexander J. McNeil:  mcneil@math.ethz.ch */
/* SU: Move these files into subfolder of our src so the user doesn't need gsl installed */
/*
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_psi.h>
*/

#include "gsl/specfunc/gsl_sf_bessel.h"
#include "gsl/specfunc/gsl_sf_psi.h"

#include "QRMbessel.h"

void besselM3(int *lambda, double *x, int *n, double *err, double *y, int *logvalue)
{
  gsl_sf_result result;
  int i, nint;
  double nu;

  if (*logvalue == 0){
    nint = (int)*lambda;
    for (i = 0; i < *n; i++){
      gsl_sf_bessel_Kn_e(nint, *(x+i), &result);
      *(y+i) = result.val;
      *(err+i) = result.err;
    }
  }
  if (*logvalue ==1){
    nu = (double)*lambda;
    for (i = 0; i < *n; i++){
      gsl_sf_bessel_lnKnu_e(nu, *(x+i), &result);
      *(y+i) = result.val;
      *(err+i) = result.err;
    }
  }

}

void besselM3f(double *lambda, double *x, int *n, double *err, double *y, int *logvalue)
{
  gsl_sf_result result;
  int i;

  if (*logvalue == 0){
    for (i = 0; i < *n; i++){
      gsl_sf_bessel_Knu_e(*lambda, *(x+i), &result);
      *(y+i) = result.val;
      *(err+i) = result.err;
    }
  }
  if (*logvalue == 1){
    for (i = 0; i < *n; i++){
      gsl_sf_bessel_lnKnu_e(*lambda, *(x+i), &result);
      *(y+i) = result.val;
      *(err+i) = result.err;
    }
  }
 }


void besselM3z(double *x, int *n, double *err, double *y)
{
  gsl_sf_result result;
  int i;

  for (i = 0; i < *n; i++){
    gsl_sf_bessel_K0_e(*(x+i), &result);
    *(y+i) = result.val;
    *(err+i) = result.err;
  }

}


void psifunc(double *x, int *n, double *err, double *y)
{
  gsl_sf_result result;
  int i;

  for (i = 0; i < *n; i++){
    gsl_sf_psi_e(*(x+i), &result);
    *(y+i) = result.val;
    *(err+i) = result.err;
  }
}

