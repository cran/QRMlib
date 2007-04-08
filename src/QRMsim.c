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
/* Modified by Scott Ulman, 05/23/2006 to work in R rather than S-Plus */
/*  #include <S.h>     */
#include <R.h>
#include <R_ext/Random.h>
#include "QRMsim.h"

void frank(long *n, double *theta, double *output)
{
 /* Modified by Scott Ulman, 05/23/2006 to work in R rather than S-Plus */
  /*  S_EVALUATOR         */

  int i, k;
  double U, alpha, thetaval, p;
  thetaval = *theta;
  alpha = 1-exp(-thetaval);
  RANDIN; 

  for( i=0; i<*n; i++ ) {
    k = 1;
    /* Modified by Scott Ulman, 05/23/2006 to work in R rather than S-Plus */
    /*  U = UNIF;       */
    U = unif_rand();
    p = alpha/thetaval;
    while(U > p) {
      k++;
      p = p + pow(alpha,k)/(k*thetaval);
    }
    *(output+i)=k;
  }

  RANDOUT;
}


void rgig(long *n, double *r, double *s, double *p, double *k1, double *k2, double *lambda, double *chi, double *psi, double *s1, double *s2, double *xsim)
{
 /* Modified by Scott Ulman, 05/23/2006 to work in R rather than S-Plus */
  /* S_EVALUATOR     */

  int i=0;
  long count =0;
  double U, Ustar, level, x;
  RANDIN;

  while (i < *n){
   /* Modified by Scott Ulman, 05/23/2006 to work in R rather than S-Plus */
  /* U = UNIF;
    Ustar = UNIF;    */
    U = unif_rand();
    Ustar = unif_rand();

    count++;
    if (U <= *r){
      x = log(1+(*s)*U/(*k1))/(*s);  
      level = log(ef(x, *lambda, *chi,*psi+2*(*s))/(*s1));
      if (log(Ustar) <= level){
	*(xsim+i) = x;
	i++;
      }
    }
    else {
      x = -log((*p)*(1-U)/(*k2))/(*p);
      level = log(ef(x, *lambda, *chi,*psi-2*(*p))/(*s2));
      if (log(Ustar) <= level){
	*(xsim+i) = x;
	i++;
      }
    }

  }
  *n = count;
  RANDOUT;
}

double ef(double x, double lambda, double chi, double psi)
{
  double result;

  result = pow(x,lambda-1.0)*exp(-0.5*(psi*x + chi/x));
  return result;
}


