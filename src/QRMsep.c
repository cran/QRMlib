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

#include <math.h>
#include "QRMsep.h"


void SEprocExciteFunc(int*n, double *times, int *nmarks, double *marktimes, double *marks, double *beta, int *model, double *result)
{
  int i=0, j;
  double thetime, gamma, delta, rho=0.0, tmp;

  gamma = *beta;
  delta = 0.0;
  if (*model ==2 )
    /* Hawkes with mark influence */
    delta = *(beta+1);
  if (*model == 3)
    /* ETAS without mark influence */
    rho = *(beta+1);
  if (*model == 4){
    /* ETAS with mark influence */
    rho = *(beta+1);
    delta = *(beta+2);
  }
  while (i < *n){
    tmp = 0.0;
    thetime = *(times+i);
    j = 0;
    while ((*(marktimes+j) < thetime) & (j < *nmarks)){
      if (*model == 1)
	tmp = tmp + contribH((thetime-*(marktimes+j)),0.0,gamma,delta);
      if (*model == 2)
	tmp = tmp + contribH((thetime-*(marktimes+j)),*(marks+j),gamma,delta);
      if (*model == 3)
	tmp = tmp + contribE((thetime-*(marktimes+j)),0.0,gamma,rho,delta);
      if (*model == 4)
	tmp = tmp + contribE((thetime-*(marktimes+j)),*(marks+j),gamma,rho,delta);
      j++;
    }
    *(result+i) = tmp;
    i++;
  }
}


double contribH(double s, double y, double gammaval, double deltaval)
{
  double result;
  result = (1+deltaval*y)*exp(-gammaval*s);
  return result;
}

double contribE(double s, double y, double gammaval, double rhoval, double deltaval)
{
  double result;
  result = (1+deltaval*y)/pow(1+s/gammaval,(1.0+rhoval));
  return result;
}

