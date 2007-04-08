/* QRMlib: version 1.4 */
/* this file is a component of QRMlib */

/* Copyright (C) 2005-06 Alexander McNeil */
/* This file has been modified  by Scott Ulman, FDIC to work in R rather than SPlus */

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

/* definitions for random numbers in C under SPlus */
/* Modified by Scott Ulman, 05/23/2006 to work in R rather than S-Plus */
/* All the following must be commented out: 
# define RANDIN seed_in((long *)NULL, S_evaluator);
# define RANDOUT seed_out((long *)NULL, S_evaluator);
# define UNIF unif_rand(S_evaluator);
*/

/* Modified by Scott Ulman, 05/23/2006 to work in R rather than S-Plus */
# define RANDIN GetRNGstate();
# define RANDOUT PutRNGstate();
/* # define UNIF unif_rand();     */


/* Mixing distribution for Frank copula */

void frank(long *n, double *theta, double *output);

/* Simulating generalised inverse Gaussian */
void rgig(long *n, double *r, double *s, double *p, double *k1, double *k2, double *lambda, double *chi, double *psi, double *s1, double *s2, double *xsim);



double ef(double x, double lambda, double chi, double psi);



