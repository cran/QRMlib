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
/* Bessel function of 3rd kind */
void besselM3( int *lambda, double *x, int *n, double *err, double *y, int *logvalue);

void besselM3f( double *lambda, double *x, int *n, double *err, double *y, int *logvalue);

void besselM3z( double *x, int *n, double *err, double *y);

/* psi or digamma function */
void psifunc(double *x, int *n, double *err, double *y);

