/* gsl_pow_int.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2004 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef __GSL_POW_INT_H__
#define __GSL_POW_INT_H__

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* Changes in GCC 4.2 and 4.3 relative to GCC 3.1 compiler:
*When compiling with -std=c99 or -std=gnu99, the extern inline keywords changes meaning. 
*GCC 4.3 conforms to the ISO C99 specification, where extern inline is very different 
*thing than the GNU extern inline extension. For the following code compiled with -std=c99, 
*extern inline int
*foo()
*{ return 5; }
*will result in a function definition for foo being emitted [produced] in the subsequent 
*object file, whereas previously there was none. As a result, files that use this extension 
*and compile in the C99 dialect will see many errors of the form: 
*       multiple definition of `foo'
*first defined here.
*To correct this error in GCC 4.2 and 4.3, we can add the following attribute:
*  __attribute__((__gnu_inline__)) 
* after extern inline to return to the old desired behavior.
*/


#ifdef HAVE_INLINE
/* add __attribute__((__gnu_inline__)) after "extern inline" to get old behaviour: */
extern inline __attribute__((__gnu_inline__)) double gsl_pow_2(const double x);
extern inline __attribute__((__gnu_inline__)) double gsl_pow_3(const double x);
extern inline __attribute__((__gnu_inline__)) double gsl_pow_4(const double x);
extern inline __attribute__((__gnu_inline__)) double gsl_pow_5(const double x);
extern inline __attribute__((__gnu_inline__)) double gsl_pow_6(const double x);
extern inline __attribute__((__gnu_inline__)) double gsl_pow_7(const double x);
extern inline __attribute__((__gnu_inline__)) double gsl_pow_8(const double x);
extern inline __attribute__((__gnu_inline__)) double gsl_pow_9(const double x);

extern inline __attribute__((__gnu_inline__)) double gsl_pow_2(const double x) { return x*x;   }
extern inline __attribute__((__gnu_inline__)) double gsl_pow_3(const double x) { return x*x*x; }
extern inline __attribute__((__gnu_inline__)) double gsl_pow_4(const double x) { double x2 = x*x;   return x2*x2;    }
extern inline __attribute__((__gnu_inline__)) double gsl_pow_5(const double x) { double x2 = x*x;   return x2*x2*x;  }
extern inline __attribute__((__gnu_inline__)) double gsl_pow_6(const double x) { double x2 = x*x;   return x2*x2*x2; }
extern inline __attribute__((__gnu_inline__)) double gsl_pow_7(const double x) { double x3 = x*x*x; return x3*x3*x;  }
extern inline __attribute__((__gnu_inline__)) double gsl_pow_8(const double x) { double x2 = x*x;   double x4 = x2*x2; return x4*x4; }
extern inline __attribute__((__gnu_inline__)) double gsl_pow_9(const double x) { double x3 = x*x*x; return x3*x3*x3; }
#else
double gsl_pow_2(const double x);
double gsl_pow_3(const double x);
double gsl_pow_4(const double x);
double gsl_pow_5(const double x);
double gsl_pow_6(const double x);
double gsl_pow_7(const double x);
double gsl_pow_8(const double x);
double gsl_pow_9(const double x);
#endif

double gsl_pow_int(double x, int n);

__END_DECLS

#endif /* __GSL_POW_INT_H__ */
