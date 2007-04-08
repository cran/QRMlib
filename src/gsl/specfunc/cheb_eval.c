/* Modified by Scott Ulman, 05/24/2006 to work in R rather than S-Plus */
/* Added four header files */
/* SU: Move these files into subfolder of our src so the user doesn't need gsl installed */
/* 
#include <gsl/gsl_errno.h>   //needed for GSL_SUCCESS enum 
#include <gsl/gsl_machine.h> //needed for GSL_DBL_EPSILON 
#include <gsl/gsl_sf_result.h> //needed for gsl_sf_result 
#include <gsl/chebyshev.h>     //needed to define cheb_series  
*/

#include "../gsl_errno.h"  //needed for GSL_SUCCESS enum 
#include "../gsl_machine.h" /*needed for GSL_DBL_EPSILON */
#include "gsl_sf_result.h"  /*needed for gsl_sf_result */
/* Moved chebyshev.h to general header file folder along wtih cheb_eval.c since */
/* I now call cheb_eval.c via #include <gsl/cheb_eval.c> */
#include "chebyshev.h"     /*needed to define cheb_series   */


#include <R.h>     /*needed for R_INLINE */

/* Modified by Scott Ulman, 05/24/2006 to work in R rather than S-Plus */
/* Tested to see whether in-lining supported by R platform and compiler using method  */
/* on p. 79 (section 5.14) of R-exts.pdf.  Replaced next line with the one below*/
/*static inline int  */
static R_INLINE int
cheb_eval_e(const cheb_series * cs,
            const double x,
            gsl_sf_result * result)
{
  int j;
  double d  = 0.0;
  double dd = 0.0;

  double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
  double y2 = 2.0 * y;

  double e = 0.0;

  for(j = cs->order; j>=1; j--) {
    double temp = d;
    d = y2*d - dd + cs->c[j];
    e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
    dd = temp;
  }

  { 
    double temp = d;
    d = y*d - dd + 0.5 * cs->c[0];
    e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
  }

  result->val = d;
  result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

  return GSL_SUCCESS;
}

