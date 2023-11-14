/*							chdtrf.c
 *
 *	Chi-square distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * float df, x, y, chdtrf();
 *
 * y = chdtrf( df, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the left hand tail (from 0 to x)
 * of the Chi square probability density function with
 * v degrees of freedom.
 *
 *
 *                                  inf.
 *                                    -
 *                        1          | |  v/2-1  -t/2
 *  P( x | v )   =   -----------     |   t      e     dt
 *                    v/2  -       | |
 *                   2    | (v/2)   -
 *                                   x
 *
 * where x is the Chi-square variable.
 *
 * The incomplete gamma integral is used, according to the
 * formula
 *
 *	y = chdtr( v, x ) = igam( v/2.0, x/2.0 ).
 *
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       5000       3.2e-5      5.0e-6
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * chdtrf domain  x < 0 or v < 1        0.0
 */
/*							chdtrcf()
 *
 *	Complemented Chi-square distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * float v, x, y, chdtrcf();
 *
 * y = chdtrcf( v, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the right hand tail (from x to
 * infinity) of the Chi square probability density function
 * with v degrees of freedom:
 *
 *
 *                                  inf.
 *                                    -
 *                        1          | |  v/2-1  -t/2
 *  P( x | v )   =   -----------     |   t      e     dt
 *                    v/2  -       | |
 *                   2    | (v/2)   -
 *                                   x
 *
 * where x is the Chi-square variable.
 *
 * The incomplete gamma integral is used, according to the
 * formula
 *
 *	y = chdtr( v, x ) = igamc( v/2.0, x/2.0 ).
 *
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       5000       2.7e-5      3.2e-6
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * chdtrc domain  x < 0 or v < 1        0.0
 */
/*							chdtrif()
 *
 *	Inverse of complemented Chi-square distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * float df, x, y, chdtrif();
 *
 * x = chdtrif( df, y );
 *
 *
 *
 *
 * DESCRIPTION:
 *
 * Finds the Chi-square argument x such that the integral
 * from x to infinity of the Chi-square density is equal
 * to the given cumulative probability y.
 *
 * This is accomplished using the inverse gamma integral
 * function and the relation
 *
 *    x/2 = igami( df/2, y );
 *
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       10000      2.2e-5      8.5e-7
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * chdtri domain   y < 0 or y > 1        0.0
 *                     v < 1
 *
 */

/*								chdtr() */


/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>
#include <cmath>


namespace oyrke { namespace algorithm { namespace sse { namespace math {
    namespace {
        typedef float4v (*gamma_fn)(const float4v&, const float4v&);

        float4v chdtr_impl(const float4v& df, const float4v& x, const mask4v& x_invalid, 
                           gamma_fn func, float gamma_2nd_arg_scale) {
            float4v y = select(isnan(df), df);
            y = select(isnan(x), x);

            mask4v nanarg = isnan(df) || isnan(x);
            mask4v domainerr = x_invalid || df < 1.0f;
            // no special assign to y if domainerr, should return 0
            mask4v openarg = !(nanarg || domainerr);

            float4v yn = openarg.any_true()
                       ? (*func)(0.5f * df, gamma_2nd_arg_scale* x)
                       : float4v::zeros();
            y = select(openarg, yn, y);
            return y;
    }

    }


    float4v chdtrc(const float4v& df, const float4v& x) {
        float4v y = chdtr_impl(df, x, x<0.0f, igamc, 0.5f);
        return y;
    }

    
    
    float4v chdtr(const float4v& df, const float4v& x) {
        float4v y = chdtr_impl(df, x, x<0.0f, igam, 0.5f);
        return y;
    }



    float4v chdtri(const float4v& df, const float4v& x) {
        float4v y = chdtr_impl(df, x, x<0.0f || x>1.0f, igam, 1.0f);
        return y;
    }
}}}}
