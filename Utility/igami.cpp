/*							igamif()
 *
 *      Inverse of complemented imcomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * float a, x, y, igamif();
 *
 * x = igamif( a, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Given y, the function finds x such that
 *
 *  igamc( a, x ) = y.
 *
 * Starting with the approximate value
 *
 *         3
 *  x = a t
 *
 *  where
 *
 *  t = 1 - d - ndtri(y) sqrt(d)
 * 
 * and
 *
 *  d = 1/9a,
 *
 * the routine performs up to 10 Newton iterations to find the
 * root of igamc(a,x) - y = 0.
 *
 *
 * ACCURACY:
 *
 * Tested for a ranging from 0 to 100 and x from 0 to 1.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,100         5000       1.0e-5      1.5e-6
 *
 */

/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    float4v igami(const float4v& a, const float4v& x) {

        /* approximation to inverse function */
        float4v d = 1.0f/(9.0f*a);
        float4v t = (1.0f - d - ndtri(x) * sqrt(d));
        float4v x0 = a * t * t * t;

        float4v lgm = lgam(a);
        mask4v nanarga = isnan(a);
        mask4v nanargx = isnan(x);

        float4v y = select(nanarga, a);
        y = select(nanargx, x, y);

        mask4v  openargs = !(nanarga || nanargx);

        for (int i=0; i<10 && openargs.any_true(); ++i) {
            mask4v good = x0 > 0.0f;
            y = select(good, y, 0.0f);
            openargs = openargs && good;

            float4v r = igamc(a, x0);
            /* compute the derivative of the function at this point */
	        float4v s = (a - 1.0f) * log(x0) - x0 - lgm;
            good = s >= -MAXLOG;
            y = select(good, y, x0);
            openargs = openargs && good;

	        s = -exp(s);
            /* compute the step to the next approximation of x */
            good = s != 0.0f;
            y = select(good, y, x0);
            openargs = openargs && good;

            s = (r - x)/s;
	        x0 = x0 - s;
	
            good = (int4v(i) < 3) || abs(s/x0) >= 2.0f*MACHEPF;
            y = select(good, y, x0);
            openargs = openargs && good;
        }

        return y;
	}
}}}}
