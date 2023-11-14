/*							gdtrf.c
 *
 *	Gamma distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * float a, b, x, y, gdtrf();
 *
 * y = gdtrf( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from zero to x of the gamma probability
 * density function:
 *
 *
 *                x
 *        b       -
 *       a       | |   b-1  -at
 * y =  -----    |    t    e    dt
 *       -     | |
 *      | (b)   -
 *               0
 *
 *  The incomplete gamma integral is used, according to the
 * relation
 *
 * y = igam( b, ax ).
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       5000       5.8e-5      3.0e-6
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * gdtrf domain        x < 0            0.0
 *
 */
/*							gdtrcf.c
 *
 *	Complemented gamma distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * float a, b, x, y, gdtrcf();
 *
 * y = gdtrcf( a, b, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the integral from x to infinity of the gamma
 * probability density function:
 *
 *
 *               inf.
 *        b       -
 *       a       | |   b-1  -at
 * y =  -----    |    t    e    dt
 *       -     | |
 *      | (b)   -
 *               x
 *
 *  The incomplete gamma integral is used, according to the
 * relation
 *
 * y = igamc( b, ax ).
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       5000       9.1e-5      1.5e-5
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * gdtrcf domain        x < 0            0.0
 *
 */

/*							gdtr()  */


/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    namespace {
        typedef float4v (*gamma_fn)(const float4v&, const float4v&);

        float4v gdtr_impl(const float4v& a, const float4v& b, const float4v& x, gamma_fn func) {
            float4v t = a+b+x;
            mask4v nanarg = isnan(t);
            float4v y = select(nanarg, t);

            mask4v domainerr = x < 0.0f;
            mask4v openarg = !(nanarg || domainerr);

            if (openarg.any_true()) {
                float4v x_mod = select(openarg, x, float4v::NaNs());
                float4v yn = (*func)(b, a*x_mod);
                y = blend(y, select(openarg, yn));
            }

            return y;
        }
    }


    
    
    float4v gdtr(const float4v& a, const float4v& b, const float4v& x) {
        return gdtr_impl(a, b, x, igam);
    }

        
    
    float4v gdtrc(const float4v& a, const float4v& b, const float4v& x) {
        return gdtr_impl(a, b, x, igamc);
    }
 }}}}