/*							k0f.c
 *
 *	Modified Bessel function, third kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, k0f();
 *
 * y = k0f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of the third kind
 * of order zero of the argument.
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at 2000 random points between 0 and 8.  Peak absolute
 * error (relative when K0 > 1) was 1.46e-14; rms, 4.26e-15.
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       7.8e-7      8.5e-8
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 *  K0 domain          x <= 0          MAXNUM
 *
 */
/*							k0ef()
 *
 *	Modified Bessel function, third kind, order zero,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, k0ef();
 *
 * y = k0ef( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of the third kind of order zero of the argument.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       8.1e-7      7.8e-8
 * See k0().
 *
 */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    /* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
     * in the interval [0,2].  The odd order coefficients are all
     * zero; only the even order coefficients are listed.
     * 
     * lim(x->0){ K0(x) + log(x/2) I0(x) } = -EUL.
     */

    const float A[] = {
        -5.35327393233902768720E-1f,
         3.44289899924628486886E-1f,
         3.59799365153615016266E-2f,
         1.26461541144692592338E-3f,
         2.28621210311945178607E-5f,
         2.53479107902614945675E-7f,
         1.90451637722020886025E-9f
    };



    /* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
     * in the inverted interval [2,infinity].
     * 
     * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
     */

    const float B[] = {
         2.44030308206595545468E0f,
        -3.14481013119645005427E-2f,
         1.56988388573005337491E-3f,
        -1.28495495816278026384E-4f,
         1.39498137188764993662E-5f,
        -1.83175552271911948767E-6f,
         2.76681363944501510342E-7f,
        -4.66048989768794782956E-8f,
         8.57403401741422608519E-9f,
        -1.69753450938905987466E-9f,
    };

    
    namespace {
        float4v k0_noexp(const float4v& x, mask4v& smallarg, mask4v& bigarg) {
            float4v y = select(isnan(x), x);
            y = select(x <= 0.0f, float4v::infinites(), y);

            mask4v openarg = !(isnan(x) || x <= 0.0f);
            mask4v less2arg = x <= 2.0f;
            
            smallarg = less2arg && openarg;
            float4v ys, yb;
         
            if (smallarg.any_true()) {
                float4v t = x*x - 2.0f;
                ys = eval_chebychev(t, A)  -  log(0.5f*x) * i0(x);

            }

            bigarg = !less2arg && openarg;
            if (bigarg.any_true()) {
                float4v t = 8.0f/x - 2.0f;
                yb = eval_chebychev(t, B) / sqrt(x);
            }

            y = select(less2arg, blend(yb, ys), y);
            return y;
        }
    }
    
    /*							k0.c	*/
    float4v k0(const float4v& x) {
        mask4v small, big;

        float4v t = k0_noexp(x, small, big);
        float4v t_exp = big.any_true() ? exp(-x) : float4v::zeros();

        float4v y = select(big, t*t_exp, t);
        return y;
    }


    
    float4v k0e(const float4v& x) {
        mask4v small, big;

        float4v t = k0_noexp(x, small, big);
        float4v t_exp = small.any_true() ? exp(x) : float4v::zeros();

        float4v y = select(small, t*t_exp, t);
        return y;
    }
}}}}
