/*							ndtrf.c
 *
 *	Normal distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, ndtrf();
 *
 * y = ndtrf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the area under the Gaussian probability density
 * function, integrated from minus infinity to x:
 *
 *                            x
 *                             -
 *                   1        | |          2
 *    ndtr(x)  = ---------    |    exp( - t /2 ) dt
 *               sqrt(2pi)  | |
 *                           -
 *                          -inf.
 *
 *             =  ( 1 + erf(z) ) / 2
 *             =  erfc(z) / 2
 *
 * where z = x/sqrt(2). Computation is via the functions
 * erf and erfc.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -13,0        50000       1.5e-5      2.6e-6
 *
 *
 * ERROR MESSAGES:
 *
 * See erfcf().
 *
 */
/*							erff.c
 *
 *	Error function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, erff();
 *
 * y = erff( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * The integral is
 *
 *                           x 
 *                            -
 *                 2         | |          2
 *   erf(x)  =  --------     |    exp( - t  ) dt.
 *              sqrt(pi)   | |
 *                          -
 *                           0
 *
 * The magnitude of x is limited to 9.231948545 for DEC
 * arithmetic; 1 or -1 is returned outside this range.
 *
 * For 0 <= |x| < 1, erf(x) = x * P(x**2); otherwise
 * erf(x) = 1 - erfc(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -9.3,9.3    50000       1.7e-7      2.8e-8
 *
 */
/*							erfcf.c
 *
 *	Complementary error function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, erfcf();
 *
 * y = erfcf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *  1 - erf(x) =
 *
 *                           inf. 
 *                             -
 *                  2         | |          2
 *   erfc(x)  =  --------     |    exp( - t  ) dt
 *               sqrt(pi)   | |
 *                           -
 *                            x
 *
 *
 * For small x, erfc(x) = 1 - erf(x); otherwise polynomial
 * approximations 1/x P(1/x**2) are computed.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -9.3,9.3    50000       3.9e-6      7.2e-7
 *
 *
 * ERROR MESSAGES:
 *
 *   message           condition              value returned
 * erfcf underflow    x**2 > MAXLOGF              0.0
 *
 *
 */


/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {


    /* erfc(x) = exp(-x^2) P(1/x), 1 < x < 2 */
    const float P[] = {
         5.638259427386472E-001,
        -2.741127028184656E-001,
         3.404879937665872E-001,
        -4.944515323274145E-001,
         6.210004621745983E-001,
        -5.824733027278666E-001,
         3.687424674597105E-001,
        -1.387039388740657E-001,
         2.326819970068386E-002
    };

    /* erfc(x) = exp(-x^2) 1/x P(1/x^2), 2 < x < 14 */
    const float R[] = {
         5.641895067754075E-001,
        -2.820767439740514E-001,
         4.218463358204948E-001,
        -1.015265279202700E+000,
         2.921019019210786E+000,
        -7.495518717768503E+000,
         1.297719955372516E+001,
        -1.047766399936249E+001
    };

    /* erf(x) = x P(x^2), 0 < x < 1 */
    const float T[] = {
     1.128379165726710E+000,
    -3.761262582423300E-001,
     1.128358514861418E-001,
    -2.685381193529856E-002,
     5.188327685732524E-003,
    -8.010193625184903E-004,
     7.853861353153693E-005
    };


    namespace {
        __forceinline float4v 
        erfc_impl(const float4v& x, const mask4v& unprocessed) {
            float4v z = x * x;
            mask4v underflow = z > MAXLOG;
            float4v yunderflow = select(x < 0.0f, float4v(2.0f)); // else 0.0
            mask4v openargs = unprocessed && !underflow;
            float4v y = select(underflow, yunderflow, y);

            // anything left to process?
            if (openargs.any_true()) {
                z = exp(-z);
                float4v q = 1.0f/x;
                float4v q2 = q * q;

                // add branch???
                float4v p = select(x < 2.0f, eval_polynomial(q2, P), eval_polynomial(q2, R));

                float4v t = z * q * p;
                t = select(x < 0.0f, 2.0f-t, t);
                t = select(t == 0.0f, yunderflow, t);

                y = select(openargs, t, y);
            }

            return y;
        }



        __forceinline float4v 
        erf_impl(const float4v& x) {
            float4v x2 = x * x;
            float4v y = x * eval_polynomial(x, T);
            return y;
        }
    }


    float4v erfc(const float4v& xin) {
        float4v x = abs(xin);

        mask4v nanarg = isnan(xin);

        float4v y = select(nanarg, xin);
        mask4v less1arg = x < 1.0f;
        mask4v openarg = !(isnan(x) || less1arg);

        if (less1arg.any_true()) {
            float4v yless1 = 1.0f - erf_impl(x);
            y = select(less1arg, yless1);
        }

        float4v yn = erfc_impl(x, openarg);
        y = select(openarg, yn, y);
   
        return y;
    }




    float4v erf(const float4v& xin) {
        float4v xa      = abs(xin);
        mask4v bigarg   = xa > 1.0f;
        float4v y = select(isnan(xa), xin);

        if (bigarg.any_true()) {
            float4v ybig = 1.0f - erfc(xin); // TODO manip arg to avoid special path
            y = select(bigarg, ybig, y);
        }

        mask4v smallarg = xa <= 1.0f;
        if (smallarg.any_true()) {
            float4v ysmall = erf_impl(xin);
            y = blend(y, select(smallarg, ysmall));
        }

        return y;
    }
}}}}
