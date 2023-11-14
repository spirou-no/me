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
#include <cmath>

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

    float4v ndtr(const float4v& xin) {
        float4v x = xin * SQRT2_HALF;
        float4v z = abs(x);

        mask4v smallarg = z < SQRT2_HALF;
        mask4v bigarg   = z >= SQRT2_HALF;
        float4v y;

        if (smallarg.any_true()) {
            y = 1.0f +  erf(x);  // should be z?
            y = select(smallarg, y);  
        }

        if (bigarg.any_true()) {
            float4v ybig = erfc(z);
            ybig = select(x > 0.0f, 2.0f-ybig, ybig);
            y = blend(y, ybig);  // bigarg and smallarg mutually exclusive, ok to blend
        }

        y *= 0.5f;
        y = select(isnan(xin), xin, y);

        return y;
    }
}}}}
