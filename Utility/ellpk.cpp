/*							ellpkf.c
 *
 *	Complete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * float m1, y, ellpkf();
 *
 * y = ellpkf( m1 );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *            pi/2
 *             -
 *            | |
 *            |           dt
 * K(m)  =    |    ------------------
 *            |                   2
 *          | |    sqrt( 1 - m sin t )
 *           -
 *            0
 *
 * where m = 1 - m1, using the approximation
 *
 *     P(x)  -  log x Q(x).
 *
 * The argument m1 is used rather than m so that the logarithmic
 * singularity at m = 1 will be shifted to the origin; this
 * preserves maximum accuracy.
 *
 * K(0) = pi/2.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,1        30000       1.3e-7      3.4e-8
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * ellpkf domain      x<0, x>1           0.0
 *
 */

/*							ellpk.c */


/*
Cephes Math Library, Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

    const float P[] = {
        1.38629436111989062502E0,
        9.65735902811690126535E-2,
        3.08851465246711995998E-2,
        1.49380448916805252718E-2,
        8.79078273952743772254E-3,
        6.18901033637687613229E-3,
        6.87489687449949877925E-3,
        9.85821379021226008714E-3,
        7.97404013220415179367E-3,
        2.28025724005875567385E-3,
        1.37982864606273237150E-4
    };

    const float Q[] = {
        4.99999999999999999821E-1,
        1.24999999999870820058E-1,
        7.03124996963957469739E-2,
        4.88280347570998239232E-2,
        3.73774314173823228969E-2,
        3.01204715227604046988E-2,
        2.39089602715924892727E-2,
        1.54850516649762399335E-2,
        5.94058303753167793257E-3,
        9.14184723865917226571E-4,
        2.94078955048598507511E-5
    };

    const float C1 = 1.3862943611198906188E0; /* log(4) */

    float4v
    ellpk(const float4v& xin) {

        /* special cases
        < 0             0.0
        > 1             0.0
        0               inf
        */
        mask4v specialarg = xin <= 0.0f || xin > 1.0f;
        float4v yspecial = select(xin == 0.0f, float4v::infinites(), 0.0f);

        mask4v bigarg = xin > MACHEPF;
        mask4v smallarg = !(bigarg || specialarg);

        float4v xlog = log(xin);
        float4v ysmall = C1 - 0.5f*xlog;
        float4v ybig = eval_polynomial(xin, P) - xlog*eval_polynomial(xin, Q);
        float4v ynormal = select(bigarg, ybig, ysmall);

        float4v y = select(specialarg, yspecial, ynormal);
        y = select(isnan(xin), xin, y);

        return y;
    }
 }}}}
