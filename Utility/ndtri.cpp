/*							ndtrif.c
 *
 *	Inverse of Normal distribution function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, ndtrif();
 *
 * x = ndtrif( y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the argument, x, for which the area under the
 * Gaussian probability density function (integrated from
 * minus infinity to x) is equal to y.
 *
 *
 * For small arguments 0 < y < exp(-2), the program computes
 * z = sqrt( -2.0 * log(y) );  then the approximation is
 * x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
 * There are two rational functions P/Q, one for 0 < y < exp(-32)
 * and the other for y up to exp(-2).  For larger arguments,
 * w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain        # trials      peak         rms
 *    IEEE     1e-38, 1        30000       3.6e-7      5.0e-8
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition    value returned
 * ndtrif domain      x <= 0        -MAXNUM
 * ndtrif domain      x >= 1         MAXNUM
 *
 */


/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>
#include <cmath>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    /* sqrt(2pi) */
    const float SQRT_2PI = 2.50662827463100050242;

/* approximation for 0 <= |y - 0.5| <= 3/8 */
    const float P0[] = {
        -1.23916583867381258016E0,
         1.39312609387279679503E1,
        -5.66762857469070293439E1,
         9.80010754185999661536E1,
        -5.99633501014107895267E1,
    };
    const float Q0[] = {
        -1.18331621121330003142E0,
         1.59056225126211695515E1,
        -8.20372256168333339912E1,
         2.00260212380060660359E2,
        -2.25462687854119370527E2,
         8.63602421390890590575E1,
         4.67627912898881538453E0,
         1.95448858338141759834E0,
         1.00000000000000000000E0
    };

    /* Approximation for interval z = sqrt(-2 log y ) between 2 and 8
     * i.e., y between exp(-2) = .135 and exp(-32) = 1.27e-14.
     */
    const float P1[] = {
        -8.57456785154685413611E-4,
        -3.50424626827848203418E-2,
        -1.40256079171354495875E-1,
         2.18663306850790267539E0,
         1.46849561928858024014E1,
         4.40805073893200834700E1,
         5.71628192246421288162E1,
         3.15251094599893866154E1,
         4.05544892305962419923E0
    };

    const float Q1[] = {
        -9.33259480895457427372E-4,
        -3.80806407691578277194E-2,
        -1.42182922854787788574E-1,
         2.50464946208309415979E0,
         1.50425385692907503408E1,
         4.13172038254672030440E1,
         4.53907635128879210584E1,
         1.57799883256466749731E1,
         1.00000000000000000000E0,
    };


    /* Approximation for interval z = sqrt(-2 log y ) between 8 and 64
     * i.e., y between exp(-32) = 1.27e-14 and exp(-2048) = 3.67e-890.
     */

    const float P2[] = {
        6.23974539184983293730E-9,
        2.65806974686737550832E-6,
        3.01581553508235416007E-4,
        1.23716634817820021358E-2,
        2.01485389549179081538E-1,
        1.33303460815807542389E0,
        3.93881025292474443415E0,
        6.91522889068984211695E0,
        3.23774891776946035970E0,
    };

    const float Q2[] = {
        6.79019408009981274425E-9,
        2.89247864745380683936E-6,
        3.28014464682127739104E-4,
        1.34204006088543189037E-2,
        2.16236993594496635890E-1,
        1.37702099489081330271E0,
        3.67983563856160859403E0,
        6.02427039364742014255E0,
        1.00000000000000000000E0,
    };

    
    
    float4v ndtri(const float4v& y0) {

        float4v y = select(isnan(y0), y0);
        y = select(y0 <= 0.0f, -float4v::infinites(), y);
        y = select(y0 >= 1.0f,  float4v::infinites(), y);

        mask4v openarg = !(isnan(y0) || y0 <= 0.0f || y0 >= 1.0f);

        const float EXP_MINUS_TWO = 0.13533528323661269189; // exp(-2)
        mask4v negate = y0 > (1.0f - EXP_MINUS_TWO);
        float4v y0a = select(negate, 1.0f-y0, y0);

        mask4v largearg = y0a > EXP_MINUS_TWO && openarg;
        if (largearg.any_true()) {
	        float4v t = y0a - 0.5;
	        float4v t2 = t * t;
	        float4v p = t + t * (t2 * eval_polynomial(t2, P0)/eval_polynomial(t2, Q0));
	        p *= SQRT_2PI;
            y = select(largearg, p, y);
	    }
        openarg |= !largearg;

        if (openarg.any_true()) {
            float4v x = sqrt(-2.0f * log(y0a) );
            float4v x0 = x - log(x)/x;
            float4v z = 1.0f/x;

            float4v x1small, x1big;
            mask4v less8arg = x < 8.0f;
            if (less8arg.any_true()) {
                x1small = z * eval_polynomial(z, P1) / eval_polynomial(z, Q1);
            }
            if (less8arg.any_false()) {
                x1big = z = eval_polynomial(z, P2) / eval_polynomial(z, Q2);
            }
            float4v x1 = select(less8arg, x1small, x1big);
            x = x0 - x1;
            x = changesign(x, negate);

            y = select(openarg, x, y);
        }
        
        return y;
    }
}}}}
