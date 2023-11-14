/*							j0f.c
 *
 *	Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, j0f();
 *
 * y = j0f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 2] and
 * (2, infinity). In the first interval the following polynomial
 * approximation is used:
 *
 *
 *        2         2         2
 * (w - r  ) (w - r  ) (w - r  ) P(w)
 *       1         2         3   
 *
 *            2
 * where w = x  and the three r's are zeros of the function.
 *
 * In the second interval, the modulus and phase are approximated
 * by polynomials of the form Modulus(x) = sqrt(1/x) Q(1/x)
 * and Phase(x) = x + 1/x R(1/x^2) - pi/4.  The function is
 *
 *   j0(x) = Modulus(x) cos( Phase(x) ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 2        100000      1.3e-7      3.6e-8
 *    IEEE      2, 32       100000      1.9e-7      5.4e-8
 *
 */
/*							y0f.c
 *
 *	Bessel function of the second kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, y0f();
 *
 * y = y0f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind, of order
 * zero, of the argument.
 *
 * The domain is divided into the intervals [0, 2] and
 * (2, infinity). In the first interval a rational approximation
 * R(x) is employed to compute
 *
 *                  2         2         2
 * y0(x)  =  (w - r  ) (w - r  ) (w - r  ) R(x)  +  2/pi ln(x) j0(x).
 *                 1         2         3   
 *
 * Thus a call to j0() is required.  The three zeros are removed
 * from R(x) to improve its numerical stability.
 *
 * In the second interval, the modulus and phase are approximated
 * by polynomials of the form Modulus(x) = sqrt(1/x) Q(1/x)
 * and Phase(x) = x + 1/x S(1/x^2) - pi/4.  Then the function is
 *
 *   y0(x) = Modulus(x) sin( Phase(x) ).
 *
 *
 *
 *
 * ACCURACY:
 *
 *  Absolute error, when y0(x) < 1; else relative error:
 *
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,  2       100000      2.4e-7      3.4e-8
 *    IEEE      2, 32       100000      1.8e-7      5.3e-8
 *
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/math/constants.h>
#include <cmath>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    const float MO[] = {
        7.978845717621440E-001f,
       -3.355424622293709E-006f,
       -4.969382655296620E-002f,
       -3.560281861530129E-003f,
        1.197549369473540E-001f,
       -2.145007480346739E-001f,
        1.864949361379502E-001f,
       -6.838999669318810E-002f
    };


    const float PH[] = {
       -1.249992184872738E-001f,
        6.490598792654666E-002f,
       -1.939906941791308E-001f,
        1.001973420681837E+000f,
       -4.974978466280903E+000f,
        1.756221482109099E+001f,
       -3.630592630518434E+001f,
        3.242077816988247E+001f
    };


    const float YP[] = {
        1.707584643733568E-001f,
       -1.584289289821316E-002f,
        5.344486707214273E-004f,
       -9.413212653797057E-006f,
        9.454583683980369E-008f
    };

    const float YZ1 =  0.43221455686510834878f;
    const float YZ2 = 22.401876406482861405f;
    const float YZ3 = 64.130620282338755553f;

    const float DR1 =  5.78318596294678452118f;
    /*
    static float DR2 = 30.4712623436620863991;
    static float DR3 = 74.887006790695183444889;
    */

    const float JP[] = {
        -1.729150680240724E-001f,
         1.332913422519003E-002f,
        -3.969646342510940E-004f,
         6.388945720783375E-006f,
        -6.068350350393235E-008f
    };


    
    float4v j0(const float4v& xin) {

        float4v x = abs(x);
        mask4v smallarg = x < 2.0f;
        mask4v bigarg   = x > 2.0f;
        float4v ybig, ysmall;

        if (smallarg.any_true()) {
            float4v z = x*x;
            ysmall = select(x < 1.0e-3f, 1.0f - 0.25f*z, (z-DR1)*eval_polynomial(z, JP));
        }

        if (bigarg.any_true()) {
            float4v q = 1.0f/x;
            float4v w = sqrt(q);
            float4v p = w * eval_polynomial(q, MO);

            w = q*q;
            float4v xn = q * eval_polynomial(w, PH) - PI_QUART;
            ybig = p * cos(xn + x);
        }

        float4v y = blend(ybig, ysmall);
        y = select(isnan(x), xin, y);

        return y;
	}


/*							y0() 2	*/
/* Bessel function of second kind, order zero	*/

/* Rational approximation coefficients YP[] are used for x < 6.5.
 * The function computed is  y0(x)  -  2 ln(x) j0(x) / pi,
 * whose value at x = 0 is  2 * ( log(0.5) + EUL ) / pi
 * = 0.073804295108687225 , EUL is Euler's constant.
 */

    const float TWO_INV_PI =  0.636619772367581343075535f; /* 2/pi */

    float4v y0(const float4v& x) {
        
        float4v ysmall, ybig;
        mask4v smallarg = x <= 2.0f;
        mask4v bigarg   = x > 2.0f;

        if (smallarg.any_true()) {
            float4v z = x*x;
            float4v w = (z-YZ1) * eval_polynomial(z, YP);
            w += TWO_INV_PI * log(x) * j0(x);
            ysmall = select(x <= 0.0f, float4v::NaNs(), w);
        }

        if (bigarg.any_true()) {
            float4v q = 1.0f/x;
            float4v w = sqrt(q);
            float4v p = w * eval_polynomial(q, MO);

            w = q*q;
            float4v xn = q * eval_polynomial(w, PH) - PI_QUART;
            ybig = p * sin(xn + x);
        }

        float4v y = blend(ysmall, ybig);
        y = select(isnan(x), x, y);

        return y;
    }
}}}}
