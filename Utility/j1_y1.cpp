/*							j1f.c
 *
 *	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, j1f();
 *
 * y = j1f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 2] and
 * (2, infinity). In the first interval a polynomial approximation
 *        2 
 * (w - r  ) x P(w)
 *       1  
 *                     2 
 * is used, where w = x  and r is the first zero of the function.
 *
 * In the second interval, the modulus and phase are approximated
 * by polynomials of the form Modulus(x) = sqrt(1/x) Q(1/x)
 * and Phase(x) = x + 1/x R(1/x^2) - 3pi/4.  The function is
 *
 *   j0(x) = Modulus(x) cos( Phase(x) ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak       rms
 *    IEEE      0,  2       100000       1.2e-7     2.5e-8
 *    IEEE      2, 32       100000       2.0e-7     5.3e-8
 *
 *
 */
/*							y1.c
 *
 *	Bessel function of second kind of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y1();
 *
 * y = y1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind of order one
 * of the argument.
 *
 * The domain is divided into the intervals [0, 2] and
 * (2, infinity). In the first interval a rational approximation
 * R(x) is employed to compute
 *
 *                  2
 * y0(x)  =  (w - r  ) x R(x^2)  +  2/pi (ln(x) j1(x) - 1/x) .
 *                 1
 *
 * Thus a call to j1() is required.
 *
 * In the second interval, the modulus and phase are approximated
 * by polynomials of the form Modulus(x) = sqrt(1/x) Q(1/x)
 * and Phase(x) = x + 1/x S(1/x^2) - 3pi/4.  Then the function is
 *
 *   y0(x) = Modulus(x) sin( Phase(x) ).
 *
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    IEEE      0,  2       100000       2.2e-7     4.6e-8
 *    IEEE      2, 32       100000       1.9e-7     5.3e-8
 *
 * (error criterion relative when |y1| > 1).
 *
 */


/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    const float JP[] = {
        -3.405537384615824E-002f,
         1.937383947804541E-003f,
        -4.541343896997497E-005f,
         6.009061827883699E-007f,
        -4.878788132172128E-009f
    };

    const float YP[] = {
         4.202369946500099E-002f,
        -2.641785726447862E-003f,
         6.719543806674249E-005f,
        -9.496460629917016E-007f,
         8.061978323326852E-009f
    };

    const float MO1[] = {
         7.978845453073848E-001f,
         4.976029650847191E-006f,
         1.493389585089498E-001f,
         5.435364690523026E-003f,
        -2.102302420403875E-001f,
         3.138238455499697E-001f,
        -2.284801500053359E-001f,
         6.913942741265801E-002f
    };

    const float PH1[] = {
         3.749989509080821E-001f,
        -1.637986776941202E-001f,
         3.503787691653334E-001f,
        -1.544842782180211E+000f,
         7.222973196770240E+000f,
        -2.485774108720340E+001f,
         5.073465654089319E+001f,
        -4.497014141919556E+001f
    };

    static float YO1 =  4.66539330185668857532f;
    static float Z1 = 1.46819706421238932572E1f;

    static float THPIO4F =  2.35619449019234492885f;    /* 3*pi/4 */
    static float TWOOPI =  0.636619772367581343075535f; /* 2/pi */



    float4v j1(const float4v& xin) {

        mask4v nanarg = isnan(xin);
        float4v y = select(nanarg, xin);

        float4v x = abs(xin);

        mask4v less2arg = x <= 2.0f;

        if (less2arg.any_true()) {
            float4v x2 = x*x;
            float4v p = (x2-Z1) * x * eval_polynomial(x2, JP);
            y = select(less2arg, p);
        }

        mask4v openarg = !(nanarg || less2arg);
        if (openarg.any_true()) {
            float4v q = 1.0f/x;
            float4v w = sqrt(q);
            float4v p = w * eval_polynomial(q, MO1);
            float4v xn = q * eval_polynomial(q*q, PH1) - THPIO4F;
            float4v yn = p * cos(xn+x);

            y = select(openarg, yn, y);
        }

        return y;
    }


    
    float4v y1(const float4v& x) {

        mask4v nanarg = isnan(x);
        float4v y = x;

        mask4v negarg = x <= 0.0f;
        y = select(negarg, float4v::infinites(), y);

        mask4v openarg = !(nanarg || negarg);
        mask4v less2arg = x <= 2.0f  && openarg;

        if (openarg.any_true()) {
            float4v x2 = x*x;
            float4v w = (x2 - YO1) * x * eval_polynomial(x2, YP);
            w += TWOOPI * (j1(x) * log(x)  -  1.0f/x);
            y = select(openarg, w, y);
            openarg = openarg && !less2arg;
        }

        if (openarg.any_true()) {
            // TODO consolidate same steps as j1, cos instead of sin
            float4v q = 1.0f/x;
            float4v w = sqrt(q);
            float4v p = w * eval_polynomial(q, MO1);
            float4v xn = q * eval_polynomial(q*q, PH1) - THPIO4F;
            float4v yn = p * sin(xn+x);

            y = select(openarg, yn, y);
        }

        return y;
    }
}}}}
