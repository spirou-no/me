/*							fresnlf.c
 *
 *	Fresnel integral
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, S, C;
 * void fresnlf();
 *
 * fresnlf( x, _&S, _&C );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the Fresnel integrals
 *
 *           x
 *           -
 *          | |
 * C(x) =   |   cos(pi/2 t**2) dt,
 *        | |
 *         -
 *          0
 *
 *           x
 *           -
 *          | |
 * S(x) =   |   sin(pi/2 t**2) dt.
 *        | |
 *         -
 *          0
 *
 *
 * The integrals are evaluated by power series for small x.
 * For x >= 1 auxiliary functions f(x) and g(x) are employed
 * such that
 *
 * C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
 * S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
 *
 *
 *
 * ACCURACY:
 *
 *  Relative error.
 *
 * Arithmetic  function   domain     # trials      peak         rms
 *   IEEE       S(x)      0, 10       30000       1.1e-6      1.9e-7
 *   IEEE       C(x)      0, 10       30000       1.1e-6      2.0e-7
 */

/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/
#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    /* S(x) for small x */
    const float sn[] = {
         5.235987735681432E-001,
        -9.228055941124598E-002,
         7.244727626597022E-003,
        -3.120693124703272E-004,
         8.424748808502400E-006,
        -1.522754752581096E-007,
         1.647629463788700E-009
    };

    /* C(x) for small x */
    const float cn[] = {
         9.999999760004487E-001,
        -2.467398198317899E-001,
         2.818489036795073E-002,
        -1.604381798862293E-003,
         5.387223446683264E-005,
        -1.157231412229871E-006,
         1.416802502367354E-008
    };


    /* Auxiliary function f(x) */
    const float fn[] = {
         2.999401847870011E+000,
        -1.032877601091159E+002,
         8.560515466275470E+003,
        -8.732356681548485E+005,
         7.343848463587323E+007,
        -4.158143148511033E+009,
         1.355942388050252E+011,
        -1.903009855649792E+012
    };

    /* Auxiliary function g(x) */
    const float gn[] = {
         9.999841934744914E-001,
        -1.493439396592284E+001,
         8.602931494734327E+002,
        -7.787789623358162E+004,
         6.492611570598858E+006,
        -3.779387713202229E+008,
         1.278350673393208E+010,
        -1.860843997624650E+011
    };



    void fresnlf(const float4v& xin, float4v& ssa, float4v& cca) {

        float4v x = abs(xin);
        float4v x2 = x*x;
        float4v ss, cc;

        mask4v nanarg = isnan(xin);
        float4v y = select(nanarg, xin);
        mask4v openarg = nanarg;

        mask4v smallarg = (x2 < 2.5625f) && openarg;
        if (smallarg.any_true()) {
            float4v t = x2 * x2;
            ss = x * x2 * eval_polynomial(t, sn);
            cc = x * eval_polynomial(t, cn);

            openarg = openarg && !smallarg;
            ss = select(smallarg, ss);
            cc = select(smallarg, cc);
        }

        if (openarg.any_true()) {
            /*		Asymptotic power series auxiliary functions
            *		for large argument
            */
            float4v t = PI * x2;
            float4v u = 1.0f/(t*t);    // 1.0/(t*t)
            t = 1.0f/t;              // 1.0/t
            float4v f = 1.0f  -  u * eval_polynomial(u, fn);
            float4v g = t * eval_polynomial(u, gn);

            float4v s, c;
            sincos(PI_HALF * x2, s, c);

            t = PI * x;
	        float4v cct = 0.5f  +  (f * s  -  g * c)/t;
	        float4v sst = 0.5f  -  (f * c  +  g * s)/t;

            ss = blend(ss, select(openarg, sst));
            cc = blend(cc, select(openarg, cct));

            mask4v hugearg = x > 36974.0f;
            float4v half = 0.5f;

            ss = select(hugearg, half, ss);
            cc = select(hugearg, half, cc);
	    }
    
        cca = copysign(cc, xin);
        ssa = copysign(ss, xin);
    }
}}}}
