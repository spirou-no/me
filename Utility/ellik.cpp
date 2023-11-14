/*							ellikf.c
 *
 *	Incomplete elliptic integral of the first kind
 *
 *
 *
 * SYNOPSIS:
 *
 * float phi, m, y, ellikf();
 *
 * y = ellikf( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *
 *                phi
 *                 -
 *                | |
 *                |           dt
 * F(phi\m)  =    |    ------------------
 *                |                   2
 *              | |    sqrt( 1 - m sin t )
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random points with phi in [0, 2] and m in
 * [0, 1].
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,2         10000       2.9e-7      5.8e-8
 *
 *
 */


/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/*	Incomplete elliptic integral of first kind	*/

#include <utility/math.h>
#include <utility/math/constants.h>
#include <cmath>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    float4v
    ellik(const float4v& phi_in, float ma) {

        float4v phi = abs(phi_in);

        if (ma == 0.0f) {
	        return phi_in;
        }

        if (ma == 1.0f) {
	        return log(tan(0.5f*(PI_HALF + phi)));
        }

        float a = 1.0f;
        float b = std::sqrt(1.0f - ma);
        float c = std::sqrt(ma);
        int   d = 1;
        float4v t = tan(phi);
        float4v mod = trunc((phi + PI_HALF)/PI);  //  INV_PI;  OR INT???

        // NaN???
        while (std::abs(c/a) > MACHEPF) {
            float temp = b/a;
        	phi += atan(t*temp) + mod * PI;
	        mod = trunc((phi + PI_HALF)/PI);
	        t = t * (1.0f + temp)/(1.0 - temp*t*t);
	        c = (a - b) * 0.5f;
	        temp = std::sqrt(a * b);
	        a = (a + b) * 0.5f;
	        b = temp;
	        d += d;
        }

        float4v y = (atan(t) + mod * PI)/(a * float4v(float(d)));
        y = copysign(y, phi_in);
        y = select(isnan(phi), phi_in, y);

        return y;
    }

    float4v
    ellik(const float4v& /*phia*/, const float4v& /*ma*/ ) { return float4v::NaNs(); }

 }}}}
