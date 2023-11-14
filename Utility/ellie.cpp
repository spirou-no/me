/*							ellief.c
 *
 *	Incomplete elliptic integral of the second kind
 *
 *
 *
 * SYNOPSIS:
 *
 * float phi, m, y, ellief();
 *
 * y = ellief( phi, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Approximates the integral
 *
 *
 *                phi
 *                 -
 *                | |
 *                |                   2
 * E(phi\m)  =    |    sqrt( 1 - m sin t ) dt
 *                |
 *              | |    
 *               -
 *                0
 *
 * of amplitude phi and modulus m, using the arithmetic -
 * geometric mean algorithm.
 *
 *
 *
 * ACCURACY:
 *
 * Tested at random arguments with phi in [0, 2] and m in
 * [0, 1].
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,2        10000       4.5e-7      7.4e-8
 *
 *
 */


/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/*	Incomplete elliptic integral of second kind	*/

#include <utility/math.h>
#include <utility/math/constants.h>
#include <cmath>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    float4v 
    ellie(const float4v& phi_in, float ma) {

        float4v phi = abs(phi_in);

        if (ma == 0.0f) {
            return phi_in;
        }

        if (ma == 1.0f) {
	        return sin(phi_in);
        }

        float a = 1.0f;
        float b = std::sqrt(1.0f - ma);
        float c = std::sqrt(ma);
        float4v e = 0.0f;
        float4v t = tan(phi);
        float4v mod = trunc((phi + PI_HALF)/PI);
        int     d = 1;

        while((abs(c/a) > MACHEPF).any_true()) {
	        float temp = b/a;
	        phi += atan(t*temp) + float4v(mod) * PI;
	        mod = trunc((phi + PI_HALF)/PI);
	        t = t * (1.0f + temp)/(1.0 - temp*t*t);
	        c = 0.5f * (a - b);
	        temp = std::sqrt(a * b);
	        a = 0.5f * (a + b);
	        b = temp;
	        d += d;
	        e += c * sin(phi);
	    }

        b = 1.0f - ma;
        float4v y = ellpe(b)/ellpk(b);
        y *= (atan(t) + mod * PI)/(float(d) * a);
        y += e;
        y = copysign(y, phi_in);
        y = select(isnan(phi_in), phi_in, y);

        return y;
    }
 }}}}
