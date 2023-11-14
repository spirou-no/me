/*							coshf.c
 *
 *	Hyperbolic cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, coshf();
 *
 * y = coshf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns hyperbolic cosine of argument in the range MINLOGF to
 * MAXLOGF.
 *
 * cosh(x)  =  ( exp(x) + exp(-x) )/2.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     +-MAXLOGF    100000      1.2e-7      2.8e-8
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * coshf overflow  |x| > MAXLOGF       MAXNUMF
 *
 *
 */

/*							cosh.c */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1985, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <Utility/math.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

float4v coshf(const float4v& xin) {

	float4v x = abs(xin);
	float4v y = exp(x);
	y += 1.0f/y;
	y *= 0.5f;

	return y;
}

}}}}
