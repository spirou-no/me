/*							asinhf.c
 *
 *	Inverse hyperbolic sine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, asinhf();
 *
 * y = asinhf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns inverse hyperbolic sine of argument.
 *
 * If |x| < 0.5, the function is approximated by a rational
 * form  x + x**3 P(x)/Q(x).  Otherwise,
 *
 *     asinh(x) = log( x + sqrt(1 + x*x) ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -3,3        100000       2.4e-7      4.1e-8
 *
 */

/*						asinh.c	*/

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision inverse hyperbolic sine
 * test interval: [-0.5, +0.5]
 * trials: 10000
 * peak relative error: 8.8e-8
 * rms relative error: 3.2e-8
 */

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

	float4v asinh(const float4v& xin) {
	
		float4v x = abs(xin);
		float4v ybig;

		const float asinhcof[] = {
			-1.6666288134E-1,
			 7.4847586088E-2,
			-4.2699340972E-2,
			 2.0122003309E-2
		};

		float4v z = x * x;
		float4v ysmall = eval_polynomial(z, asinhcof) * z*x + x;
		mask4v smallarg = x < 0.5;

		if (smallarg.any_false()) {
			// y = log(x) + log(2)  if x > 1500
			// y = log(x+sqrt(x*x + 1)) if 0.5<x<1500
			mask4v hugearg = x > 1500.0;
			float4v arg = select(hugearg, x, x+sqrt(z+1.0));
			ybig = log(arg);
			ybig += select(hugearg, float4v(LOG2_BASEE));
		}

		float4v y = select(smallarg, ysmall, ybig);
		y = copysign(y, xin);

		return y;
	}
}}}}
