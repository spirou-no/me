/*							tanhf.c
 *
 *	Hyperbolic tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, tanhf();
 *
 * y = tanhf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns hyperbolic tangent of argument in the range MINLOG to
 * MAXLOG.
 *
 * A polynomial approximation is used for |x| < 0.625.
 * Otherwise,
 *
 *    tanh(x) = sinh(x)/cosh(x) = 1  -  2/(exp(2x) + 1).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -2,2        100000      1.3e-7      2.6e-8
 *
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision hyperbolic tangent
 * test interval: [-0.625, +0.625]
 * trials: 10000
 * peak relative error: 7.2e-8
 * rms relative error: 2.6e-8
 */

#include <utility/math.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	float4v tanh(const float4v& xin) {

		float4v x = abs(xin);
		mask4v hugearg = x > float(0.5*MAXLOG);
		// result is +1 or -1
		float4v yhuge = copysign(1.0, xin);

		mask4v smallarg = x < 0.625;
		float4v ybig;
		float4v ysmall;

		if (smallarg.any_true()) {
			const float tanhcof[] = {
			   -3.33332819422E-1,
				1.33314422036E-1,
			   -5.37397155531E-2,
				2.06390887954E-2,
			   -5.70498872745E-3
			};

			float4v z = x * x;
			ysmall = eval_polynomial(z, tanhcof) *z*xin + xin;
		}

		if (smallarg.any_false()) {
			ybig = exp(x+x);
			ybig = 1.0 - 2.0/(ybig + 1.0);
			ybig = changesign(ybig, xin);  // negate if xin is negative
		}

		// process results, make output
		float4v y = select(smallarg, ysmall, ybig);
		y = select(hugearg, yhuge, y);
		y = select(isnan(xin), xin, y);  // needed? NaNs should propagate thru equations

		return y;
	}
}}}}
