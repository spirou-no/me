/*							atanhf.c
 *
 *	Inverse hyperbolic tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, atanhf();
 *
 * y = atanhf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns inverse hyperbolic tangent of argument in the range
 * MINLOGF to MAXLOGF.
 *
 * If |x| < 0.5, a polynomial approximation is used.
 * Otherwise,
 *        atanh(x) = 0.5 * log( (1+x)/(1-x) ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -1,1        100000      1.4e-7      3.1e-8
 *
 */

/*						atanh.c	*/


/*
Cephes Math Library Release 2.2:  June, 1992
Copyright (C) 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision inverse hyperbolic tangent
 * test interval: [-0.5, +0.5]
 * trials: 10000
 * peak relative error: 8.2e-8
 * rms relative error: 3.0e-8
 */
#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

	float4v atanh(const float4v& xin) {

		float4v x = xin;
		mask4v smallarg = abs(x) < 0.5;
		float4v ysmall;
		float4v ybig;
		if (smallarg.any_true()) {
			const float atanhcof[] = {
				3.33337300303E-1,
				1.99782164500E-1,
				1.46691431730E-1,
				8.24370301058E-2,
				1.81740078349E-1
			};

			float4v xx = x*x;
			ysmall = eval_polynomial(xx, atanhcof) * xx * x   + x;
		}

		if (smallarg.any_false()) {
			ybig = 0.5 * log((1.0+x) / (1.0-x));
			ybig = select(x == 1.0, float4v::infinites(), ybig);
			ybig = select(x >  1.0, float4v::NaNs(), ybig);
		}

		float4v y = select(smallarg, ysmall, ybig);

		return y;
	}
}}}}