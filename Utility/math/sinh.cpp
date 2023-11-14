/*							sinhf.c
 *
 *	Hyperbolic sine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, sinhf();
 *
 * y = sinhf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns hyperbolic sine of argument in the range MINLOGF to
 * MAXLOGF.
 *
 * The range is partitioned into two segments.  If |x| <= 1, a
 * polynomial approximation is used.
 * Otherwise the calculation is sinh(x) = ( exp(x) - exp(-x) )/2.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     +-MAXLOG     100000      1.1e-7      2.9e-8
 *
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision hyperbolic sine
 * test interval: [-1, +1]
 * trials: 10000
 * peak relative error: 9.0e-8
 * rms relative error: 3.0e-8
 */

#include <utility/math.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

float4v sinhf(const float4v& xin) {

	float4v x = abs(xin);
	float4v zbig;
	// no domain check, handled by exp
	mask4v bigarg = x > 1.0f;
	if (bigarg.any_true()) {
		zbig = exp(x);
		zbig = 0.5f*zbig - 0.5f/zbig;
		zbig = changesign(zbig, xin.reinterpret_bits());
	}

	float4v zsmall;
	if (bigarg.any_false()) {
		float4v z = x*x;
		static const float sinhcof[] = {
			1.66667160211E-1,
			8.33028376239E-3,
			2.03721912945E-4
		};

		zsmall = eval_polynomial(z, sinhcof) *z*x + x;
	}

	return select(bigarg, zbig, zsmall);
}

}}}}