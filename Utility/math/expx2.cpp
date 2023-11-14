/*							expx2f.c
 *
 *	Exponential of squared argument
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, expx2f();
 *
 * y = expx2f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes y = exp(x*x) while suppressing error amplification
 * that would ordinarily arise from the inexactness of the argument x*x.
 * 
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic    domain     # trials      peak       rms
 *   IEEE      -9.4, 9.4      10^7       1.7e-7     4.7e-8
 *
 */

/*
Cephes Math Library Release 2.9:  June, 2000
Copyright 2000 by Stephen L. Moshier
*/

#include <utility/math.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	

	float4v expx2f(const float4v& xin) {

		float4v x = abs(xin);

		/* Represent x as an exact multiple of 1/32 plus a residual.  */
		float4v m = .03125f * floor(32.0f * x + 0.5f);
		x -= m;

		/* x**2 = m**2 + 2mf + f**2 */
		float4v u = m * m;
		float4v u1 = (m+m + x)*x; //2 * m * x  +  x * x;

		// no overflow handling here, exp() will deal with it
		/* u is exact, u1 is small.  */
		float4v y = exp(u) * exp(u1);

		return y;
	}
}}}}


