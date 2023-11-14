/*							asinf.c
 *
 *	Inverse circular sine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, asinf();
 *
 * y = asinf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle between -pi/2 and +pi/2 whose sine is x.
 *
 * A polynomial of the form x + x**3 P(x**2)
 * is used for |x| in the interval [0, 0.5].  If |x| > 0.5 it is
 * transformed by the identity
 *
 *    asin(x) = pi/2 - 2 asin( sqrt( (1-x)/2 ) ).
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -1, 1       100000       2.5e-7       5.0e-8
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * asinf domain        |x| > 1           0.0
 *
 */
/*							acosf()
 *
 *	Inverse circular cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, acosf();
 *
 * y = acosf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle between -pi/2 and +pi/2 whose cosine
 * is x.
 *
 * Analytically, acos(x) = pi/2 - asin(x).  However if |x| is
 * near 1, there is cancellation error in subtracting asin(x)
 * from pi/2.  Hence if x < -0.5,
 *
 *    acos(x) =	 pi - 2.0 * asin( sqrt((1+x)/2) );
 *
 * or if x > +0.5,
 *
 *    acos(x) =	 2.0 * asin(  sqrt((1-x)/2) ).
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -1, 1      100000       1.4e-7      4.2e-8
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * acosf domain        |x| > 1           0.0
 */

/*							asin.c	*/

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision circular arcsine
 * test interval: [-0.5, +0.5]
 * trials: 10000
 * peak relative error: 6.7e-8
 * rms relative error: 2.5e-8
 */

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	float4v
	asin(const float4v& xin) {

		float4v x = abs(xin);
		mask4v validmask = x <= 1.0f;

		// no specail handling of small x<1e-4 args.
		// it is an optimization, polynomial will handle it.

		mask4v bigmask = x > 0.5f;
		float4v zbig = 0.5f * (1.0f - x);
		float4v xbig = sqrt(zbig);
		float4v zsmall = x*x;

		x = select(bigmask, xbig, x);
		float4v z = select(bigmask, zbig, zsmall);

		static float asincof[] = {
			1.6666752422E-1,
			7.4953002686E-2,
			4.5470025998E-2,
			2.4181311049E-2,
			4.2163199048E-2
		};

		float4v y = eval_polynomial(x, asincof)*z*x + x;

		float4v t = PI_HALF - (y+y);
		y = select(bigmask, t, y);
		y = changesign(y, xin.reinterpret_bits());
		y = select(validmask, y, float4v::NaNs());
	
		return y;
	}


	float4v acos(const float4v& xin) {

		// no valid range check here, asinf will do it
		int4v  signbits = xin.reinterpret_bits();  // TODO not used!!

		float4v xarg   = sqrt(0.5f*(1.0f-abs(xin)));
		mask4v  bigarg = xin > 0.5f;
		mask4v  smallarg = xin < -0.5f;

		xarg = select(bigarg || smallarg, xarg, xin);

		float4v k1 = select(smallarg, float4v(PI));
		float4v k2 = select(!(bigarg || smallarg), float4v(PI_HALF));
		float4v k  = blend(k1, k2);

		float4v yasin   = asin(xarg);
		float4v a       = select(bigarg || smallarg, float4v(2.0f), 1.0f);
		a               = changesign(a, !bigarg);

		float4v y = k + a * yasin;

		return y;
	}

}}}}