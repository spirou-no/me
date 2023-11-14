/*							tanf.c
 *
 *	Circular tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, tanf();
 *
 * y = tanf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular tangent of the radian argument x.
 *
 * Range reduction is modulo pi/4.  A polynomial approximation
 * is employed in the basic interval [0, pi/4].
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     +-4096        100000     3.3e-7      4.5e-8
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * tanf total loss   x > 2^24              0.0
 *
 */
/*							cotf.c
 *
 *	Circular cotangent
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, cotf();
 *
 * y = cotf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the circular cotangent of the radian argument x.
 * A common routine computes either the tangent or cotangent.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     +-4096        100000     3.0e-7      4.5e-8
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition          value returned
 * cot total loss   x > 2^24                0.0
 * cot singularity  x = 0                  MAXNUMF
 *
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision circular tangent
 * test interval: [-pi/4, +pi/4]
 * trials: 10000
 * peak relative error: 8.7e-8
 * rms relative error: 2.8e-8
 */
#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	namespace {
		// same constants as sincos
		const float DP1 = 0.78515625;
		const float DP2 = 2.4187564849853515625e-4;
		const float DP3 = 3.77489497744594108e-8;
		const float PI_QUART_INVERSE = 1.27323954473516;  /* 4/pi */
		const float lossth = 8192.;
		const float T24M1 = 16777215.;


		void tancot(const float4v& xin, float4v& yout, int4v& octant_out) {
		
			float4v x = abs(xin);
			mask4v hugearg = xin > T24M1;

			/* compute x mod PIO4 */
			int4v octant = trunc_int(PI_QUART_INVERSE * x); /* integer part of x/(PI/4) */

			/* map zeros and singularities to origin */
			octant = ((octant + 1) >> 1) << 1;  // same as octant += octant&1

			const float tancotcof[] = {
				3.33331568548E-1,
				1.33387994085E-1,
				5.34112807005E-2,
				2.44301354525E-2,
				3.11992232697E-3,
				9.38540185543E-3
			};

			float4v y(octant);
			float4v z = ((x - y * DP1) - y * DP2) - y * DP3;
			float4v zz = z * z;
			y = eval_polynomial(zz, tancotcof) * zz * z  + z;

			octant_out = octant;
			yout = y;
		}
	}



	float4v tan(const float4v& xin) {
		
		float4v y;
		int4v octant;
		tancot(xin, y, octant);  // cotflg==0

		y = select((octant & 0x2) == 0, y, -1.0/y);
		
		return copysign(y, xin);
	}

	float4v cot(const float4v& xin) {

		float4v y;
		int4v octant;
		tancot(xin, y, octant);  // cotflg==0

		y = select((octant & 0x2) == 0, 1.0/y, -y);
		y = select(xin == 0.0, float4v::infinites(), y);  // is this needed?  check return from tancot
		
		return copysign(y, xin);
	}
}}}}
