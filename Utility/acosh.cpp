/*							acoshf.c
 *
 *	Inverse hyperbolic cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, acoshf();
 *
 * y = acoshf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns inverse hyperbolic cosine of argument.
 *
 * If 1 <= x < 1.5, a polynomial approximation
 *
 *	sqrt(z) * P(z)
 *
 * where z = x-1, is used.  Otherwise,
 *
 * acosh(x)  =  log( x + sqrt( (x-1)(x+1) ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      1,3         100000      1.8e-7       3.9e-8
 *    IEEE      1,2000      100000                   3.0e-8
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * acoshf domain      |x| < 1            0.0
 *
 */

/*							acosh.c	*/

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision inverse hyperbolic cosine
 * test interval: [1.0, 1.5]
 * trials: 10000
 * peak relative error: 1.7e-7
 * rms relative error: 5.0e-8
 *
 * Copyright (C) 1989 by Stephen L. Moshier.  All rights reserved.
 */

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

	float4v
	acosh(const float4v& xin) {

		float4v x = xin;
		mask4v validmask = x >= 1.0;
		mask4v hugearg   = x > 1500.0;

		float4v yhuge;
		float4v ynothuge;

		if (hugearg.any_true()) {
			yhuge = log(x) + LOG2_BASEE;
		}

		if (hugearg.any_false()) {
			float4v ymed;
			float4v ysmall;
			float4v z = x - 1.0;
			mask4v smallarg  = x < 0.5;

			if (smallarg.any_true()) {
				const float acoshcof[] = {
					 1.4142135263E0,
					-1.1784741703E-1,
					+2.6454905019E-2,
					-7.5272886713E-3,
					 1.7596881071E-3
				};

				ysmall = eval_polynomial(z, acoshcof) * sqrt(z);
			}

			if (smallarg.any_false()) {
				ymed = log(x + sqrt(z*(x+1.0)));  // log(x + sqrt((x+1)*(x-1))
			}
			ynothuge = select(smallarg, ysmall, ymed);
		}

		float4v y = select(hugearg, yhuge, ynothuge);
		y = select(x < 1.0, float4v::NaNs(), y);  // domain error?

		return y;
	}
}}}}
