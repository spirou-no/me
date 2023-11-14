

#include <utility/math.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	/*							expf.c
 *
 *	Exponential function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, expf();
 *
 * y = expf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns e (2.71828...) raised to the x power.
 *
 * Range reduction is accomplished by separating the argument
 * into an integer k and fraction f such that
 *
 *     x    k  f
 *    e  = 2  e.
 *
 * A polynomial is used to approximate exp(f)
 * in the basic range [-0.5, 0.5].
 *
 *
 * ACCURACY:	
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      +- MAXLOG   100000      1.7e-7      2.8e-8
 *
 *
 * Error amplification in the exponential function can be
 * a serious matter.  The error propagation involves
 * exp( X(1+delta) ) = exp(X) ( 1 + X*delta + ... ),
 * which shows that a 1 lsb error in representing X produces
 * a relative error of X times 1 lsb in the function.
 * While the routine gives an accurate result for arguments
 * that are exactly represented by a double precision
 * computer number, the result contains amplified roundoff
 * error for large arguments not exactly represented.
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * expf underflow    x < MINLOGF         0.0
 * expf overflow     x > MAXLOGF         MAXNUMF
 *
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision exponential function.
 * test interval: [-0.5, +0.5]
 * trials: 80000
 * peak relative error: 7.6e-8
 * rms relative error: 2.8e-8
 */


	float4v exp(const float4v& xin) {
		const float LOGE2A =   0.693359375;
		const float LOGE2B =  -2.12194440e-4;

		/* Express e**x = e**g 2**n
		 *   = e**g e**( n loge(2) )
		 *   = e**( g + n loge(2) )
		 */
		float4v z = round(LOG2_BASEE * xin);
		float4v x = (xin - z * LOGE2A) - z*LOGE2B;  // reduction through extended precision arithmetic
		int4v   n = round_int(x);

		/* Theoretical peak relative error in [-0.5, +0.5] is 4.2e-9. */
		const float expcof[] = {
			1.000000000,
			1.000000000,
			5.0000001201E-1,
			1.6666665459E-1,
			4.1665795894E-2,
			8.3334519073E-3,
			1.3981999507E-3,
			1.9875691500E-4
		};
		//z = x * x;
		//float4v y = eval_polynomial(x, expcof)*z + x + 1.0f;
		float4v y = eval_polynomial(x, expcof);
		y = ldexp(y, n);

		mask4v overflow = xin > MAXLOG;
		mask4v underflow= xin < MINLOG;

		float4v inf_zero  = select(overflow,  float4v::infinites());
		// float4v zero = select(underflow, float4v::zeros()); // don't need this
		y = select(overflow || underflow, inf_zero, y);

		return y;
	}



	/*							exp2f.c
	 *
	 *	Base 2 exponential function
	 *
	 *
	 *
	 * SYNOPSIS:
	 *
	 * float x, y, exp2f();
	 *
	 * y = exp2f( x );
	 *
	 *
	 *
	 * DESCRIPTION:
	 *
	 * Returns 2 raised to the x power.
	 *
	 * Range reduction is accomplished by separating the argument
	 * into an integer k and fraction f such that
	 *     x    k  f
	 *    2  = 2  2.
	 *
	 * A polynomial approximates 2**x in the basic range [-0.5, 0.5].
	 *
	 *
	 * ACCURACY:
	 *
	 *                      Relative error:
	 * arithmetic   domain     # trials      peak         rms
	 *    IEEE     -127,+127    100000      1.7e-7      2.8e-8
	 *
	 *
	 * See exp.c for comments on error amplification.
	 *
	 *
	 * ERROR MESSAGES:
	 *
	 *   message         condition      value returned
	 * exp underflow    x < -MAXL2        0.0
	 * exp overflow     x > MAXL2         MAXNUMF
	 *
	 * For IEEE arithmetic, MAXL2 = 127.
	 */


	float4v 
	exp2(const float4v& xin) {
		static const float exp2cof[] = {
			1.000000000000000,
			6.931472028550421E-001,
			2.402264791363012E-001,
			5.550332471162809E-002,
			9.618437357674640E-003,
			1.339887440266574E-003,
			1.535336188319500E-004,
		};
		const float MAXL2 = 127.0;
		const float MINL2 =-127.0;

		/* The following is necessary because range reduction blows up: */
		// really needed? 
		//if( x == 0 )
		//	return(1.0);

		/* separate into integer and fractional parts */
		float4v px	= round(xin);
		int4v   n	= round_int(px);
		float4v x	= xin - px;

		/* rational approximation
		 * exp2(x) = 1.0 +  xP(x)
		 */
		float4v y = eval_polynomial(x, exp2cof);
		y = ldexp(y, n);

		mask4v overflow = xin > MAXL2;
		mask4v underflow= xin < MINL2;

		float4v inf_zero  = select(overflow,  float4v::infinites());
		// float4v zero = select(underflow, float4v::zeros()); // don't need this
		y = select(overflow || underflow, inf_zero, y);

		return y;
	}




		/*							exp10f.c
	 *
	 *	Base 10 exponential function
	 *      (Common antilogarithm)
	 *
	 *
	 *
	 * SYNOPSIS:
	 *
	 * float x, y, exp10f();
	 *
	 * y = exp10f( x );
	 *
	 *
	 *
	 * DESCRIPTION:
	 *
	 * Returns 10 raised to the x power.
	 *
	 * Range reduction is accomplished by expressing the argument
	 * as 10**x = 2**n 10**f, with |f| < 0.5 log10(2).
	 * A polynomial approximates 10**f.
	 *
	 *
	 *
	 * ACCURACY:
	 *
	 *                      Relative error:
	 * arithmetic   domain     # trials      peak         rms
	 *    IEEE      -38,+38     100000      9.8e-8      2.8e-8
	 *
	 * ERROR MESSAGES:
	 *
	 *   message         condition      value returned
	 * exp10 underflow    x < -MAXL10        0.0
	 * exp10 overflow     x > MAXL10       MAXNUM
	 *
	 * IEEE single arithmetic: MAXL10 = 38.230809449325611792.
	 *
	 */

	/*
	Cephes Math Library Release 2.2:  June, 1992
	Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
	Direct inquiries to 30 Frost Street, Cambridge, MA 02140
	*/

	float4v exp10(const float4v& xin) {
		const float LOG210 = 3.32192809488736234787e0;
		const float LG102A = 3.00781250000000000000E-1;
		const float LG102B = 2.48745663981195213739E-4;
		const float MAXL10 = 38.230809449325611792;
		const float MINL10 = -MAXL10; // why not log10(1e-45?), i.e. smallest denorm?

		/* Express e**x = e**g 2**n
		 *   = e**g e**( n loge(2) )
		 *   = e**( g + n loge(2) )
		 */
		float4v z = round(LOG210 * xin);
		float4v x = (xin - z * LG102A) - z*LG102B;  // reduction through extended precision arithmetic
		int4v   n = round_int(z);

		/* Theoretical peak relative error in [-0.5, +0.5] is 4.2e-9. */
		const float exp10cof[] = {
			1.0000000000E+0,
			5.0000001201E-1,
			1.6666665459E-1,
			4.1665795894E-2,
			8.3334519073E-3,
			1.3981999507E-3,
			1.9875691500E-4
		};
		float4v y = eval_polynomial(x, exp10cof);
		y = ldexp(y, n);

		mask4v overflow = xin > MAXL10;
		mask4v underflow= xin < MINL10;

		float4v inf_zero  = select(overflow,  float4v::infinites());
		// float4v zero = select(underflow, float4v::zeros()); // don't need this
		y = select(overflow || underflow, inf_zero, y);

		return y;
	}

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


	float4v expm1(const float4v& x) {
		// TODO reimplement for x near 0.
		return exp(x) - 1.0f;
	}

}}}}