/*							i1f.c
 *
 *	Modified Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, i1f();
 *
 * y = i1f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order one of the
 * argument.
 *
 * The function is defined as i1(x) = -i j1( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       100000      1.5e-6      1.6e-7
 *
 *
 */
/*							i1ef.c
 *
 *	Modified Bessel function of order one,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, i1ef();
 *
 * y = i1ef( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of order one of the argument.
 *
 * The function is defined as i1(x) = -i exp(-|x|) j1( ix ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       30000       1.5e-6      1.5e-7
 * See i1().
 *
 */

/*							i1.c 2		*/


/*
Cephes Math Library Release 2.0:  March, 1987
Copyright 1985, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/
#include <utility/math.h>
#include <utility/sse.h>
//#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	/* Chebyshev coefficients for exp(-x) I1(x) / x
	 * in the interval [0,8].
	 *
	 * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
	 */

	static const float A[] = {
		 2.52587186443633654823E-1f,
		-1.76416518357834055153E-1f,
		 1.02643658689847095384E-1f,
		-5.29459812080949914269E-2f,
		 2.47264490306265168283E-2f,
		-1.05640848946261981558E-2f,
		 4.15642294431288815669E-3f,
		-1.51357245063125314899E-3f,
		 5.12285956168575772895E-4f,
		-1.61760815825896745588E-4f,
		 4.78156510755005422638E-5f,
		-1.32731636560394358279E-5f,
		 3.47025130813767847674E-6f,
		-8.56872026469545474066E-7f,
		 2.00329475355213526229E-7f,
		-4.44505912879632808065E-8f,
		 9.38153738649577178388E-9f
	};


	/* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
	 * in the inverted interval [8,infinity].
	 *
	 * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
	 */

	static const float B[] = {
		 7.78576235018280120474E-1f
		-9.76109749136146840777E-3f,
		-1.10588938762623716291E-4f,
		-3.88256480887769039346E-6f,
		-2.51223623787020892529E-7f,
		-2.63146884688951950684E-8f,
		-3.83538038596423702205E-9f
	 };

/*							i1.c	*/

	float4v 
	i1(const float4v& xin) {
		
		float4v x = abs(xin);
		float4v y = exp(x) * i1e(x);

		return y;
	}

	
	
	/*	i1e()	
	i1e has same body as i0e, only difference is A & B coefficients
	*/

	float4v 
	i1e(const float4v& xin) {

		float4v x = abs(x);
		mask4v smallarg = x <= 8.0;
		mask4v bigarg = x > 8.0;
		
		float4v zsmall;
		float4v zbig;

		if (smallarg.any_true()) {
			float4v y = 0.5f*x - 2.0f;
			zsmall = eval_chebychev(y, A);
		}

		if (bigarg.any_true()) {
			float4v y = 32.0f/x - 2.0f;
			zbig = eval_chebychev(y, B) / sqrt(x);
		}

		float4v result = blend(zsmall, zbig);
		result = copysign(result, xin);
		result = select(isnan(x), xin, result);

		return result;
	}
 }}}}