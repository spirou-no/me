#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {


/*							logf.c
 *
 *	Natural logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, logf();
 *
 * y = logf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of x.
 *
 * The argument is separated into its exponent and fractional
 * parts.  If the exponent is between -1 and +1, the logarithm
 * of the fraction is approximated by
 *
 *     log(1+x) = x - 0.5 x**2 + x**3 P(x)
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0.5, 2.0    100000       7.6e-8     2.7e-8
 *    IEEE      1, MAXNUMF  100000                  2.6e-8
 *
 * In the tests over the interval [1, MAXNUM], the logarithms
 * of the random arguments were uniformly distributed over
 * [0, MAXLOGF].
 *
 * ERROR MESSAGES:
 *
 * logf singularity:  x = 0; returns MINLOG
 * logf domain:       x < 0; returns MINLOG
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision natural logarithm
 * test interval: [sqrt(2)/2, sqrt(2)]
 * trials: 10000
 * peak relative error: 7.1e-8
 * rms relative error: 2.7e-8
 */

	__forceinline float4v
	normalize_arg(const float4v& xin, int4v& expnt) {
		float4v x = frexp(xin, expnt);

		// expnt-- if x<sqrt(2)/2
		expnt -= select(x < SQRT2_HALF, int4v(1));
		x += select(x < SQRT2_HALF, x);
		x -= 1.0f;

		return x;
	}

	
	
	float4v 
	log(const float4v& xin) {
		
		int4v expnt;
		float4v x = normalize_arg(xin, expnt);

		const float logcof[] = {
			 3.3333331174E-1,
			-2.4999993993E-1,
			 2.0000714765E-1,
			-1.6668057665E-1,
			 1.4249322787E-1,
			-1.2420140846E-1,
			 1.1676998740E-1,
			-1.1514610310E-1,
			 7.0376836292E-2
		};

		float4v z = x*x;
		float4v y = eval_polynomial(x, logcof) * z;

		y -= -2.12194440e-4 * float4v(expnt);
		y -= 0.5 * z;				// y - 0.5 x^2 
		y += x;
		y += 0.693359375 * float4v(expnt);

		y = select(isinf(xin), xin, y);
		y = select(xin == 0.0f, -float4v::infinites(), y);  // singularity at x==0  --> -infinity
		y = select(xin < 0.0f , float4v::NaNs(), y);				// domain error for x<0 --> NaN

		return y;
	}


/*							log10f.c
 *
 *	Common logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, log10f();
 *
 * y = log10f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns logarithm to the base 10 of x.
 *
 * The argument is separated into its exponent and fractional
 * parts.  The logarithm of the fraction is approximated by
 *
 *     log(1+x) = x - 0.5 x**2 + x**3 P(x).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0.5, 2.0    100000      1.3e-7      3.4e-8
 *    IEEE      0, MAXNUMF  100000      1.3e-7      2.6e-8
 *
 * In the tests over the interval [0, MAXNUM], the logarithms
 * of the random arguments were uniformly distributed over
 * [-MAXL10, MAXL10].
 *
 * ERROR MESSAGES:
 *
 * log10f singularity:  x = 0; returns -MAXL10
 * log10f domain:       x < 0; returns -MAXL10
 * MAXL10 = 38.230809449325611792
 */

/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


/* Coefficients for log(1+x) = x - x**2/2 + x**3 P(x)/Q(x)
 * 1/sqrt(2) <= x < sqrt(2)
 */
	const float log10cof[] = {
		3.3333331174E-1,
	   -2.4999993993E-1,
		2.0000714765E-1,
	   -1.6668057665E-1,
		1.4249322787E-1,
	   -1.2420140846E-1,
		1.1676998740E-1,
	   -1.1514610310E-1,
		7.0376836292E-2
	};

	const float L102A = 3.0078125E-1;
	const float L102B = 2.48745663981195213739E-4;
	const float L10EA = 4.3359375E-1;
	const float L10EB = 7.00731903251827651129E-4;

	const float MAXL10 = 38.230809449325611792;



	float4v 
	log10(const float4v& xin) {

		int4v expnt;
		float4v x = normalize_arg(xin, expnt);
		float4v z = x*x;
		float4v y = x * (z * eval_polynomial(x, log10cof));
		y = y - 0.5 * z;   /*  y - 0.5 * x**2  */

		/* multiply log of fraction by log10(e)
		 * and base 2 exponent by log10(2)
		 */
		z = (x + y) * L10EB;  /* accumulate terms in order of size */
		z += y * L10EA;
		z += x * L10EA;
		x = float4v(expnt);
		z += x * L102B;
		z += x * L102A;

		y = select(xin == 0.0f, -float4v::infinites(), z);  // singularity at x==0  --> -infinity
		y = select(xin < 0.0f , float4v::NaNs(), z);				// domain error for x<0 --> NaN

		return y;
	}

/*							log2f.c
 *
 *	Base 2 logarithm
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, log2f();
 *
 * y = log2f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base 2 logarithm of x.
 *
 * The argument is separated into its exponent and fractional
 * parts.  If the exponent is between -1 and +1, the base e
 * logarithm of the fraction is approximated by
 *
 *     log(1+x) = x - 0.5 x**2 + x**3 P(x)/Q(x).
 *
 * Otherwise, setting  z = 2(x-1)/x+1),
 * 
 *     log(x) = z + z**3 P(z)/Q(z).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      exp(+-88)   100000      1.1e-7      2.4e-8
 *    IEEE      0.5, 2.0    100000      1.1e-7      3.0e-8
 *
 * In the tests over the interval [exp(+-88)], the logarithms
 * of the random arguments were uniformly distributed.
 *
 * ERROR MESSAGES:
 *
 * log singularity:  x = 0; returns MINLOGF/log(2)
 * log domain:       x < 0; returns MINLOGF/log(2)
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


/* Coefficients for log(1+x) = x - x**2/2 + x**3 P(x)
 * 1/sqrt(2) <= x < sqrt(2)
 */


	float4v 
	log2(const float4v& xin) {
		const float LOG2EA = 0.44269504088896340735992;
		const float log2cof[] = {
			 3.3333331174E-1,
			-2.4999993993E-1,
			 2.0000714765E-1,
			-1.6668057665E-1,
			 1.4249322787E-1,
			-1.2420140846E-1,
			 1.1676998740E-1,
			-1.1514610310E-1,
			 7.0376836292E-2,
		};

		int4v expnt;
		float4v x = normalize_arg(xin, expnt);

		float4v z = x*x;
		float4v y = x * (z * eval_polynomial(x, log2cof));
		y = y - 0.5 * z;   /*  y - 0.5 * x**2  */

		/* Multiply log of fraction by log2(e)
		 * and base 2 exponent by 1
		 *
		 * ***CAUTION***
		 *
		 * This sequence of operations is critical and it may
		 * be horribly defeated by some compiler optimizers.
		 */
		z = y * LOG2EA;
		z += x * LOG2EA;
		z += y;
		z += x;
		z += float4v(expnt);

		y = select(xin == 0.0f, -float4v::infinites(), z);  // singularity at x==0  --> -infinity
		y = select(xin < 0.0f , float4v::NaNs(), z);				// domain error for x<0 --> NaN

		return y;
	}



    //! Calculate log(1+x), preserving accuracy for small x
    float4v log1p(const float4v& x) {
        // Beware of compiler optimization...
        return (log(1.0+x) * x) / ((1.0f - x) - 1.0f);
    }

}}}}