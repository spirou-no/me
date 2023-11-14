/*							sinf.c
 *
 *	Circular sine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, sinf();
 *
 * y = sinf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of pi/4.  The reduction
 * error is nearly eliminated by contriving an extended precision
 * modular arithmetic.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the sine is approximated by
 *      x  +  x**3 P(x**2).
 * Between pi/4 and pi/2 the cosine is represented as
 *      1  -  x**2 Q(x**2).
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain      # trials      peak       rms
 *    IEEE    -4096,+4096   100,000      1.2e-7     3.0e-8
 *    IEEE    -8192,+8192   100,000      3.0e-7     3.0e-8
 * 
 * ERROR MESSAGES:
 *
 *   message           condition        value returned
 * sin total loss      x > 2^24              0.0
 *
 * Partial loss of accuracy begins to occur at x = 2^13
 * = 8192. Results may be meaningless for x >= 2^24
 * The routine as implemented flags a TLOSS error
 * for x >= 2^24 and returns 0.0.
 */

/*							cosf.c
 *
 *	Circular cosine
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, cosf();
 *
 * y = cosf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Range reduction is into intervals of pi/4.  The reduction
 * error is nearly eliminated by contriving an extended precision
 * modular arithmetic.
 *
 * Two polynomial approximating functions are employed.
 * Between 0 and pi/4 the cosine is approximated by
 *      1  -  x**2 Q(x**2).
 * Between pi/4 and pi/2 the sine is represented as
 *      x  +  x**3 P(x**2).
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain      # trials      peak         rms
 *    IEEE    -8192,+8192   100,000      3.0e-7     3.0e-8
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1985, 1987, 1988, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


/* Single precision circular sine
 * test interval: [-pi/4, +pi/4]
 * trials: 10000
 * peak relative error: 6.8e-8
 * rms relative error: 2.6e-8
 */

#include <utility/math.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	/* Note, these constants are for a 32-bit significand: */
	/*
	static float DP1 =  0.7853851318359375;
	static float DP2 =  1.30315311253070831298828125e-5;
	static float DP3 =  3.03855025325309630e-11;
	static float lossth = 65536.;
	*/

	/* These are for a 24-bit significand: */
	namespace {
		const float PI_QUART_INVERSE = 1.27323954473516;
		const float DP1			= 0.78515625;
		const float DP2			= 2.4187564849853515625e-4;
		const float DP3			= 3.77489497744594108e-8;
		const float lossth	= 8192.;
		const float MAX_ARG	= 16777215.;

		const float sincof[] = {
			-1.6666654611E-1,
			 8.3321608736E-3,
			-1.9515295891E-4
		};
	
		const float coscof[] = {
			 4.166664568298827E-002,
			-1.388731625493765E-003,
			 2.443315711809948E-005
		};


		__forceinline void
		normalize_args(float4v& x, int4v& sign, int4v& quadrant) {
			int4v   octant   = trunc_int(PI_QUART_INVERSE * x);		// integer part of x/(PI/4) 
			quadrant = (octant+1) >> 1;  // round j upwards to next multiple of 2
		
			float4v y(quadrant+quadrant);
			quadrant &= 0x3;   /* quadrant modulo 360 degrees */

			/* reflect in x axis */
			mask4v reflect = quadrant > 1;
			quadrant -= int4v(reflect) & 2;
			sign ^= reflect;  // negate all bits in signmask (also mantissa), then and with signmask to keep only sign

			//mask4v loss_precision = x > lossth;
			//float4v lost_prec = x - y*PI_QUART;
			/* Extended precision modular arithmetic */
			float4v ext_prec  = ((x - y * DP1) - y * DP2) - y * DP3;
	
			// TODO why not always extended precision?
			x = ext_prec; // select(loss_precision, lost_prec, ext_prec);
		}

		__forceinline float4v calc_cosm1(const float4v& z) {
			return (eval_polynomial(z, coscof)*z - 0.5f) * z;
		}

		__forceinline float4v calc_cos(const float4v& z) {
			return calc_cosm1(z) + 1.0f;
		}

		__forceinline float4v calc_sinc(const float4v& z) {
			return eval_polynomial(z, sincof) *z + 1.0f;
		}

    
        struct eval_cos {
            __forceinline float4v operator()(const float4v& x) const { return calc_cos(x*x); }
        };

        struct eval_sin {
            __forceinline float4v operator()(const float4v& x) const { return x*calc_sinc(x*x); }
        };

        struct eval_cosc {
            __forceinline float4v operator()(const float4v& x) const { return calc_cos(x*x)/x; }
        };

        struct eval_sinc {
            __forceinline float4v operator()(const float4v& x) const { return calc_sinc(x*x); }
        };

        struct eval_cosm1 {
            __forceinline float4v operator()(const float4v& x) const { return calc_cosm1(x*x); }
        };

        struct eval_sinm1 {
            __forceinline float4v operator()(const float4v& x) const { return x*calc_sinc(x*x) - 1.0f; }
        };
	}


#if 0
	float4v
	sin(const float4v& xin) {
		int4v  signmask  = xin.reinterpret_bits(); 						// keep sign of input, non-sign bits are don't care. thus copy input x ok.
		float4v x        = abs(xin);								
		mask4v  validmask= x < MAX_ARG;									// which inputs are valid
	
		int4v quadrant;
		normalize_args(x, signmask, quadrant);
	
		mask4v  use_cos_mask = quadrant == 1;
		float4v z = x * x;
		float4v ycos = calc_cos(z);
		float4v ysin = x * calc_sinc(z); 

		float4v y = select(use_cos_mask, ycos, ysin);
		y = changesign(y, signmask);
		y = select(validmask, y, float4v::NaNs());

		return y;
	}


	float4v cos(const float4v& xin) {
		int4v  signbits;									// default positive correct for cos
		float4v x        = abs(xin);								
		mask4v  validmask= x < MAX_ARG;						// which inputs are valid
	
		int4v quadrant;
		normalize_args(x, signbits, quadrant);
		signbits ^= quadrant > 0;							// negate again if in quadrant 1

		mask4v  use_sin_mask = quadrant == 1;
		float4v z = x * x;
		float4v ycos = calc_cos(z);
		float4v ysin = x * calc_sinc(z); 

		float4v y = select(use_sin_mask, ysin, ycos);
		y = changesign(y, signbits);
		y = select(validmask, y, float4v::NaNs());

		return y;
	}



	/*!
	Calculate cos(x)-1, avoiding cancellation error for small x.
	*/
	float4v cosm1(const float4v& xin) {
		mask4v smallarg = abs(xin) < PI_QUART;
		float4v ysmall;
		float4v ybig;

		if (smallarg.any_true()) {
			ysmall = calc_cosm1(xin*xin);
		}

		if (smallarg.any_false()) {
			ybig = cos(xin) - 1.0f;
		}

		float4v y = select(smallarg, ysmall, ybig);
		return y;
	}


	/*!
	Calculate y=sin(x)/x, 
	*/
	float4v sinc(const float4v& xin) {
		mask4v smallarg = abs(xin) < PI_QUART;
		float4v ysmall;
		float4v ybig;

		if (smallarg.any_true()) {
			ysmall = calc_sinc(xin*xin);
		}

		if (smallarg.any_false()) {
			ybig = sin(xin) / xin;
		}

		float4v y = select(smallarg, ysmall, ybig);
		return y;
	}



	void sincos(const float4v& xin, float4v& sine, float4v& cosine) {
		int4v  signbits;
		float4v x        = abs(xin);								
		mask4v  validmask= x < MAX_ARG;								// which inputs are valid
	
		int4v quadrant;
		normalize_args(x, signbits, quadrant);
	
		mask4v  use_mask = quadrant == 1;
		float4v z = x * x;
		float4v ycos = calc_cos(z);
		float4v ysin = calc_sin(x, z); 

		{
			// build sine result
			float4v y = select(use_mask, ycos, ysin);
			int4v sinsign = signbits ^ xin.reinterpret_bits(); 		// negate result if input negative
			y = changesign(y, sinsign);
			sine = select(validmask, y, float4v::NaNs());
		}

		{
			// and cosine result
			int4v cossign = signbits ^ (quadrant > 0); // cos sign = -sin sign if quadrant 1
			float4v y = select(use_mask, ysin, ycos);
			y = changesign(y, cossign);
			cosine = select(validmask, y, float4v::NaNs());
		}
	}
#else
    template <typename SIN, typename COS>
	__forceinline float4v
	sin_impl(const float4v& xin, const SIN& sinfunc, const COS& cosfunc) {
		int4v  signmask  = xin.reinterpret_bits(); 						// keep sign of input, non-sign bits are don't care. thus copy input x ok.
		float4v x        = abs(xin);								
		mask4v  validmask= x < MAX_ARG;									// which inputs are valid
	
		int4v quadrant;
		normalize_args(x, signmask, quadrant);
	
		mask4v  use_cos_mask = quadrant == 1;
		//float4v z = x * x;
		float4v ycos = cosfunc(x);
		float4v ysin = sinfunc(x); 

		float4v y = select(use_cos_mask, ycos, ysin);
		y = changesign(y, signmask);
		y = select(validmask, y, float4v::NaNs());

		return y;
    }


    
    template <typename SIN, typename COS>
	__forceinline float4v
	cos_impl(const float4v& xin, const SIN& sinfunc, const COS& cosfunc) {
		int4v  signbits;									// default positive correct for cos
		float4v x        = abs(xin);								
		mask4v  validmask= x < MAX_ARG;						// which inputs are valid
	
		int4v quadrant;
		normalize_args(x, signbits, quadrant);
		signbits ^= quadrant > 0;							// negate again if in quadrant 1

		mask4v  use_sin_mask = quadrant == 1;
		//float4v z = x * x;
		float4v ycos = cosfunc(x);
		float4v ysin = sinfunc(x); 

		float4v y = select(use_sin_mask, ysin, ycos);
		y = changesign(y, signbits);
		y = select(validmask, y, float4v::NaNs());

		return y;
	}


    float4v sin(const float4v& x) {
        return sin_impl(x, eval_sin(), eval_cos());
    }


    float4v cos(const float4v& x) {
        return cos_impl(x, eval_sin(), eval_cos());
    }


	/*!
	Calculate cos(x)-1, avoiding cancellation error for small x.
	*/
	float4v cosm1(const float4v& xin) {
        return cos_impl(xin, eval_sinm1(), eval_cosm1());
    }


	/*!
	Calculate y=sin(x)/x, 
	*/
	float4v sinc(const float4v& xin) {
        return sin_impl(xin, eval_sinc(), eval_cosc());
	}



	void sincos(const float4v& xin, float4v& sine, float4v& cosine) {
		int4v  signbits;
		float4v x        = abs(xin);								
		mask4v  validmask= x < MAX_ARG;								// which inputs are valid
	
		int4v quadrant;
		normalize_args(x, signbits, quadrant);
	
		mask4v  use_mask = quadrant == 1;
		float4v z = x * x;
		float4v ycos = calc_cos(z);
		float4v ysin = x*calc_sinc(z); 

		{
			// build sine result
			float4v y = select(use_mask, ycos, ysin);
			int4v sinsign = signbits ^ xin.reinterpret_bits(); 		// negate result if input negative
			y = changesign(y, sinsign);
			sine = select(validmask, y, float4v::NaNs());
		}

		{
			// and cosine result
			int4v cossign = signbits ^ (quadrant > 0); // cos sign = -sin sign if quadrant 1
			float4v y = select(use_mask, ysin, ycos);
			y = changesign(y, cossign);
			cosine = select(validmask, y, float4v::NaNs());
		}
	}

#endif
}}}}
