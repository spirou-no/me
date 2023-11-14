/*							powf.c
 *
 *	Power function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, z, powf();
 *
 * z = powf( x, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Computes x raised to the yth power.  Analytically,
 *
 *      x**y  =  exp( y log(x) ).
 *
 * Following Cody and Waite, this program uses a lookup table
 * of 2**-i/16 and pseudo extended precision arithmetic to
 * obtain an extra three bits of accuracy in both the logarithm
 * and the exponential.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 *  arithmetic  domain     # trials      peak         rms
 *    IEEE     -10,10      100,000      1.4e-7      3.6e-8
 * 1/10 < x < 10, x uniformly distributed.
 * -10 < y < 10, y uniformly distributed.
 *
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * powf overflow     x**y > MAXNUMF     MAXNUMF
 * powf underflow   x**y < 1/MAXNUMF      0.0
 * powf domain      x<0 and y noninteger  0.0
 *
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/


#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>
#include <cmath>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

  /*
  Approximate pow http://martin.ankerl.com/2012/01/25/optimized-approximative-pow-in-c-and-cpp/
  */
	namespace {
		/* log2(e) - 1 */
		const float LOG2EA = 0.44269504088896340736F;

		/* 2^(-i/16)
		* The decimal values are rounded to 24-bit precision
		*/
		const float A[] = {
			1.00000000000000000000E0,
			9.57603275775909423828125E-1,
			9.17004048824310302734375E-1,
			8.78126084804534912109375E-1,
			8.40896427631378173828125E-1,
			8.05245161056518554687500E-1,
			7.71105408668518066406250E-1,
			7.38413095474243164062500E-1,
			7.07106769084930419921875E-1,
			6.77127778530120849609375E-1,
			6.48419797420501708984375E-1,
			6.20928883552551269531250E-1,
			5.94603538513183593750000E-1,
			5.69394290447235107421875E-1,
			5.45253872871398925781250E-1,
			5.22136867046356201171875E-1,
			5.00000000000000000000E-1
		};
		/* continuation, for even i only
		 * 2^(i/16)  =  A[i] + B[i/2]
		 */
		const float B[] = {
			 0.00000000000000000000E0,
			-5.61963907099083340520586E-9,
			-1.23776636307969995237668E-8,
			 4.03545234539989593104537E-9,
			 1.21016171044789693621048E-8,
			-2.00949968760174979411038E-8,
			 1.89881769396087499852802E-8,
			-6.53877009617774467211965E-9,
			 0.00000000000000000000E0
		};

		/* 1 / A[i]
		 * The decimal values are full precision
		 */
		const float Ainv[] = {
			 1.00000000000000000000000E0,
			 1.04427378242741384032197E0,
			 1.09050773266525765920701E0,
			 1.13878863475669165370383E0,
			 1.18920711500272106671750E0,
			 1.24185781207348404859368E0,
			 1.29683955465100966593375E0,
			 1.35425554693689272829801E0,
			 1.41421356237309504880169E0,
			 1.47682614593949931138691E0,
			 1.54221082540794082361229E0,
			 1.61049033194925430817952E0,
			 1.68179283050742908606225E0,
			 1.75625216037329948311216E0,
			 1.83400808640934246348708E0,
			 1.91520656139714729387261E0,
			 2.00000000000000000000000E0
		};
	}
//#define MEXP 2048.0
//#define MNEXP -2400.0


	namespace {
		/* Find a multiple of 1/16 that is within 1/16 of x. */
		__forceinline float4v reduce(const float4v& x) { return 0.0625 * floor(16.0 * x); }

//#define PERF_TEST


		__forceinline float4v
		load_index(const float data[], const int4v& ix) {
#ifdef PERF_TEST
				return float4v(data, sse::unaligned);
#else
				return float4v(data[ix.at<0>()], data[ix.at<1>()], data[ix.at<2>()], data[ix.at<3>()]);
#endif
		}


		void
		find_data_index(const float4v& x, const float data[], int N, float4v& value, int4v& ix) {
			int4v found_ix;
			float4v found_value(data[0]);
      mask4v found_mask = mask4v::falses();
#ifdef PERF_TEST
			for (int i=1; i<1; ++i) {
#else
			for (int i=1; i<N; ++i) {
#endif
				mask4v greater = x >= data[i];
        mask4v update = greater && !found_mask;
				found_value = select(update, data[i], found_value);
				found_ix    = select(update, i      , found_ix);
        found_mask |= update;
				// break if less.all_zero()
			}

			value=found_value;
			ix = found_ix;
		}

#if 1
		void
		find_data_index2(const float4v& x, const float data[], int, float4v& value, int4v& ix) {
			const float powixcof[] = { 16.66, -43.245, 21.5 }; //1.7298e-2, 3.44e-6}
			// coeffs valid for w = 2x-1
			float4v float_ix = eval_polynomial(x+x-1.0, powixcof);
			int4v   ix0    = trunc_int(float_ix);  // 0..16 since we truncate, no need to check >=0
			float4v value0 = load_index(data, ix0);
			float4v value1 = load_index(data, ix0+1);

			mask4v less = x <= value1;
			value = select(less, value1, value0);
			ix    = select(less, ix0+1, ix0);
		}
#endif
	}


  template <typename T>
	T
	pow_int(const T& xin, int N) {
		T result = T::ones();
		T x = xin;
		while (N != 0) {
			if ((N & 1) != 0) {
				result *= x;
			}
			x *= x;
			N >>= 1;
		}
    
    return result;
	}


  template <typename T>
	T
	pow_int(const T& xin, const int4v& N) {
    T ones(1);
		T result(1);
		T x = xin;
    int4v n = N;
		while ((n > 0).any_true()) {
			result *= select((n & 1) == 0, ones, x);
			x *= x;
			n >>= 1;
		}

		return result;
	}


  float4v
  pow(const float4v& xin, int N) {
    float4v y = pow_int(xin, std::abs(N));
    return N < 0 ? 1.0f/y : y;
  }


  float4v
  pow(const float4v& xin, const int4v& N) {
    float4v y = pow_int(xin, abs(N));
		return select(N < 0, 1.0f/y, y);
  }


  int4v
  pow(const int4v& xin, int N) {
    int4v y = pow_int(xin, std::abs(N));
    y = select(abs(xin) < 2   ||  int4v(N) > -2, y);
    return y;
    //return select(abs(xin) > 1   &&   N < -1, 0, y ? 1/y : y;
  }


  int4v
  pow(const int4v& xin, const int4v& N) {
    int4v y = pow_int(xin, abs(N));
    y = select(abs(xin) < 2   ||  N > -2, y);
    return y;
		//return select(N < 0, 1/y, y);
  }



  /*! Calculate x**y-1, special treatment to reduce cancellation if x**y is near 1
  */
  float4v
  powm1(const float4v& xin, const float4v& yin) {
    return pow(xin, yin) - 1.0f;
  }

	float4v 
	pow(const float4v& xin, const float4v& yin) {

		int4v nflg;

		float4v x = xin;
		float4v y = yin;
		float4v w = floor(y);

		nflg = 0;	/* flag = 1 if x<0 raised to integer power */
		float4v z = abs(w);

		mask4v intpowermask = w == y;
		if ((intpowermask  &&  z < 32768.0).all_true()) {
			return pow(xin, round_int(y));
		}

		mask4v x_neg = x < 0.0f;
		mask4v fracpowermask= w != y;
		mask4v invalidmask = x_neg && fracpowermask;
		x = abs(x);
		//if( x <= 0.0F )
		//	{
		//	if( x == 0.0 )
		//		{
		//		if( y == 0.0 )
		//			return( 1.0 );  /*   0**0   */
		//		else  
		//			return( 0.0 );  /*   0**y   */
		//		}
		//	else
		//		{
		//		if( w != y )
		//			{ /* noninteger power of negative number */
		//			mtherr( fname, DOMAIN );
		//			return(0.0);
		//			}
		//		nflg = 1;
		//		if( x < 0 )
		//			x = -x;
		//		}
		//	}

		/* separate significand from exponent */
		int4v expnt;
		x = frexp(x, expnt);

		/* find significand in antilog table A[] */
		//i = 1;
		//if( x <= A[9] )
		//	i = 9;
		//if( x <= A[i+4] )
		//	i += 4;
		//if( x <= A[i+2] )
		//	i += 2;
		//if( x >= A[1] )
		//	i = -1;
		//i += 1;

		int4v index;
		float4v a;
		find_data_index(x, A, sizeof(A)/sizeof(A[0]), a, index);
		x -= a;
		x -= load_index(B, index>>1);
		x *= load_index(Ainv, index);

		/* Find (x - A[i])/A[i]
		 * in order to compute log(x/A[i]):
		 *
		 * log(x) = log( a x/a ) = log(a) + log(x/a)
		 *
		 * log(x/a) = log(1+v),  v = x/a - 1 = (x-a)/a
		 */
		//x -= A[i];
		//x -= B[ i >> 1 ];
		//x *= Ainv[i];


		/* rational approximation for log(1+v):
		 *
		 * log(1+v)  =  v  -  0.5 v^2  +  v^3 P(v)
		 * Theoretical relative error of the approximation is 3.5e-11
		 * on the interval 2^(1/16) - 1  > v > 2^(-1/16) - 1
		 */
		{
			float4v xx = x*x;
			const float4v powwcof[] = {
				+ 0.3333331095506474,
				- 0.2500006373383951,
				+ 0.2003770364206271,
				- 0.1663883081054895
			};
			w = eval_polynomial(x, powwcof)*x*xx - 0.5*xx;
		}
		//w = (((-0.1663883081054895  * x
		//			+ 0.2003770364206271) * x
		//			- 0.2500006373383951) * x
		//			+ 0.3333331095506474) * x * z;
		//w -= 0.5 * z;

		/* Convert to base 2 logarithm:
		 * multiply by log2(e)
		 */
		w = w + LOG2EA * w;
		/* Note x was not yet added in
		 * to above rational approximation,
		 * so do it now, while multiplying
		 * by log2(e).
		 */
		z = w + LOG2EA * x;
		z = z + x;

		/* Compute exponent term of the base 2 logarithm. */
		w = -index;
		w *= 0.0625;  /* divide by 16 */
		w += float4v(expnt);
		/* Now base 2 log of x is w + z. */

		/* Multiply base 2 log by y, in extended precision. */

		/* separate y into large part ya
		 * and small part yb less than 1/16
		 */
		{
			float4v ya = reduce(y);
			float4v yb = y - ya;
		
			float4v F = z * y  +  w * yb;
			float4v Fa = reduce(F);
			float4v Fb = F - Fa;

			float4v G = Fa + w * ya;
			float4v Ga = reduce(G);
			float4v Gb = G - Ga;

			float4v H = Fb + Gb;
			float4v Ha = reduce(H);
			w = 16 * (Ga + Ha);

			//F = z * y  +  w * yb;
			//Fa = reduce(F);
			//Fb = F - Fa;

			//G = Fa + w * ya;
			//Ga = reduce(G);
			//Gb = G - Ga;

			//H = Fb + Gb;
			//Ha = reduce(H);
			//w = 16 * (Ga + Ha);

			///* Test the power of 2 for overflow */
			//if( w > MEXP )
			//	{
			//	mtherr( fname, OVERFLOW );
			//	return( MAXNUMF );
			//	}

			//if( w < MNEXP )
			//	{
			//	mtherr( fname, UNDERFLOW );
			//	return( 0.0 );
			//	}

			expnt = trunc_int(w);
			float4v Hb = H - Ha;
			//e = w;
			//Hb = H - Ha;

			expnt += select(Hb > 0.0, int4v(1));
			Hb    -= select(Hb > 0.0, float4v(0.0625));
			//if( Hb > 0.0 ) {
			//	e += 1;
			//	Hb -= 0.0625;
			//}

			/* Now the product y * log2(x)  =  Hb + e/16.0.
			 *
			 * Compute base 2 exponential of Hb,
			 * where -0.0625 <= Hb <= 0.
			 * Theoretical relative error of the approximation is 2.8e-12.
			 */
			/*  z  =  2**Hb - 1    */
			// these coefficients are nearly the same as for exp2
			// slight difference since Hb is small (< 1/16)
			const float4v powHbcof[] = {  
				6.931471791490764E-001,
				2.402262883964191E-001,
				5.549356188719141E-002,
				9.416993633606397E-003
			};

			z = eval_polynomial(Hb, powHbcof) * Hb;
		}
		//z = ((( 9.416993633606397E-003 * Hb
		//			+ 5.549356188719141E-002) * Hb
		//			+ 2.402262883964191E-001) * Hb
		//			+ 6.931471791490764E-001) * Hb;

		/* Express e/16 as an integer plus a negative number of 16ths.
		 * Find lookup table entry for the fractional power of 2.
		 */
		int4v i = select(expnt < 0, -(-expnt >> 4), (expnt >> 4) + 1);

		//if( e < 0 )
		//	i = -( -e >> 4 );
		//else
		//	i = (e >> 4) + 1;

		expnt = (i << 4) - expnt;
		// e = (i << 4) - e;

		w = load_index(A, expnt);
		//w = A[e];
		
		z = w + w * z;      /*    2**-e * ( 1 + (2**Hb-1) )    */
		z = ldexp(z, i);  /* multiply by integer power of 2 */

		//z = w + w * z;      /*    2**-e * ( 1 + (2**Hb-1) )    */
		//z = ldexpf( z, i );  /* multiply by integer power of 2 */

		w = floor(0.5*w);
		w = w + w;
		z = changesign(z, w != y && intpowermask);
		//if( nflg )
		//	{
		///* For negative x,
		// * find out if the integer exponent
		// * is odd or even.
		// */
		//	w = 2 * floorf( (float) 0.5 * w );
		//	if( w != y )
		//		z = -z; /* odd exponent */
		//	}

		// special cases
		// x == 0													0
		// x == 0 && y == 0								1
		// x < 0  && y not integer				NaN
		// x NaN or y NaN									NaN
		// 
		z = select(isnan(x) || invalidmask, float4v::NaNs(), z);

		return z;
	}

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////


	float4v 
	exp2_normalized(const float4v& xin) {
		static const float exp2cof[] = {
			1.000000000000000E+000,
			6.931472028550421E-001,
			2.402264791363012E-001,
			5.550332471162809E-002,
			9.618437357674640E-003,
			1.339887440266574E-003,
			1.535336188319500E-004,
		};

		// rational approximation exp2(x) = 1.0 +  xP(x)
		float4v y = eval_polynomial(xin, exp2cof);

		return y;
	}

		
	float4v 
	log2_normalized(const float4v& xin, int4v& expnt) {

		float4v x = xin;

		// expnt-- if x<sqrt(2)/2
		expnt -= select(x < SQRT2_HALF, int4v(1));
		x += select(x < SQRT2_HALF, x);
		x -= 1.0f;

		const float log2cof[] = {
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
		float LOG2EA = 0.44269504088896340735992;
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

		return z;
	}


	// from ieee754
	const int4v exp_bias(126);
	const int exp_pos = 23;
	const mask4v exponent_mask(0xff << exp_pos);
	const mask4v sign_mask(float4v::signmask());
	const mask4v mantissa_mask = ~(exponent_mask | sign_mask);  // remaining bits to mantissa
	const int4v exp_bias_shifted = exp_bias << exp_pos;

	__forceinline float4v 
	ldexp_local(const float4v& xin, const int4v& exp_in) {
		
		// Exponent limit when to call the standard ldexp
		// i.e. handle small exponents guaranteed to not over- or underflow
		// here quickly
		const int		ldexp_limit 			= 64;
		
		float4v y;
		float4v x = xin;
		int4v expnt = exp_in;
		mask4v large = abs(expnt) > ldexp_limit;
		if (large.any_true()) {
			y = ldexp(xin, exp_in);
		}
		else {
			// don't care about NaN, inf or 0.  They are handled elsewhere
			int4v xbits(xin.reinterpret_bits());

			int4v ebits = (expnt + 1) << exp_pos;  // add exponent bias in ieee format, shift to right bits
			y = float4v::create_reinterpreted(xbits + ebits);
		}

		return y;
	}



	__forceinline float4v 
	frexp_local(const float4v& xin, int4v& expnt_out) {

		float4v y;

		//if (isnormal(xin).all_true()) {
			// fastlane if no special handling needed.
			int4v xbits(xin.reinterpret_bits());   // absolute value, as bits
			xbits  = (xbits & (sign_mask | mantissa_mask)) | exp_bias_shifted; // (exp_bias << exp_pos);

			y = float4v::create_reinterpreted(xbits);
			expnt_out = (abs(xin).reinterpret_bits() >> exp_pos) - exp_bias;
#if 0
		}
		else {
			// handle all special cases
			y = frexp(xin, expnt_out);
		}
#endif

		return y;
	}


		__forceinline mask4v
		isnonzerovalue(const float4v& x) {
			return isfinite(x) && x.reinterpret_bits() != 0;
		}



				
		// special cases
		// y == 0  >>> 1
		// x == 0  >>> 0
		// x < 0 AND y not int >>> NaN
		// x < 0 AND y int >>> normal compute, sign = signum(x)^y
		/*
		Special Values
		pow ( ±0, y ) returns ±∞ and raises the divide-by-zero floating-point exception for y < 0 and an odd integer.
		pow ( ±0, y ) returns +∞ and raises the divide-by-zero floating-point exception for y < 0 and not an odd integer.
		pow ( ±0, y ) returns ±0 for y an odd integer > 0.
		pow ( ±0, y ) returns +0 for y > 0 and not an odd integer.
		pow ( -1, ±∞ ) returns 1.
		pow ( 1, y ) returns 1 for any y, even a NaN.
		pow ( x, ±0 ) returns 1 for any x, even a NaN.
		pow ( x, y ) returns a NaN and raises the invalid floating-point exception for finite x < 0 and finite non-integer y.
		pow ( x, -∞ ) returns +∞ for |x| < 1.
		pow ( x, -∞ ) returns +0 for |x| > 1.
		pow ( x, +∞ ) returns +0 for |x| < 1.
		pow ( x, +∞ ) returns +∞ for |x| > 1.
		pow ( -∞, y ) returns -0 for y an odd integer < 0.
		pow ( -∞, y ) returns +0 for y < 0 and not an odd integer.
		pow ( -∞, y ) returns -∞ for y an odd integer > 0.
		pow ( -∞, y ) returns +∞ for y > 0 and not an odd integer.
		pow ( +∞, y ) returns +0 for y < 0.
		pow ( +∞, y ) returns +∞ for y > 0.
		A domain error occurs if x is finite and negative and y is finite and not an integer.
		A domain error can occur if x is 0 and y less than or equal to 0.
		http://www.codecogs.com/reference/computing/c/math.h/pow.php


		pow ( ±0, y )  =  ±∞   y < 0 and an odd integer.									1/x         													A
		pow ( ±0, y )  =  +∞   y < 0 and not an odd integer.							1/abs(x) or inf												D
		pow ( ±0, y )  =  ±0   y an odd integer > 0.											x																			B

		pow ( ±0, y )  =  +0   y > 0 and not an odd integer.							0 or abs(x)														C
		pow ( -1, ±∞ ) =  1.																							1 or abs(x)														C
		pow ( 1, y )   =  1 for any y, even a NaN.												x or 1																B

		pow ( x, ±0 )  =  1 for any x, even a NaN.												1																			E
		pow ( x, y )   =  NaN finite x < 0 and finite non-integer y.			NaN																		F
		pow ( x, -∞ )  =  +∞ for |x| < 1.																	as is, check trunc_int handles it		  0

		pow ( x, -∞ )  =  +0 for |x| > 1.																	--"--																	0
		pow ( x, +∞ )  =  +0 for |x| < 1.																	--"--																	0
		pow ( x, +∞ )  =  +∞ for |x| > 1.																	--"--																	0

		pow ( -∞, y )  =  -0 for y an odd integer < 0.										1/x																		A
		pow ( -∞, y )  =  +0 for y < 0 and not an odd integer.						-1/x or 1/abs(x)											D
		pow ( -∞, y )  =  -∞ for y an odd integer > 0.										x																			B

		pow ( -∞, y )  =  +∞ for y > 0 and not an odd integer.						-x or abs(x)													C
		pow ( +∞, y )  =  +0 for y < 0.																		1/x or 0, == pow( -∞, y )             A
		pow ( +∞, y )  =  +∞ for y > 0.																		x																			B

		pow ( NaN, y ) = NaN for y != 0																		x																			B
		pow ( x, NaN ) = NaN for x != 1																		NaN or y															F
		*/
		void
		calc_specials(const float4v& x, const float4v& y, float4v& result, mask4v& mask) {
			mask4v y_int       = y == trunc(y);
			mask4v y_odd       = y_int && (trunc_int(y) & 0x1) == 0x1;
			mask4v y_frac      = y != trunc(y);
			mask4v y_neg       = y < 0.0f;
			mask4v y_pos       = y > 0.0f;
			mask4v y_zero      = y == 0.0f;
			mask4v y_inf       = isinf(y);
			mask4v y_nan       = isnan(y);

			mask4v x_neg       = x < 0.0f;
			mask4v x_pos       = x > 0.0f;
			mask4v x_zero      = x == 0.0f;
			mask4v x_m1        = x == -1.0f;
			mask4v x_p1        = x == 1.0f;
			mask4v x_minf      = x == -float4v::infinites();
			mask4v x_pinf      = x == float4v::infinites();
			mask4v x_nan       = isnan(x);

			mask4v invert_x   = x_zero && y_neg && y_odd  ||  x_minf && y_neg && y_odd  || x_pinf && y_neg;			// A
			mask4v copy_x     = y_pos && y_odd            ||  x_p1 
											 || x_minf && y_odd && y_pos  ||  x_pinf && y_pos  || x_nan && y != 0.0f;						// B
			mask4v copy_absx  = x_zero && y_pos && !y_odd ||  x_m1 && y_inf    || x_minf && y_pos && !y_odd;		// C
			mask4v invert_absx= x_zero && y_neg && !y_odd ||  x_minf && y_neg && !y_odd;												// D
			mask4v set_ones   = y_zero;																																					// E
			mask4v set_nan    = isnan(y) && !x_p1;						// F

			// Cases A, B, C, D are mutually exclusive
			float4v specials = blend(select(invert_x, 1.0/x), select(copy_x, x), select(copy_absx, abs(x)), select(invert_absx, 1.0/abs(x)));
			specials = select(set_nan, float4v::NaNs(), specials);
			specials = select(set_ones, 1.0, specials);

			mask = invert_x || copy_x || copy_absx || invert_absx || set_ones || set_nan;
			result = specials;
		}




	float4v 
	pow_v2(const float4v& xin, const float4v& yin) {

		float4v x = abs(xin);
		float4v y = yin;

		/* separate significand from exponent */
		int4v expnt;
		x = frexp_local(x, expnt);

		// rational approximation for log(1+v)
		float4v z = log2_normalized(x, expnt);
		float4v w = float4v(expnt);

		/* Multiply base 2 log by y, in extended precision. 
			 separate y into large part ya and small part yb less than 1
		 */
		{
			float4v ya = round(y);
			float4v yb = y - ya;
		
			float4v F = z * y  +  w * yb;
			float4v Fa = round(F);
			float4v Fb = F - Fa;

			float4v G = Fa + w * ya;
			float4v Ga = round(G);
			float4v Gb = G - Ga;

			float4v H = Fb + Gb;
			float4v Ha = round(H);
			w = Ga + Ha;

			expnt = trunc_int(w);
			float4v Hb = H - Ha;
			float4v adjust = round(Hb);
			Hb -= adjust;					// Hb is now in the range -0.5..+0.5
			expnt += trunc_int(adjust);

			/* Now the product y * log2(x)  =  Hb + expnt
			 * Theoretical relative error of the approximation is 2.8e-12 (???)
			 */
			z = exp2_normalized(Hb);  // TODO should be exp2m1_normalized???
		}

		// For powm1, calc z=exp2m1_normalized instead
		// then z = select(expnt == 0, z, ldexp(z+1.0, expnt)
		z = ldexp_local(z, expnt);		/* multiply by integer power of 2 */

		// change sign of result?
		{
			mask4v y_int       = y == trunc(y);
			int4v  a           = trunc_int(y);
			mask4v y_odd       = a != ((a >> 1) << 1);

			// mask4(xin...) is a less intuitive but quicker way to say xin<0.0
			mask4v altersign = mask4v(xin.reinterpret_bits().value()) && y_odd && y_int;
			z = changesign(z, altersign);
		}

		// prepare result if x or y is 0, inf or nan
		{
			float4v special_result;
			mask4v  special_mask;
			if ((isnonzerovalue(x) && isnonzerovalue(y)).any_false()) {
				calc_specials(xin, yin, special_result, special_mask);
			}

			mask4v invalid = xin < 0.0f && isfinite(x) && isfinite(y) && (y != trunc(y));
			special_result = select(invalid, float4v::NaNs(), special_result);
			special_mask |= invalid;
			z = select(special_mask, special_result, z);
		}

		return z;
	}

}}}}