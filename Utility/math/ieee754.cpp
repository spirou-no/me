/*		
 *							frexpf()
 *							ldexpf()
 *							signbitf()  -- inlined
 *							isnanf()    -- inlined
 *							isfinitef() -- inlined
 *							nextafter()
 *							fpclassify()
 *
 *	Single precision floating point numeric utilities
 *
 *
 *
 * SYNOPSIS:
 *
 * y = floorf(x);
 * y = ceilf(x);
 * y = frexpf( x, &expnt );
 * y = ldexpf( x, n );
 * n = signbit(x);
 * n = isnan(x);
 * n = isfinite(x);
 *
 *
 *
 * DESCRIPTION:
 *
 * All four routines return a single precision floating point
 * result.
 *
 * sfloor() returns the largest integer less than or equal to x.
 * It truncates toward minus infinity.
 *
 * sceil() returns the smallest integer greater than or equal
 * to x.  It truncates toward plus infinity.
 *
 * sfrexp() extracts the exponent from x.  It returns an integer
 * power of two to expnt and the significand between 0.5 and 1
 * to y.  Thus  x = y * 2**expn.
 *
 * ldexpf() multiplies x by 2**n.
 *
 * signbit(x) returns 1 if the sign bit of x is 1, else 0.
 *
 * These functions are part of the standard C run time library
 * for many but not all C compilers.  The ones supplied are
 * written in C for either DEC or IEEE arithmetic.  They should
 * be used only if your compiler library does not already have
 * them.
 *
 * The IEEE versions assume that denormal numbers are implemented
 * in the arithmetic.  Some modifications will be required if
 * the arithmetic has abrupt rather than gradual underflow.
 */


/*
Cephes Math Library Release 2.1:  December, 1988
Copyright 1984, 1987, 1988 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <float.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {
	/* For description of single precision float see http://en.wikipedia.org/wiki/Single_precision
	31			sign
	30..23  exponent
	22..0   mantissa

	normal values have exp != 0 && exp != 0xff
	exp		mantissa
	0			0							really zero, +0 or -0 from sign bit
	0			<>0						denormal number
	x			x							normal number, 1.mantissa * 2^(2-126)
	ff		0							+- infinity
	ff		<>0						NaN
	
	The minimum positive (subnormal) value is 2e−149 ≈ 1.4 × 10e−45. 
	The minimum positive normal value is 2e−126 ≈ 1.18 × 10e−38. 
	The maximum representable value is (2−2e−23) × 2e127 ≈ 3.4 × 10e38.
	*/
	const int4v exp_bias(126);
	const int exp_pos = 23;
	const mask4v exponent_mask(0xff << exp_pos);
	const mask4v sign_mask(float4v::signmask());
	const mask4v mantissa_mask = ~(exponent_mask | sign_mask);  // remaining bits to mantissa
	const int4v  inf_or_nan(0xff << exp_pos);
	const int denorm_shift = 25;
	const float pow_2_denorm = float(1<<denorm_shift);




	float4v 
	frexp(const float4v& xin, int4v& expnt_out) {

		float4v x = abs(xin);
		int4v xbits(x.reinterpret_bits());   // absolute value, as bits

		mask4v keepmask = (xbits >= inf_or_nan) || (xbits == 0);
		//mask4v keepmask = ((xbits >> exp_pos) == inf_or_nan) || (xbits == 0);
		int4v  expnt    = int4v::zeros();

		// 1<<exp_pos is the smallest normalized number possible (exp bits==1, mantissa=0)
		// only prevent denorm handling if xin=0.0 (the logic with keepmask)
		mask4v denormalized = (xbits < (1<<exp_pos)) && !keepmask;  
		if (denormalized.any_true()) {
			x *= select(denormalized, float4v(pow_2_denorm), 1.0f);
			xbits = int4v(x.reinterpret_bits());
			expnt = select(denormalized, int4v(-denorm_shift));
		}
    
		expnt += (xbits >> exp_pos) - exp_bias;
		//xbits  = (xbits & 0x007fffff ) | 0x3f000000; 
		// -- was (xbits & 0x807fffff ) | 0x3f000000; but xbits is unsigned, so no point in anding with signmask here
		xbits  = (xbits & mantissa_mask) | (exp_bias << exp_pos);

		float4v y = float4v::create_reinterpreted(xbits);
		y = copysign(y, xin);
		y = select(keepmask, xin, y);
		expnt_out = select(keepmask, 0, expnt);

		return y;
	}


	float4v 
	ldexp(const float4v& xin, const int4v& exp_in) {
		
		// the absolute largest range of expin is going from the smallest denormalized number
		// to the largest float, i.e. 255 + 23 bits precision in mantissa.
		const int		scale_up_pow 			= 127;
		//const float scale_up_factor		= 0.1701411834604692e39;   // 0x1p127
		const int		scale_down_pow 		= -126;
		//const float scale_down_factor = 0.1175494350822288e-37;  // 0x1p-126

		const float4v scale_up   = float4v::create_reinterpreted(mask4v(0xfe << exp_pos));
		const float4v scale_down = float4v::create_reinterpreted(mask4v(0x01 << exp_pos));

		// don't change NaN, infinity or 0.0
		//mask4v keepmask = xbits >= 0x7f800000 || xbits == 0;
		int4v xbits(xin.reinterpret_bits());
		mask4v keepmask = xbits >= inf_or_nan || xbits == 0;
		int4v exp = clip(exp_in, -300, 300);  // ca 279 powers of 2 is the span from smallest denorm to largest float.x

		float4v x = xin;
		mask4v large = exp > scale_up_pow;
		while (large.any_true()) {
			exp -= select(large, int4v(scale_up_pow));
			x *= select(large, scale_up, 1.0f);
			large = exp > scale_up_pow;
		}

		mask4v small = exp < scale_down_pow;
		while (small.any_true()) {
			exp -= select(small, int4v(scale_down_pow));  // recall scale_down_pow is negative
			x *= select(small, scale_down, 1.0f);
			small = exp < scale_down_pow;
		}

		// -126 <= exp <= 127; convert n to float scale factor.                      *
		// 1<<exp_pos is the smallest normalized number possible (exp bits==1, mantissa=0)
		// only prevent denorm handling if xin=0.0 (the logic with keepmask)

		int4v ebits = (exp + exp_bias + 1) << exp_pos;  // add exponent bias in ieee format, shift to right bits
      
		float4v exp_scale = float4v::create_reinterpreted(ebits);
		float4v y = x * exp_scale;
    
		y = select(keepmask, xin, y);
		
		return y;
	}




	//! Calculate x with n bits of precision truncated.
	// Used to calculate constants in extended precision argument reduction
	float4v truncate_prec(const float4v& xin, int trunc_prec) {
		int4v xbits(xin.reinterpret_bits());
		int4v mantissa = xbits & int4v(mantissa_mask);
		mantissa = ((mantissa + 1) >> trunc_prec) << trunc_prec;

		mask4v overflow = (mantissa & int4v(mantissa_mask)) != 0;
		int4v  addexp   = select(overflow, int4v::ones());
		int4v  exponent = xbits & exponent_mask;
		exponent += addexp;

		overflow    = (exponent & exponent_mask) != 0;
		int4v ybits = mantissa | exponent;
		float4v   y = copysign(float4v::create_reinterpreted(ybits), xin);

		return y;
	}

#if 0
	float4v nextafter(const float4v& x, const float4v& y) {
		float fp0 = _nextafter(x.at<0>(), y.at<0>());
		float fp1 = _nextafter(x.at<1>(), y.at<0>());
		float fp2 = _nextafter(x.at<2>(), y.at<0>());
		float fp3 = _nextafter(x.at<3>(), y.at<0>());

		return float4v(fp0, fp1, fp2, fp3);
	}
#endif
	int4v fpclass(const float4v& x) {
		int fp0 = _fpclass(x.at<0>()); 
		int fp1 = _fpclass(x.at<1>());
		int fp2 = _fpclass(x.at<2>());
		int fp3 = _fpclass(x.at<3>());

		return int4v(fp0, fp1, fp2, fp3);
	}

}}}}