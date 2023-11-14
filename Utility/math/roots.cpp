/*							cbrtf.c
 *
 *	Cube root
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, cbrtf();
 *
 * y = cbrtf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the cube root of the argument, which may be negative.
 *
 * Range reduction involves determining the power of 2 of
 * the argument.  A polynomial of degree 2 applied to the
 * mantissa, and multiplication by the cube root of 1, 2, or 4
 * approximates the root to within about 0.1%.  Then Newton's
 * iteration is used to converge to an accurate result.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,1e38      100000      7.6e-8      2.7e-8
 *
 */
/*							cbrt.c  */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1988, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {


	/* Cephes is quite slow, 100+ns
	Try this instead http://metamerist.com/cbrt/cbrt.htm
	*/
	float4v cbrt(const float4v& xin) {
		
		const float CBRT2 = 1.25992104989487316477;
		const float CBRT4 = 1.58740105196819947475;
		const float INV_CBRT2 = float(1.0/CBRT2);
		const float INV_CBRT4 = float(1.0/CBRT4);

		/* extract power of 2, leaving
		 * mantissa between 0.5 and 1
		 */
		int4v expnt;
		float4v x = frexp(xin, expnt);

		/* Approximate cube root of number between .5 and 1,
		 * peak relative error = 9.2e-6
		 */
		const float cbrtcof[] = {
			  0.40238979564544752126924,
			  1.1399983354717293273738,
			 -0.95438224771509446525043,
			  0.54664601366395524503440,
			 -0.13466110473359520655053
		};

		float4v y = eval_polynomial(x, cbrtcof);

		// really integer divide by 3.  But no int divide.
		// Use float div instead, add 0.5 to avoid result is .99999
		int4v exp_div3 = trunc_int((float4v(abs(expnt))+0.5f) / 3.0f);
		int4v rem = abs(expnt) - 3*exp_div3; 

		float4v part_p1 = select(expnt >= 0 && rem == 1, float4v(CBRT2));
		float4v part_p2 = select(expnt >= 0 && rem == 2, float4v(CBRT4));
		float4v part_m1 = select(expnt <  0 && rem == 1, float4v(INV_CBRT2));
		float4v part_m2 = select(expnt <  0 && rem == 2, float4v(INV_CBRT4));
		float4v part_00 = select(rem == 0, float4v(1.0));

		float4v factor = blend(part_00, blend(part_p1, part_p2, part_m1, part_m2));
		expnt = copysign(exp_div3, expnt);

		y *= factor;

		/* multiply by power of 2 */
		y = ldexp(y, expnt);

		/* Newton iteration */
		y -= ( y - (abs(xin)/(y*y)) ) * 0.333333333333;

		y = copysign(y, xin);

		y = select(isfinite(xin), y, xin);

		return y;
	}



/*							cbrt_halleys.c
 *
 *	Cube root
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, cbrtf();
 *
 * y = cbrtf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the cube root of the argument, which may be negative.
 *
 * The method is based in Kahans/Turkowski which is based on the fact the bits of
 * an ieee float is very similar to log2(x)*2^23 when reinterpreted as an integer.
 * This similarity is used to quickly calculate an initial seed (5-6 bit precision).
 * The seed is then refined twice using Halleys method (converges faster than Newton-Raphson).
 * 
 * See http://metamerist.com/cbrt/cbrt.htm
 */

	namespace {
		// cube root approximation using bit hack for 32-bit float
		// http://metamerist.com/cbrt/CubeRoot.cpp
		__forceinline float4v cbrt_5(const float4v& x) {
			int4v xbits = int4v(x.reinterpret_bits());
			xbits = trunc_int(float4v(xbits)*0.333333333333f) + 709921077;

			float4v y = float4v::create_reinterpreted(xbits); 
			return y;
		}


		__forceinline float4v nth_root_5(const float4v& x, int n) {
			const int ebits = 8;
			const int fbits = 23;
			const int bias = (1 << (ebits-1))-1;

			int4v xbits(x.reinterpret_bits());
			xbits = trunc_int(float4v(xbits - (bias << fbits)) / float4v(float(n))) + (bias << fbits);

			float4v y = float4v::create_reinterpreted(xbits);
			return y;
		}


		// iterative cube root approximation using Halley's method (float)
		__forceinline float4v cbrta_halley(const float4v& a, const float4v& R) {
			const float4v a3 = a*a*a;
			const float4v b = a * ((a3 + R + R) / (a3 + a3 + R));
			return b;
		}
	}


	float4v cbrt_halleys(const float4v& xin) {
#if 1
		// in its direct form, Halleys method converges slowly for very small x (1e-30).
		// to avoid problems, scale input by 2^3n, and divide result by 2^n
		const int exp_pos = 23;
		const int exp_bias = 126;
		const int scale_up = 96;							// multiply small args by 2^96
		const int scale_down = scale_up/3;		// then divide result by 2^(96/3)

		const float4v scale = float4v::create_reinterpreted((exp_bias + scale_up + 1) << exp_pos);
		const float4v inv_scale = float4v::create_reinterpreted((exp_bias - scale_down + 1) << exp_pos);

		float4v x = abs(xin);
		mask4v smallarg = x < 1e-18;
		float4v x_scaled = x * select(smallarg, scale, 1.0f);
		
		float4v y_scaled = cbrt_5(x_scaled);
		y_scaled = cbrta_halley(y_scaled, x_scaled);
		y_scaled = cbrta_halley(y_scaled, x_scaled);
		// x = cbrta_halley(x, xin);

		float4v y = y_scaled * select(smallarg, inv_scale, 1.0f);
		y = copysign(y, xin);

		// only keep result if input is finite, otherwise copy NaN/infinite
		y = select(isfinite(xin), y, xin);  

		return y;
#else
		float4v x = abs(xin);
		
		float4v y = cbrt_5(x);
		y = cbrta_halley(y, x);
		y = cbrta_halley(y, x);
		// x = cbrta_halley(x, xin);

		y = copysign(y, xin);

		// only keep result if input is finite, otherwise copy NaN/infinite
		y = select(isfinite(xin), y, xin);  

		return y;
#endif
	}

	/*! Calculate y=sqrt(x) - 1, calculated with high precision for x close to 1
	*/
	float4v
	sqrtm1(const float4v& xin) {
		return sqrt(xin) - 1.0;  // todo add series expansion to avoid cancellation
	}



    /*! Find roots of quadratic equation ax*x + b*x + c.
    The return roots are unordered.
    For complex roots, the actual roots are r0+i*r1 and r0-i*r1
    */
    void
	quadratic_roots(const float4v& a, const float4v& b, const float4v& c,
					float4v& r0,        //!< root 1 if real root, real part if complex root
                    float4v& r1,        //!< root 2 if real root, imaginary part if complex root
                    mask4v&  realroots  //!< 0xff.. for real roots, 0x00... for complex roots
    ) {

        float4v t2 = b*b - 4.0f*a*c;
        float4v t = sqrt(t2);
        
        realroots = t2 >= 0.0f;
        float4v r_large = select(realroots, -b-copysign(t, b), -b);
        r_large /= a+a;
        
        float4v r_small = select(realroots, r_large*c/a, t/(a+a));

        r0 = r_large;
        r1 = r_small;
    }
  
    
    /*! Find roots of quadratic equation ax*x + b*x + c, return only real part
    */
	void
	quadratic_real_roots(const float4v& a, const float4v& b, const float4v& c,
						 float4v& r0, float4v& r1) {
        float4v t2 = b*b - 4.0f*a*c;
        float4v t = sqrt(t2);
        mask4v realroots = t2 >= 0.0f;
        float4v r_large = select(realroots, -b-copysign(t, b), -b);
        r_large /= a+a;
        float4v r_small = select(realroots, r_large*c/a, r_large);

        r0 = r_large;
        r1 = r_small;
    }


}}}}
