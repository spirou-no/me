/*							gammaf.c
 *
 *	Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, gammaf();
 * extern int sgngamf;
 *
 * y = gammaf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns gamma function of the argument.  The result is
 * correctly signed, and the sign (+1 or -1) is also
 * returned in a global (extern) variable named sgngamf.
 * This same variable is also filled in by the logarithmic
 * gamma function lgam().
 *
 * Arguments between 0 and 10 are reduced by recurrence and the
 * function is approximated by a polynomial function covering
 * the interval (2,3).  Large arguments are handled by Stirling's
 * formula. Negative arguments are made positive using
 * a reflection formula.  
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,-33      100,000     5.7e-7      1.0e-7
 *    IEEE       -33,0      100,000     6.1e-7      1.2e-7
 *
 *
 */
/*							lgamf()
 *
 *	Natural logarithm of gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, lgamf();
 * extern int sgngamf;
 *
 * y = lgamf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the absolute
 * value of the gamma function of the argument.
 * The sign (+1 or -1) of the gamma function is returned in a
 * global (extern) variable named sgngamf.
 *
 * For arguments greater than 6.5, the logarithm of the gamma
 * function is approximated by the logarithmic version of
 * Stirling's formula.  Arguments between 0 and +6.5 are reduced by
 * by recurrence to the interval [.75,1.25] or [1.5,2.5] of a rational
 * approximation.  The cosecant reflection formula is employed for
 * arguments less than zero.
 *
 * Arguments greater than MAXLGM = 2.035093e36 return MAXNUM and an
 * error message.
 *
 *
 *
 * ACCURACY:
 *
 *
 *
 * arithmetic      domain        # trials     peak         rms
 *    IEEE        -100,+100       500,000    7.4e-7       6.8e-8
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 * The routine has low relative error for positive arguments.
 *
 * The following test used the relative error criterion.
 *    IEEE    -2, +3              100000     4.0e-7      5.6e-8
 *
 */

/*							gamma.c	*/
/*	gamma function	*/

/*
Cephes Math Library Release 2.7:  July, 1998
Copyright 1984, 1987, 1989, 1992, 1998 by Stephen L. Moshier
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>
#include <math.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	/* define MAXGAM 34.84425627277176174 */

	/* Stirling's formula for the gamma function
	 * gamma(x) = sqrt(2 pi) x^(x-.5) exp(-x) ( 1 + 1/x P(1/x) )
	 * .028 < 1/x < .1
	 * relative error < 1.9e-11
	 */
	const float MAXSTIR = 26.77;
	const float SQTPI   = 2.50662827463100050242; /* sqrt( 2 pi ) */

	int sgngamf = 0;
	extern int sgngamf;

	/* Gamma function computed by Stirling's formula,
	 * sqrt(2 pi) x^(x-.5) exp(-x) (1 + 1/x P(1/x))
	 * The polynomial STIR is valid for 33 <= x <= 172.
	 */
	float4v
	stirling(const float4v& xin) {
		const float stircof[] = {
			 8.333331788340907E-002,
			 3.473255786154910E-003,
			-2.705194986674176E-003,
		};

		float4v x = xin;
		float4v w = 1.0/x;
		w = 1.0 + w * eval_polynomial(w, stircof);
		float4v y = exp(-x);
		
		// multiply exponent in pow by 0.5 to avoid overflow, square result after.
		float4v v = pow(x, 0.5*(x-0.5));
		
		y *= v;
		y *= v;  
		y *= SQTPI * w;

		return y;
	}


	/* gamma(x+2), 0 < x < 1 */
	static float P[] = {
	 1.536830450601906E-003,
	 5.397581592950993E-003,
	 4.130370201859976E-003,
	 7.232307985516519E-002,
	 8.203960091619193E-002,
	 4.117857447645796E-001,
	 4.227867745131584E-001,
	 9.999999822945073E-001,
	};


#if 0
	// TODO C++20
	float4v 
	gamma(const float4v& xin) {
		
		float4v x = xin;
		mask4v   signgam;
		mask4v  negative = x < 0.0;
		float4v nz			 = 0.0;
		mask4v  infmask;

		if (negative.any_true())	{
			float4v q = -x;
			float4v p = floor(q);
			infmask = p == q;

			int4v i = trunc_int(p);
			signgam = (i & 0x1) == 0;
			nz = q - p;

			p += round(nz);
			nz = q - p;
			
			nz = q * sin(PI * nz);
			infmask |= nz == 0.0;
			
			nz = abs(nz);
			x = q;
		}
		

		mask4v hugearg = x >= 10.0;
		float4v zhuge;
		float4v zsmall;

		if(hugearg.any_true())	{
			zhuge = stirling(x);
		}
		
		if ((x < 10.0).any_true()) {
			mask4v direction = x < 2.0;
		
			float4v z = 1.0;
			mask4v geq3 = x >= 3.0; 
			while(geq3.any_true())	{
				float4v tx = x - 1.0;
				float4v tz = z * tx;
				z = select(geq3, tz, z);
				x = select(geq3, tx, x);
				geq3 = x >= 3.0;
			}
		
			mask4v nottiny = x >= 1.e-4;
			mask4v less2 = x < 2.0;

			while (less2.any_true()) {
				float4v tz = z * x;
				float4v tx = x + 1.0;
				z = select(less2 && nottiny, tz, z);
				x = select(less2 && nottiny, tx, x);
				less2 = x < 2.0;
			}

			mask4v tinyarg = x < 1e-4; 
			float4v ptiny = z / ((1.0 + EULER * x) * x); EULER
			infmask &= x == 0.0;

			z = select(direction, 1.0/z, z);

			mask4v  arg2          = x == 2.0;
			float4v ztwo          = select(arg2, z);
			mask4v  inf_or_2_mask = infmask || arg2;
			float4v inf_or_2      = blend(ztwo, float4v::infinites());

			x -= 2.0f;
			p = select(tinyarg, ptiny, z * eval_polynomial(x, P));

			float4v pneg = changesign(PI/(nz * p), signgam);
			
			p = select(negative, pneg, p);
			zsmall = select(inf_or_2_mask, inf_or_2, p);
		}

		return select(hugearg, zhuge, zsmall);
	}
#endif



	/* log gamma(x+2), -.5 < x < .5 */
	const float B[] = {
	 6.055172732649237E-004,
	-1.311620815545743E-003,
	 2.863437556468661E-003,
	-7.366775108654962E-003,
	 2.058355474821512E-002,
	-6.735323259371034E-002,
	 3.224669577325661E-001,
	 4.227843421859038E-001
	};

	/* log gamma(x+1), -.25 < x < .25 */
	const float C[] = {
	 1.369488127325832E-001,
	-1.590086327657347E-001,
	 1.692415923504637E-001,
	-2.067882815621965E-001,
	 2.705806208275915E-001,
	-4.006931650563372E-001,
	 8.224670749082976E-001,
	-5.772156501719101E-001
	};

	/* log( sqrt( 2*pi ) ) */
	static float LS2PI  =  0.91893853320467274178;
	#define MAXLGM 2.035093e36
	static float PIINV =  0.318309886183790671538;

	/* Logarithm of gamma function */


#if 0
	// TODO C++20
	float4v 
	lgam(const float4v& xin) {
		float p, q, w, z, x;
		float nx, tx;
		int i, direction; 

		sgngamf = 1;

		x = xx;
		if( x < 0.0 ) {
			q = -x;
			w = lgam(q); /* note this modifies sgngam! */
			p = floorf(q);
			if( p == q )
				goto loverf;
			i = p;
			if( (i & 1) == 0 )
				sgngamf = -1;
			else
				sgngamf = 1;
			z = q - p;
			if( z > 0.5 ) {
				p += 1.0;
				z = p - q;
			}
			
			z = q * sinf( PIF * z );
			if( z == 0.0 )
				goto loverf;
			z = -logf( PIINV*z ) - w;
			return( z );
		}

		if( x < 6.5 ) {
			direction = 0;
			z = 1.0;
			tx = x;
			nx = 0.0;
			if( x >= 1.5 ) {
				while( tx > 2.5 ) {
					nx -= 1.0;
					tx = x + nx;
					z *=tx;
				}
				x += nx - 2.0;
		iv1r5:
				p = x * polevlf( x, B, 7 );
				goto cont;
			}
			if( x >= 1.25 )	{
				z *= x;
				x -= 1.0; /* x + 1 - 2 */
				direction = 1;
				goto iv1r5;
			}
			if( x >= 0.75 )	{
				x -= 1.0;
				p = x * polevlf( x, C, 7 );
				q = 0.0;
				goto contz;
			}
			while( tx < 1.5 ) {
				if( tx == 0.0 )
					goto loverf;
				z *=tx;
				nx += 1.0;
				tx = x + nx;
			}
			direction = 1;
			x += nx - 2.0;
			p = x * polevlf( x, B, 7 );

		cont:
			if( z < 0.0 ) {
				sgngamf = -1;
				z = -z;
			}
			else {
				sgngamf = 1;
			}
			q = logf(z);
			if( direction )
				q = -q;
		contz:
			return( p + q );
		}

		if( x > MAXLGM ) {
		loverf:
			mtherr( "lgamf", OVERFLOW );
			return( sgngamf * MAXNUMF );
		}

		/* Note, though an asymptotic formula could be used for x >= 3,
		 * there is cancellation error in the following if x < 6.5.  */
		q = LS2PI - x;
		q += ( x - 0.5 ) * logf(x);

		if( x <= 1.0e4 ) {
			z = 1.0/x;
			p = z * z;
			q += ((    6.789774945028216E-004 * p
				 - 2.769887652139868E-003 ) * p
				+  8.333316229807355E-002 ) * z;
		}

		return( q );
	}
#endif
}}}}