/*							atanf.c
 *
 *	Inverse circular tangent
 *      (arctangent)
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, atanf();
 *
 * y = atanf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle between -pi/2 and +pi/2 whose tangent
 * is x.
 *
 * Range reduction is from four intervals into the interval
 * from zero to  tan( pi/8 ).  A polynomial approximates
 * the function in this basic interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10, 10     100000      1.9e-7      4.1e-8
 *
 */
/*							atan2f()
 *
 *	Quadrant correct inverse circular tangent
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, z, atan2f();
 *
 * z = atan2f( y, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns radian angle whose tangent is y/x.
 * Define compile time symbol ANSIC = 1 for ANSI standard,
 * range -PI < z <= +PI, args (y,x); else ANSIC = 0 for range
 * 0 to 2PI, args (x,y).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      -10, 10     100000      1.9e-7      4.1e-8
 * See atan.c.
 *
 */

/*							atan.c */


/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

/* Single precision circular arcsine
 * test interval: [-tan(pi/8), +tan(pi/8)]
 * trials: 10000
 * peak relative error: 7.7e-8
 * rms relative error: 2.9e-8
 */

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

	float4v atan(const float4v& xin) {
	
		float4v x = abs(xin);

		mask4v bigarg = x > 2.414213562373095;				// tan 3pi/8
		mask4v medarg = x > 0.4142135623730950;				// tan pi/8 -- real4v(ly x>0.414.. AND !bigarg

		float4v y = select(bigarg, PI_HALF, select(medarg, float4v(PI_QUART))); 
		x = select(bigarg, -1.0/x, select(medarg, (x-1.0)/(x+1.0), x));

		const float atancof[] = {
			-3.33329491539E-1,
			 1.99777106478E-1,
			-1.38776856032E-1,
			 8.05374449538e-2
		};

		float4v z = x*x;
		y += eval_polynomial(z, atancof) * z*x + x;
		y = copysign(y, xin);

		return y;
	}

	
	
	
	float4v atan2(const float4v& yin, const float4v& xin) {
	
		mask4v x_neg = xin < 0.0;
		mask4v x_zero= xin == 0.0;
		mask4v y_neg = yin < 0.0;
		mask4v y_zero= yin == 0.0;

		/* Special cases
				x   y   
				0  >0     pi/2
				0  <0    -pi/2
				0   0     0
			 >0   0     0
			 <0   0     pi
		*/
		float4v z_x0 = select(x_zero, float4v(PI_HALF));
		z_x0 = changesign(z_x0, y_neg);   // also pi/2 if y==0
		z_x0 = select(y_zero, float4v::zeros(), z_x0);  // set to 0.0 if y=0

		float4v z_y0 = select(y_zero, select(x_neg, float4v(PI)));
		float4v z_special = blend(z_x0, z_y0);

		float4v w = select(x_neg, float4v(PI));
		w = changesign(w, y_neg);

		// if xin is 0.0, result will be NaN or infinite.
		// will ignore that result anyway, and replace with z_special
		float4v z = atan(yin/xin);
		z += w;
		z = select(x_zero || y_zero, z_special, z);

		return z;

		/*
		float4v z_x0 = if_then(x_zero, float4v(PI_HALF));
		z_x0 = changesign(z_x0, y_neg);   // also pi/2 if y==0
		z_x0 = if_then_else(y_zero, float4v::zeros(), z_x0);  // set to 0.0 if y=0

		float4v z_y0 = if_then(y_zero, if_then_else(x_neg, float4v(PI)));
		float4v z_special = blend(z_x0, z_y0);

		float4v w = if_then(x_neg, float4v(PI));
		*/
	}

}}}}
