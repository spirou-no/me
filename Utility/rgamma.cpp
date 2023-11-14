/*						rgammaf.c
 *
 *	Reciprocal gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, rgammaf();
 *
 * y = rgammaf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns one divided by the gamma function of the argument.
 *
 * The function is approximated by a Chebyshev expansion in
 * the interval [0,1].  Range reduction is by recurrence
 * for arguments between -34.034 and +34.84425627277176174.
 * 1/MAXNUMF is returned for positive arguments outside this
 * range.
 *
 * The reciprocal gamma function has no singularities,
 * but overflow and underflow may occur for large arguments.
 * These conditions return either MAXNUMF or 1/MAXNUMF with
 * appropriate sign.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE     -34,+34      100000      8.9e-7      1.1e-7
 */

/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1985, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    /* Chebyshev coefficients for reciprocal gamma function
     * in interval 0 to 1.  Function is 1/(x gamma(x)) - 1
     */

    const float R[] = {
         1.27546015610523951063E-1,
        -4.98558728684003594785E-3,
        -6.41925436109158228810E-2,
         5.06579864028608725080E-3,
         4.16609138709688864714E-4,
        -8.04814124978471142852E-5,
         2.96001177518801696639E-6,
         2.68975996440595483619E-7,
        -3.33964630686836942556E-8,
         1.08965386454418662084E-9
    };


    float4v rgamma(const float4v& x) {
        mask4v nanarg = isnan(x);
        float4v y = select(nanarg, float4v::NaNs());
        mask4v openargs = !nanarg;

        mask4v bigarg = (x > 34.84425627277176174f) && openargs;
        y = select(bigarg, 1.0f/float4v::infinites(), y);  // i.e. +0
        openargs = openargs && !bigarg;

        mask4v smallarg = (x < -34.034f) && openargs;
        if (smallarg.any_true()) {
	        float4v w = -x;
	        float4v z = sin(PI*w);
            float4v v = abs(z);
            mask4v is_zero = (z == 0.0f) && openargs;
            y = select(is_zero, 0.0f, y);
            v = select(is_zero, 1.0f, v);  // avoid strange paths for unused args
	        float4v yt = log(w * v / PI) + lgam(w);
            mask4v underflow = (yt < -MAXLOG) && openargs;
            y = select(underflow, copysign(1.0f/float4v::infinites(), z), y);
	    
            mask4v overflow = (yt > MAXLOG) && openargs;
            y = select(overflow, copysign(float4v::infinites(), z), y);
            openargs = openargs && !(underflow || overflow);
	        openargs = smallarg && openargs;
            y = select(openargs, copysign(exp(y), z), y);
        }

        if (openargs.any_true()) {
            float4v z = 1.0f;
            float4v w = x;

            mask4v reiterate = w > 1.0f;
            while ((reiterate && openargs).any_true()) {
                w -= select(reiterate, float4v::ones());  // 
                z = select(reiterate, z*w, z);
                reiterate = w > 1.0f;
            }

            reiterate = w < 0.0f;
            while ((reiterate && openargs).any_true()) {
                z = select(reiterate, z/w, z);
                w += select(reiterate, float4v::ones());  // 
                reiterate = w < 0.0f;
            }

            while ((w<0.0f).any_true()) {
                z /= w;
                w += float4v::ones();
            }

            mask4v is_zero = w == 0.0f  &&  openargs;
            y = select(is_zero, 0.0f, y);

            mask4v is_one = w == 1.0f  &&  openargs;
            y = select(is_one, 1.0f/z, y);

            openargs = openargs && !(is_zero || is_one);

            float4v yn = w * (1.0f + eval_chebychev(4.0f*w-2.0f, R)) / z;
            y = select(openargs, yn, y);
        }

        return y;
    }
 }}}}
