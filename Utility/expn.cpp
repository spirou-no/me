/*							expnf.c
 *
 *		Exponential integral En
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * float x, y, expnf();
 *
 * y = expnf( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the exponential integral
 *
 *                 inf.
 *                   -
 *                  | |   -xt
 *                  |    e
 *      E (x)  =    |    ----  dt.
 *       n          |      n
 *                | |     t
 *                 -
 *                  1
 *
 *
 * Both n and x must be nonnegative.
 *
 * The routine employs either a power series, a continued
 * fraction, or an asymptotic formula depending on the
 * relative values of n and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0, 30       10000       5.6e-7      1.2e-7
 *
 */

/*							expn.c	*/

/* Cephes Math Library Release 2.2:  July, 1992
 * Copyright 1985, 1992 by Stephen L. Moshier
 * Direct inquiries to 30 Frost Street, Cambridge, MA 02140 */

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

    const float BIG = 16777216.;

    
    float4v expn(int n, const float4v& x) {
        
        mask4v nanarg = isnan(x);
        float4v y = select(nanarg, x);
        mask4v unprocessed = !nanarg;

        mask4v domainerr = (int4v(n) < 0  || x < 0.0f || (x == 0.0f  && int4v(n) < 2))  && unprocessed;
        y = select(domainerr, float4v::infinites(), y);

        unprocessed = unprocessed && !domainerr;
        mask4v bigarg = (x > MAXLOG)  && unprocessed;
        y = select(bigarg, float4v::zeros(), y);
        unprocessed = unprocessed && !bigarg;

        mask4v zeroarg = (x == 0.0f) && unprocessed;
        y = select(zeroarg, float4v::ones()/(n-1.0f), y);
        unprocessed = unprocessed && !zeroarg;

        if (n == 0) {
	        float4v t = exp(-x)/x;
            y = select(unprocessed, t, y);
            return y;
        }

        /*							expn.c	*/
        /*		Expansion for large n		*/
        if (n > 5000) {
            float4v xk = x + float(n);
            float4v yk = 1.0f/(xk * xk);
            float4v nf = float(n);

            float4v t = yk * nf * (6.0f* x*x - 8.0f*t*x + t*t);
            t = yk * (t + nf * (nf - 2.0f*x));
            t = yk * (t + nf);
            t = (t + 1.0f) * exp(-x)/xk;

            y = select(unprocessed, t, y);
            return y;
        }


        mask4v greaterone = x > 1.0f && unprocessed;
        if (greaterone.any_true()) {
            /*		continued fraction		*/
            int k = 1;
            float4v pkm2 = 1.0;
            float4v qkm2 = x;
            float4v pkm1 = 1.0;
            float4v qkm1 = x + float(n);
            float4v t = pkm1/qkm1;
            float4v eps;

            do {
                k += 1;
                mask4v is_even = int4v(k & 1) == 0;
                float4v yk = select(is_even, x, 1.0f);
                float4v xk = select(is_even, float4v(float(k/2)), float(n+(k-1)/2));

	            float4v pk = pkm1 * yk  +  pkm2 * xk;
	            float4v qk = qkm1 * yk  +  qkm2 * xk;

                float4v r = pk/qk;
                eps = select(qk != 0.0f, abs((t-r)/r), 1.0f);
                t = select(qk != 0.0f, r, t);

	            pkm2 = pkm1;
	            pkm1 = pk;
	            qkm2 = qkm1;
	            qkm1 = qk;

                float4v scale = select(abs(pk) > BIG, MACHEPF, float4v::ones());
                pkm2 *= scale;
                pkm1 *= scale;
                qkm2 *= scale;
                qkm1 *= scale;
            } while ((eps > MACHEPF  && greaterone).any_true());

            t *= exp(-x);
            y = select(greaterone, t, y);
        }


        mask4v lessone = x <= 1.0f && unprocessed;
        if (lessone.any_true()) {
            /*		Power series expansion		*/

            float4v psi = -EULER - log(x);

            for (int i = 1; i < n; ++i) {
	            psi += 1.0f/i;
            }

            float4v z = -x;
            float4v xk = 0.0;
            float4v yk = 1.0;
            float   pk = 1.0f - n;
            float4v t = n == 1  ?  0.0f : 1.0f/pk;
            float4v eps;

            do {
	            xk += 1.0;
	            yk *= z/xk;
	            pk += 1.0;
                t += select(float4v(pk) != 0.0f, yk/pk);

                eps = select(t != 0.0f, abs(yk/t), 1.0f);
            } while ((eps > MACHEPF && lessone).any_true());

            t = pow(z, n-1) * psi / gamma(float(n)) - t;
            y = select(lessone, t, y);
        }

        return y;
    }
}}}}

