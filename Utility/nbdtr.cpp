/*							nbdtrf.c
 *
 *	Negative binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * float p, y, nbdtrf();
 *
 * y = nbdtrf( k, n, p );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms 0 through k of the negative
 * binomial distribution:
 *
 *   k
 *   --  ( n+j-1 )   n      j
 *   >   (       )  p  (1-p)
 *   --  (   j   )
 *  j=0
 *
 * In a sequence of Bernoulli trials, this is the probability
 * that k or fewer failures precede the nth success.
 *
 * The terms are not computed individually; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = nbdtr( k, n, p ) = incbet( n, k+1, p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       5000       1.5e-4      1.9e-5
 *
 */
/*							nbdtrcf.c
 *
 *	Complemented negative binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * float p, y, nbdtrcf();
 *
 * y = nbdtrcf( k, n, p );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms k+1 to infinity of the negative
 * binomial distribution:
 *
 *   inf
 *   --  ( n+j-1 )   n      j
 *   >   (       )  p  (1-p)
 *   --  (   j   )
 *  j=k+1
 *
 * The terms are not computed individually; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = nbdtrc( k, n, p ) = incbet( k+1, n, 1-p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       5000       1.4e-4      2.0e-5
 *
 */

/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {
namespace {
    float4v nbdtr_impl(const int4v& k, const int4v& n, const float4v& p, const float4v& p_incbet) {

        mask4v nanarg = isnan(p);
        float4v y = select(nanarg, p);

        mask4v validarg = !nanarg  &&  p >= 0.0f  &&  p <= 1.0f  &&  k >= 0;

        float4v p_mod = select(validarg, p_incbet, float4v::NaNs());

        float4v yn = validarg.any_true() 
                   ? incbet(float4v(k+1), float4v(n), p_mod)
                   : float4v::zeros();

        y = select(validarg, yn, y);
        return y;
    }

}

    float4v nbdtrc(const int4v& k, const int4v& n, const float4v& p) {
        float4v y = nbdtr_impl(k, n, p, 1.0f-p);
        return y;
    }



    float4v nbdtr(const int4v& k, const int4v& n, const float4v& p) {
        float4v y = nbdtr_impl(k, n, p, p);
        return y;
    }

 }}}}

