/*							bdtrf.c
 *
 *	Binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * float p, y, bdtrf();
 *
 * y = bdtrf( k, n, p );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms 0 through k of the Binomial
 * probability density:
 *
 *   k
 *   --  ( n )   j      n-j
 *   >   (   )  p  (1-p)
 *   --  ( j )
 *  j=0
 *
 * The terms are not summed directly; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = bdtr( k, n, p ) = incbet( n-k, k+1, 1-p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error (p varies from 0 to 1):
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       2000       6.9e-5      1.1e-5
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * bdtrf domain        k < 0            0.0
 *                     n < k
 *                     x < 0, x > 1
 *
 */
/*							bdtrcf()
 *
 *	Complemented binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * float p, y, bdtrcf();
 *
 * y = bdtrcf( k, n, p );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms k+1 through n of the Binomial
 * probability density:
 *
 *   n
 *   --  ( n )   j      n-j
 *   >   (   )  p  (1-p)
 *   --  ( j )
 *  j=k+1
 *
 * The terms are not summed directly; instead the incomplete
 * beta integral is employed, according to the formula
 *
 * y = bdtrc( k, n, p ) = incbet( k+1, n-k, p ).
 *
 * The arguments must be positive, with p ranging from 0 to 1.
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error (p varies from 0 to 1):
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       2000       6.0e-5      1.2e-5
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * bdtrcf domain     x<0, x>1, n<k       0.0
 */
/*							bdtrif()
 *
 *	Inverse binomial distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k, n;
 * float p, y, bdtrif();
 *
 * p = bdtrf( k, n, y );
 *
 *
 *
 * DESCRIPTION:
 *
 * Finds the event probability p such that the sum of the
 * terms 0 through k of the Binomial probability density
 * is equal to the given cumulative probability y.
 *
 * This is accomplished using the inverse beta integral
 * function and the relation
 *
 * 1 - p = incbi( n-k, k+1, y ).
 *
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error (p varies from 0 to 1):
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       2000       3.5e-5      3.3e-6
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * bdtrif domain    k < 0, n <= k         0.0
 *                  x < 0, x > 1
 *
 */

/*								bdtr() */


/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>
#include <cmath>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    float4v bdtrc(const int4v& k, const int4v& n, const float4v& p) {

        mask4v invalid = isnan(p) || p < 0.0f || p > 1.0f;
        float4v y = select(invalid, float4v::NaNs());

        mask4v openarg = !invalid;
        mask4v neg_k = k < 0;
        y = select(neg_k && openarg, 1.0f, y);
        openarg = openarg && !neg_k;

        // TODO return 0 or NaN here?
        // n == k ==> 0
        // n < k  ==> NaN
        mask4v n_less_k = n <= k;
        y = select(n_less_k && openarg, 0.0f, y);
        openarg = openarg && !n_less_k;

        if (openarg.any_true()) {
            float4v dn(n - k);

            // and zero_k with args left to process,
            // don't want expensive compute here if we're
            // ignoring the result later anyway.
            mask4v zero_k = k == 0;
            float4v dk_zero = (zero_k && openarg).any_true()
                            ? 1.0f - pow(1.0f-p, dn)
                            : float4v::zeros();
            float4v dk_pos  = (!zero_k && openarg).any_true()
                            ? incbet(float4v(k+1), dn, p)
                            : float4v::zeros();

            float4v dk = select(zero_k, dk_zero, dk_pos);
            y = select(openarg, dk, y);
        }

        return y;
    }



    float4v bdtr(const int4v& k, const int4v& n, const float4v& p) {

        // why not
        // 1.0f - bdtrc(k, n, 1-0f - p)
        mask4v invalid = isnan(p) || p < 0.0f || p > 1.0f;
        float4v y = select(invalid, float4v::NaNs());

        mask4v openarg = !invalid;
        mask4v neg_k = k < 0;
        y = select(neg_k && openarg, 1.0f, y);
        openarg = openarg && !neg_k;

        // TODO return 0 or NaN here?
        mask4v n_less_k = n < k;
        y = select(n_less_k && openarg, 0.0f, y);
        openarg = openarg && !n_less_k;

        mask4v n_equal_k = n == k;
        y = select(n_equal_k && openarg, 1.0f, y);
        openarg = openarg && !n_equal_k;

        if (openarg.any_true()) {
            float4v dn(n - k);

            // and zero_k with args left to process,
            // don't want expensive compute here if we're
            // ignoring the result later anyway.
            mask4v zero_k = k == 0;
            float4v dk_zero = (zero_k && openarg).any_true()
                            ? pow(1.0f-p, dn)
                            : float4v::zeros();
            float4v dk_pos  = (!zero_k && openarg).any_true()
                            ? incbet(dn, float4v(k+1), 1.0 - p)
                            : float4v::zeros();

            float4v dk = select(zero_k, dk_zero, dk_pos);
            y = select(openarg, dk, y);
        }

        return y;
    }




    float4v bdtri(const int4v& k, const int4v& n, const float4v& p) {
        mask4v invalid = isnan(p) || p < 0.0f || p > 1.0f;
        float4v y = select(invalid, float4v::NaNs());

        mask4v openarg = !invalid;
        mask4v invalid_k = k < 0 || n <= k;
        y = select(invalid_k && openarg, 0.0f, y);
        openarg = openarg && !invalid_k;

        if (openarg.any_true()) {
            float4v dn(n - k);

            // and zero_k with args left to process,
            // don't want expensive compute here if we're
            // ignoring the result later anyway.
            mask4v zero_k = k == 0;
            float4v p_zero = (zero_k && openarg).any_true()
                            ? 1.0f - pow(p, 1.0f/dn)
                            : float4v::zeros();
            float4v p_pos  = (!zero_k && openarg).any_true()
                            ? 1.0f - incbi(dn, float4v(k+1), p)
                            : float4v::zeros();

            float4v pp = select(zero_k, p_zero, p_pos);
            y = select(openarg, pp, y);
        }

        return y;
    }
}}}}