/*							pdtrf.c
 *
 *	Poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * float m, y, pdtrf();
 *
 * y = pdtrf( k, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the first k terms of the Poisson
 * distribution:
 *
 *   k         j
 *   --   -m  m
 *   >   e    --
 *   --       j!
 *  j=0
 *
 * The terms are not summed directly; instead the incomplete
 * gamma integral is employed, according to the relation
 *
 * y = pdtr( k, m ) = igamc( k+1, m ).
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       5000       6.9e-5      8.0e-6
 *
 */
/*							pdtrcf()
 *
 *	Complemented poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * float m, y, pdtrcf();
 *
 * y = pdtrcf( k, m );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the sum of the terms k+1 to infinity of the Poisson
 * distribution:
 *
 *  inf.       j
 *   --   -m  m
 *   >   e    --
 *   --       j!
 *  j=k+1
 *
 * The terms are not summed directly; instead the incomplete
 * gamma integral is employed, according to the formula
 *
 * y = pdtrc( k, m ) = igam( k+1, m ).
 *
 * The arguments must both be positive.
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       5000       8.4e-5      1.2e-5
 *
 */
/*							pdtrif()
 *
 *	Inverse Poisson distribution
 *
 *
 *
 * SYNOPSIS:
 *
 * int k;
 * float m, y, pdtrf();
 *
 * m = pdtrif( k, y );
 *
 *
 *
 *
 * DESCRIPTION:
 *
 * Finds the Poisson variable x such that the integral
 * from 0 to x of the Poisson density is equal to the
 * given probability y.
 *
 * This is accomplished using the inverse gamma integral
 * function and the relation
 *
 *    m = igami( k+1, y ).
 *
 *
 *
 *
 * ACCURACY:
 *
 *        Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE       0,100       5000       8.7e-6      1.4e-6
 *
 * ERROR MESSAGES:
 *
 *   message         condition      value returned
 * pdtri domain    y < 0 or y >= 1       0.0
 *                     k < 0
 *
 */

/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    namespace {
        typedef float4v (gamma_fn)(const float4v&, const float4v&);

        float4v pdtrc_impl(const int4v& k, const float4v& x, const mask4v& x_valid, 
                           gamma_fn gamma) {
            mask4v nanarg = isnan(x);
            float4v y = select(nanarg, x);
        
            // no need to modify y if domainerr, should be 0.0 anyway
            mask4v validarg = !nanarg   &&   k >= 0   &&   x_valid;
        
            // note mod'f x arg to igam to avoid unused inputs don't trigger expensive
            // compate paths in igam
            float4v y_normal = validarg.any_true()
                             ? (*gamma)(float4v(k+1), select(validarg, x, float4v::NaNs()))
                             : float4v::zeros();
            y = select(validarg, y_normal, y);

            return y;
        }
    }



    float4v pdtrc(const int4v& k, const float4v& x) {
        float4v y = pdtrc_impl(k, x, x>0.0f, igam);
        return y;
    }



    float4v pdtr(const int4v& k, const float4v& x) {
        float4v y = pdtrc_impl(k, x, x>0.0f, igamc);
        return y;
    }


    
    
    float4v pdtri(const int4v& k, const float4v& x) {
        float4v y = pdtrc_impl(k, x, x>=0.0f && x<1.0f, igami);
        return y;
    }
 }}}}
