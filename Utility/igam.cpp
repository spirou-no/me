/*							igamf.c
 *
 *	Incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * float a, x, y, igamf();
 *
 * y = igamf( a, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *                           x
 *                            -
 *                   1       | |  -t  a-1
 *  igam(a,x)  =   -----     |   e   t   dt.
 *                  -      | |
 *                 | (a)    -
 *                           0
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        20000       7.8e-6      5.9e-7
 *
 */
/*							igamcf()
 *
 *	Complemented incomplete gamma integral
 *
 *
 *
 * SYNOPSIS:
 *
 * float a, x, y, igamcf();
 *
 * y = igamcf( a, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * The function is defined by
 *
 *
 *  igamc(a,x)   =   1 - igam(a,x)
 *
 *                            inf.
 *                              -
 *                     1       | |  -t  a-1
 *               =   -----     |   e   t   dt.
 *                    -      | |
 *                   | (a)    -
 *                             x
 *
 *
 * In this implementation both arguments must be positive.
 * The integral is evaluated by either a power series or
 * continued fraction expansion, depending on the relative
 * values of a and x.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        30000       7.8e-6      5.9e-7
 *
 */

/*
Cephes Math Library Release 2.2: June, 1992
Copyright 1985, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {
    namespace {
        float4v processArgs(const float4v& a, const float4v& xin, float valueNegative, mask4v& unprocessed) {
            float4v x = select(unprocessed, xin, float4v::NaNs());
            mask4v nanarg = isnan(a);
            float4v y = select(nanarg, a);
        
            mask4v openarg = nanarg;
            nanarg = isnan(x);
            y = select(nanarg, x, y);
            openarg = !(openarg || nanarg);

            mask4v negarg = (x <= 0.0f || a <= 0.0f) && openarg;
            y = select(negarg, valueNegative, y);
            
            unprocessed = openarg && !negarg;
            return y;
        }


        float4v igam_impl(const float4v& a, const float4v& x, const mask4v& unprocessed) {
            float4v ax = a * log(x) - x - lgam(a);
            mask4v underflow = ax < -MAXLOG;
            // y = select(small, 0.0f, y); -- not needed, already 0
            mask4v openargs = unprocessed && !underflow;

            ax = exp(ax);
            /* power series */
            float4v r = a;
            float4v c = 1.0f;
            float4v ans = 1.0f;

            do {
                r += 1.0f;
                c *= x/r;
                ans += c;
            } while (((c/ans > MACHEPF) && openargs).any_true());
            return ans*ax/a;
        }


            
         float4v igamc_impl(const float4v& a, const float4v& xin, const mask4v& unprocessed) {
             /* BIG = 1/MACHEPF */
             const float BIG =  16777216.0f; // floats larger than this are ints
            
             float4v x = select(unprocessed, xin, float4v::NaNs());
             float4v ax = a * log(x) - x - lgam(a);
             mask4v underflow = (ax < -MAXLOG) && unprocessed;
             // y = select(underflow, float4v::zeros(), y); -- not needed
             mask4v openargs = unprocessed && !underflow;

             float4v y;

             if (openargs.any_true()) {
                  ax = exp(ax);
                 /* continued fraction */

                 float4v t = 1.0f - a;
                 float4v z = x + t + 1.0f;
                 float4v c;
                 float4v pkm2 = 1.0f;
                 float4v qkm2 = x;
                 float4v pkm1 = x + 1.0f;
                 float4v qkm1 = z * x;
                 float4v ans = pkm1 / qkm1;
                 float4v err = 1.0f;
                  
                 do {
                     c += 1.0f;
                     t += 1.0f;
                     z += 2.0f;
                     float4v yc = y*c;
                     float4v pk = pkm1 * z  -  pkm2 * yc;
                     float4v qk = qkm1 * z  -  qkm2 * yc;
                     mask4v qk_nzero = qk != 0.0f;
                     float4v r = pk/qk;
                     err = select(qk_nzero, (err-r)/r, 1.0f);
                     ans = select(qk_nzero, r, ans);
                      
 	                 pkm2 = pkm1;
 	                 pkm1 = pk;
	                 qkm2 = qkm1;
	                 qkm1 = qk;

                     mask4v pk_big = abs(pk) > BIG;
                     float4v scale = select(pk_big, MACHEPF, float4v::ones());
                     pkm2 *= scale;
                     pkm1 *= scale;
                     qkm2 *= scale;
                     qkm1 *= scale;
                 } while (((err > MACHEPF) && openargs).any_true());

                 y = select(openargs, ans * ax);
            }

            return y;
        }
    }



    float4v igamc(const float4v& a, const float4v& x) {

        mask4v openarg;
        float4v y = processArgs(a, x, 1.0f, openarg);
         
        mask4v process = (x < 1.0f || x < a) && openarg;
        float4v t = process.any_true() 
                  ? 1.0f - igam_impl(a, x, process)
                  : float4v::zeros();
        y = select(process, t, y);
        openarg = openarg && !process;

        if (openarg.any_true()) {
            float4v yn = igamc_impl(a, x, openarg);
            y = select(openarg, yn, y);
        }

        return y;
    }



    /* left tail of incomplete gamma function:
     *
     *          inf.      k
     *   a  -x   -       x
     *  x  e     >   ----------
     *           -     -
     *          k=0   | (a+k+1)
     *
     */

    float4v igam(const float4v& a, const float4v& x) {
        mask4v openarg;
        float4v y = processArgs(a, x, 0.0f, openarg);
         
        mask4v bigarg = (x>1.0f  && x > a) && openarg;
        float4v t = bigarg.any_true() 
                    ? 1.0f - igamc_impl(a, x, openarg)
                    : float4v::zeros();
        y = select(bigarg, t, y);
        openarg = openarg && !bigarg;

        if (openarg.any_true()) {
            float4v yn = igam_impl(a, x, openarg);
            y = select(openarg, yn, y);
        }

        return y;
    }
}}}}
