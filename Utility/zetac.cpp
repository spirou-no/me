 /*							zetacf.c
 *
 *	Riemann zeta function
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, zetacf();
 *
 * y = zetacf( x );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *
 *                inf.
 *                 -    -x
 *   zetac(x)  =   >   k   ,   x > 1,
 *                 -
 *                k=2
 *
 * is related to the Riemann zeta function by
 *
 *	Riemann zeta(x) = zetac(x) + 1.
 *
 * Extension of the function definition for x < 1 is implemented.
 * Zero is returned for x > log2(MAXNUM).
 *
 * An overflow error may occur for large negative x, due to the
 * gamma function in the reflection formula.
 *
 * ACCURACY:
 *
 * Tabulated values have full machine accuracy.
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      1,50        30000       5.5e-7      7.5e-8
 *
 *
 */

/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1989, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/math/constants.h>
//#include <cmath>
#include <limits>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    /* Riemann zeta(x) - 1
     * for integer arguments between 0 and 30.
     */
    const float azetac[] = {
        -1.50000000000000000000E0,
         std::numeric_limits<float>::infinity(), //1.70141183460469231730E38, /* infinity. */
         6.44934066848226436472E-1,
         2.02056903159594285400E-1,
         8.23232337111381915160E-2,
         3.69277551433699263314E-2,
         1.73430619844491397145E-2,
         8.34927738192282683980E-3,
         4.07735619794433937869E-3,
         2.00839282608221441785E-3,
         9.94575127818085337146E-4,
         4.94188604119464558702E-4,
         2.46086553308048298638E-4,
         1.22713347578489146752E-4,
         6.12481350587048292585E-5,
         3.05882363070204935517E-5,
         1.52822594086518717326E-5,
         7.63719763789976227360E-6,
         3.81729326499983985646E-6,
         1.90821271655393892566E-6,
         9.53962033872796113152E-7,
         4.76932986787806463117E-7,
         2.38450502727732990004E-7,
         1.19219925965311073068E-7,
         5.96081890512594796124E-8,
         2.98035035146522801861E-8,
         1.49015548283650412347E-8,
         7.45071178983542949198E-9,
         3.72533402478845705482E-9,
         1.86265972351304900640E-9,
         9.31327432419668182872E-10
    };


    /* 2**x (1 - 1/x) (zeta(x) - 1) = P(1/x)/Q(1/x), 1 <= x <= 10 */
    const float P[] = {
        2.01822444485997955865E2,
        1.51129169964938823117E4,
        5.92785467342109522998E5,
        1.60837006880656492731E7,
        3.41646073514754094281E8,
        5.15399538023885770696E9,
        4.87781159567948256438E10,
        2.57534127756102572888E11,
        5.85746514569725319540E11
    };

    const float Q[] = {
        1.96436237223387314144E2,
        1.60382976810944131506E4,
        5.66666825131384797029E5,
        1.79410371500126453702E7,
        3.39006746015350418834E8,
        5.64451517271280543351E9,
        5.22858235368272161797E10,
        3.90497676373371157516E11,
        1.00000000000000000000E0
    };

    /* log(zeta(x) - 1 - 2**-x), 10 <= x <= 50 */
    const float A[] = {
         9.26786275768927717187E16,
         7.82905376180870586444E16,
        -9.92763810039983572356E16,
        -1.98123688133907171455E15,
         5.13778997975868230192E15,
         3.31884402932705083599E14,
         2.26888156119238241487E13,
         5.29806374009894791647E11,
         2.60889506707483264896E10,
         1.76506865670346462757E8,
         8.70728567484590192539E6
    };
    
    const float B[] = {
        -1.79915597658676556828E16,
         5.71464111092297631292E16,
         5.34589509675789930199E15,
        -4.86299103694609136686E15,
        -2.96075404507272223680E14,
        -2.07820961754173320170E13,
        -4.80319584350455169857E11,
        -2.37669260975543221788E10,
        -1.60529969932920229676E8,
        -7.92625410563741062861E6,
         1.00000000000000000000E0
    };

    /* (1-x) (zeta(x) - 1), 0 <= x <= 1 */
    const float R[] = {
        -1.11578094770515181334E5,
         1.26726061410235149405E4,
         1.01050368053237678329E3,
        -2.48762831680821954401E2,
         1.55162528742623950834E1,
        -3.28717474506562731748E-1
    };

    const float S[] = {
        7.43853965136767874343E4,
        2.03665876435770579345E4,
        3.03835500874445748734E3,
        3.17710311750646984099E2,
        1.95107674914060531512E1,
        1.00000000000000000000E0,
    };


    #define MAXL2 127

    /*
     * Riemann zeta function, minus one
     */


    float4v zetac(const float4v& x) {
        
        mask4v negarg = x < 0.0f;
        mask4v openarg = mask4v::trues();
        float4v y;

        if (negarg.any_true()) {
            // make sure s have no negative elements.  these will not be used nayway, 
            // but need to avoid recursion
	        float4v s = abs(1.0f - x);  
            float4v w = zetac(s);
	        float4v b = sin(PI_HALF*x) * pow(2.0f*PI, x) * gamma(s) * (1.0 + w) / PI;
            y = select(x >= -30.8148, b - 1.0f);
            openarg = !negarg;
	    }

        // mask4v bigarg = x >= MAXL2;  -- no need to handle this case specifically

        mask4v nanarg = isnan(x);
        y = select(nanarg, x, y);
        openarg = openarg && !nanarg;

        float4v w = floor(x);
        mask4v intarg = w == x;
        
        /* Tabulated values for integer argument */
        if (intarg.any_true()) {
	        int4v ix = trunc_int(x);
            // TODOclip index to 30
            int4v clip_ix = min(ix, int(array_sizeof(azetac)-1)); //azetac_size-1);
            intarg = intarg && (clip_ix == ix);
            float4v yint = float4v::create_indexed(azetac, clip_ix);  // see create_index in pow.cpp
            y = select(intarg, yint, y);
            openarg = openarg && !intarg;
	    }

        mask4v less1arg = (x < 1.0f) && openarg;
        
        if (less1arg.any_true()) {
	        float4v w = 1.0f - x;
	        float4v yless1 = eval_polynomial(x, R) / ( (1.0f-x) * eval_polynomial(x, S));
            y = select(less1arg, yless1, y);
            openarg = openarg && !less1arg;
	    }

        mask4v equal1arg = x == 1.0f;
        float4v yequal1arg = float4v::infinites();

        mask4v less10arg = (x <= 10.0f) && openarg;
        
        if (less10arg.any_true()) {
	        float4v b = exp2(x) * (x-1.0f); // pow( 2.0, x ) * (x - 1.0);
	        float4v w = 1.0f/x;
	        float4v yless10 = (x * eval_polynomial(w, P)) / (b * eval_polynomial(w, Q));
            y = select(less10arg, yless10, y);
            openarg = openarg && !less10arg;
	    }

        mask4v less50arg = x <= 50.0f;

        if (less50arg.any_true()) {
	        float4v b = exp2(-x); //powf( 2.0, -x );
	        float4v w = eval_polynomial(x, A) / eval_polynomial(x, B);
	        float4v yless50 = exp(w) + b;
            y = select(less50arg, yless50, y);
            openarg = openarg && !less50arg;
	    }

        if (openarg.any_true()) {
            /* Basic sum of inverse powers */
            float4v s = 0.0f;
            float4v a = 1.0;
            bool not_done = false;
            do {
	            a += 2.0;
	            float4v b = pow(a, -x);
	            s += b;
                not_done = (b/s > MACHEPF).any_true();
	        } while(not_done);

            float4v b = exp2(-x); //powf( 2.0, -x );
            float4v ygreater50 = (s + b)/(1.0f-b);

            y = select(openarg, ygreater50, y);
        }

        return y;
    }
}}}}