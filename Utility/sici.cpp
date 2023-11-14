/*							sicif.c
 *
 *	Sine and cosine integrals
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, Ci, Si;
 *
 * sicif( x, &Si, &Ci );
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the integrals
 *
 *                          x
 *                          -
 *                         |  cos t - 1
 *   Ci(x) = eul + ln x +  |  --------- dt,
 *                         |      t
 *                        -
 *                         0
 *             x
 *             -
 *            |  sin t
 *   Si(x) =  |  ----- dt
 *            |    t
 *           -
 *            0
 *
 * where eul = 0.57721566490153286061 is Euler's constant.
 * The integrals are approximated by rational functions.
 * For x > 8 auxiliary functions f(x) and g(x) are employed
 * such that
 *
 * Ci(x) = f(x) sin(x) - g(x) cos(x)
 * Si(x) = pi/2 - f(x) cos(x) - g(x) sin(x)
 *
 *
 * ACCURACY:
 *    Test interval = [0,50].
 * Absolute error, except relative when > 1:
 * arithmetic   function   # trials      peak         rms
 *    IEEE        Si        30000       2.1e-7      4.3e-8
 *    IEEE        Ci        30000       3.9e-7      2.2e-8
 */

/*
Cephes Math Library Release 2.1:  January, 1989
Copyright 1984, 1987, 1989 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/
#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

    /* TODO this implementation is the same as for double precision (in Cephes).
    Should customize for single precision.
    */

    const float SN[] = {
         1.00000000000000000302E0,
        -4.13470316229406538752E-2,
         9.76945438170435310816E-4,
        -9.75759303843632795789E-6,
         4.62591714427012837309E-8,
        -8.39167827910303881427E-11,
    };


    const float SD[] = {
        9.99999999999999996984E-1,
        1.42085239326149893930E-2,
        9.96412122043875552487E-5,
        4.41827842801218905784E-7,
        1.27997891179943299903E-9,
        2.03269266195951942049E-12,
    };

    const float CN[] = {
        -1.00000000000000000080E0,
         2.89159652607555242092E-2,
        -4.74007206873407909465E-4,
         3.59325051419993077021E-6,
        -1.35249504915790756375E-8,
         2.02524002389102268789E-11,
    };

    const float CD[] = {
        4.00000000000000000080E0,
        5.10028056236446052392E-2,
        3.17442024775032769882E-4,
        1.23210355685883423679E-6,
        3.06780997581887812692E-9,
        4.07746040061880559506E-12,
    };


    const float FN4[] = {
          5.48900223421373614008E-7,
          1.08936580650328664411E-4,
          6.81020132472518137426E-3,
          1.67006611831323023771E-1,
          1.62083287701538329132E0,
          5.45937717161812843388E0,
          4.23612862892216586994E0,
    };


    const float FD4[] = {
          5.48900252756255700982E-7,
          1.10034357153915731354E-4,
          7.01710668322789753610E-3,
          1.78792052963149907262E-1,
          1.86792257950184183883E0,
          7.30828822505564552187E0,
          8.16496634205391016773E0,
          1.00000000000000000000E0,
    };


    const float FN8[] = {
          9.70507110881952024631E-14,
          9.41779576128512936592E-11,
          3.20092790091004902806E-8,
          4.86215430826454749482E-6,
          3.49556442447859055605E-4,
          1.16064229408124407915E-2,
          1.60300158222319456320E-1,
          7.13715274100146711374E-1,
          4.55880873470465315206E-1,
    };

    const float FD8[] = {
          9.70507110881952025725E-14,
          9.43720590350276732376E-11,
          3.21956939101046018377E-8,
          4.92435064317881464393E-6,
          3.58696481881851580297E-4,
          1.22253594771971293032E-2,
          1.78685545332074536321E-1,
          9.17463611873684053703E-1,
          1.00000000000000000000E0,
    };

    const float GN4[] = {
          7.82579040744090311069E-9,
          1.97963874140963632189E-6,
          1.61999794598934024525E-4,
          5.38868681462177273157E-3,
          7.48527737628469092119E-2,
          3.97180296392337498885E-1,
          6.11379109952219284151E-1,
          8.71001698973114191777E-2,
    };

    const float GD4[] = {
          7.82579218933534490868E-9,
          2.02659182086343991969E-6,
          1.73221081474177119497E-4,
          6.22396345441768420760E-3,
          9.88771761277688796203E-2,
          6.66296701268987968381E-1,
          1.64402202413355338886E0,
          1.00000000000000000000E0,
    };

    const float GN8[] = {
          3.14040098946363334640E-15,
          3.85945925430276600453E-12,
          1.70404452782044526189E-9,
          3.47131167084116673800E-7,
          3.48941165502279436777E-5,
          1.71718239052347903558E-3,
          3.84878767649974295920E-2,
          3.30410979305632063225E-1,
          6.97359953443276214934E-1,
    };

    const float GD8[] = {
          3.14040098946363335242E-15,
          3.87830166023954706752E-12,
          1.72693748966316146736E-9,
          3.57043223443740838771E-7,
          3.68475504442561108162E-5,
          1.90284426674399523638E-3,
          4.67913194259625806320E-2,
          4.87852258695304967486E-1,
          1.68548898811011640017E0,
          1.00000000000000000000E0,
    };



    void sici(const float4v& xin, float4v& Si, float4v &Ci) {
        
        float4v x = abs(xin);
        mask4v nanarg = isnan(xin);
        float4v ysi = select(nanarg, xin);
        float4v yci = ysi;

        mask4v nullarg = x == 0.0f;
        // ysi = select(nullarg, 0.0f); -- not needed, 0 already
        yci = select(nullarg, -float4v::infinites(), yci);

        mask4v hugearg = x > 1.0e9;
        if (hugearg.any_true()) {
            // TODO is it valid to be here (x>1e9) for ingle precision???
            float4v s, c;
            sincos(x, s, c);
            ysi = blend(ysi, select(hugearg, PI_HALF - c/x));
            yci = blend(yci, select(hugearg, s/x));
        }

        mask4v openargs = !(nanarg || nullarg || hugearg);

        mask4v gt4arg = x > 4.0f  &&  openargs;
        if (gt4arg.any_true()) {
            /* The auxiliary functions are:
             *
             *
             * *si = *si - PIO2;
             * c = cos(x);
             * s = sin(x);
             *
             * t = *ci * s - *si * c;
             * a = *ci * c + *si * s;
             *
             * *si = t;
             * *ci = -a;
             */
            float4v s, c;
            sincos(x, s, c);

            float4v inv_x2 = 1.0f/(x*x);

            mask4v less8arg = x < 8.0f;
            float4v f1, g1, f2, g2;

            // 4 < x < 8
            if (less8arg.any_true()) {
                f1 = eval_polynomial(inv_x2, FN4) / (x*eval_polynomial(inv_x2, FD4));
                g1 = inv_x2*eval_polynomial(inv_x2, GN4) / eval_polynomial(inv_x2, GD4);
            }

            // x > 8
            if (less8arg.any_false()) {
                f2 = eval_polynomial(inv_x2, FN8) / (x*eval_polynomial(inv_x2, FD8));
                g2 = inv_x2*eval_polynomial(inv_x2, GN8) / eval_polynomial(inv_x2, GD8);
            }

            float4v f = select(less8arg, f1, f2);
            float4v g = select(less8arg, f2, g2);
            
            float4v tsi = PI_HALF - f*c - g*s;
            tsi = copysign(tsi, xin);
            ysi = select(gt4arg, tsi, ysi);

            
            float4v tci = f * s - g * c;
            yci = select(gt4arg, tci, yci);

            openargs = openargs && !gt4arg;
        }

        if (openargs.any_true()) {
            float4v x2 = x * x;
            float4v s = x * eval_polynomial(x2, SN) / eval_polynomial(x2, SD);
            float4v c = x2 * eval_polynomial(x2, CN) / eval_polynomial(x2, CD);

            s = copysign(s, xin);
            
            ysi = select(openargs, s, ysi);
            
            float4v x_safe = select(openargs, x, 1.0f);
            float4v tci = EULER + log(x_safe) + c;	/* real part if x < 0 */
            yci = select(openargs, tci, yci);
        }

        Si = ysi;
        Ci = yci;
    }
}}}}
