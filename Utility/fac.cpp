/*							facf.c
 *
 *	Factorial function
 *
 *
 *
 * SYNOPSIS:
 *
 * float y, facf();
 * int i;
 *
 * y = facf( i );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns factorial of i  =  1 * 2 * 3 * ... * i.
 * fac(0) = 1.0.
 *
 * Due to machine arithmetic bounds the largest value of
 * i accepted is 33 in single precision arithmetic.
 * Greater values, or negative ones,
 * produce an error message and return MAXNUM.
 *
 *
 *
 * ACCURACY:
 *
 * For i < 34 the values are simply tabulated, and have
 * full machine accuracy.
 *
 */

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1984, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
#include <limits>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

    /* Factorials of integers from 0 through 33, __offset__ by 1 */
    const float factbl[] = {
        std::numeric_limits<float>::infinity(),
        1.00000000000000000000E0,
        1.00000000000000000000E0,
        2.00000000000000000000E0,
        6.00000000000000000000E0,
        2.40000000000000000000E1,
        1.20000000000000000000E2,
        7.20000000000000000000E2,
        5.04000000000000000000E3,
        4.03200000000000000000E4,
        3.62880000000000000000E5,
        3.62880000000000000000E6,
        3.99168000000000000000E7,
        4.79001600000000000000E8,
        6.22702080000000000000E9,
        8.71782912000000000000E10,
        1.30767436800000000000E12,
        2.09227898880000000000E13,
        3.55687428096000000000E14,
        6.40237370572800000000E15,
        1.21645100408832000000E17,
        2.43290200817664000000E18,
        5.10909421717094400000E19,
        1.12400072777760768000E21,
        2.58520167388849766400E22,
        6.20448401733239439360E23,
        1.55112100433309859840E25,
        4.03291461126605635584E26,
        1.0888869450418352160768E28,
        3.04888344611713860501504E29,
        8.841761993739701954543616E30,
        2.6525285981219105863630848E32,
        8.22283865417792281772556288E33,
        2.6313083693369353016721801216E35,
        8.68331761881188649551819440128E36,
        std::numeric_limits<float>::infinity()
    };

    
    float4v fac(const int4v& i) {
        // add 1 to i since we're putting infinity at 0th index to represent negative i, 
        // fac(0) is at index 1
        static const int max_index = int(array_sizeof(factbl)-1);

        int4v mod_index = clip(i+1, int4v::zeros(), max_index);
        float4v y = float4v::create_indexed(factbl, mod_index);
        return y;
    }
}}}}
