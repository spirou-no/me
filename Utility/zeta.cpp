/*							zetaf.c
 *
 *	Riemann zeta function of two arguments
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, q, y, zetaf();
 *
 * y = zetaf( x, q );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 *
 *                 inf.
 *                  -        -x
 *   zeta(x,q)  =   >   (k+q)  
 *                  -
 *                 k=0
 *
 * where x > 1 and q is not a negative integer or zero.
 * The Euler-Maclaurin summation formula is used to obtain
 * the expansion
 *
 *                n         
 *                -       -x
 * zeta(x,q)  =   >  (k+q)  
 *                -         
 *               k=1        
 *
 *           1-x                 inf.  B   x(x+1)...(x+2j)
 *      (n+q)           1         -     2j
 *  +  ---------  -  -------  +   >    --------------------
 *        x-1              x      -                   x+2j+1
 *                   2(n+q)      j=1       (2j)! (n+q)
 *
 * where the B2j are Bernoulli numbers.  Note that (see zetac.c)
 * zeta(x,1) = zetac(x) + 1.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,25        10000       6.9e-7      1.0e-7
 *
 * Large arguments may produce underflow in powf(), in which
 * case the results are inaccurate.
 *
 * REFERENCE:
 *
 * Gradshteyn, I. S., and I. M. Ryzhik, Tables of Integrals,
 * Series, and Products, p. 1073; Academic Press, 1980.
 *
 */

/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/math/constants.h>
#include <cmath>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    /* Expansion coefficients
     * for Euler-Maclaurin summation formula
     * (2k)! / B2k
     * where B2k are Bernoulli numbers
     */
    const float A[] = {
              12.0,
            -720.0,
           30240.0,
        -1209600.0,
        47900160.0,
       -1.8924375803183791606e9,    /*1.307674368e12/691*/
        7.47242496e10,
       -2.950130727918164224e12,    /*1.067062284288e16/3617*/
        1.1646782814350067249e14,   /*5.109094217170944e18/43867*/
       -4.5979787224074726105e15,   /*8.028576626982912e20/174611*/
        1.8152105401943546773e17,   /*1.5511210043330985984e23/854513*/
       -7.1661652561756670113e18    /*1.6938241367317436694528e27/236364091*/
    };
/* 30 Nov 86 -- error in third coefficient fixed */

    float4v zeta(const float4v& x, float q) {

        mask4v invalid_arg = x < 1.0f;   // --> NaN
        mask4v limit_arg   = x == 1.0f;  // --> inf
        // TODO manipulate x to avoid expensive special paths 
        /* Euler-Maclaurin summation formula */
        /*
        if( x < 25.0 )
        {
        */
        float4v s = pow(q, -x);
        float4v b;
        float4v w = q;
        const float4v ones = float4v::ones();

        bool done = false;

        for (int i=0; i<9 && !done; i++)	{
        	w += ones;
	        b = pow(w, -x);
	        s += b;
            done = (b/s < MACHEPF).all_true();
	    }

        if (!done) {
            float4v k = 0.0f;
            float4v a = 1.0f;
            s += b*w/(x-1.0f);
            s -= 0.5f * b;
        
            const int A_size = sizeof(A)/sizeof(A[0]);
            for (int i=0; i<A_size; i++) {
	            a *= x + k;
	            b /= w;
	            float4v t = a*b/A[i];
	            s = s + t;
	            t = abs(t/s);
                if ((t < MACHEPF).all_true()) {
                    break;
                }

                k += ones;
	            a *= x + k;
	            b /= w;
	            k += ones;
	        }
        }

        // TODO invalid/special arg
        return s;
    }
 }}}}

/* Basic sum of inverse powers */
/*
pseres:

s = powf( q, -x );
a = q;
do
	{
	a += 2.0;
	b = powf( a, -x );
	s += b;
	}
while( b/s > MACHEPF );

b = powf( 2.0, -x );
s = (s + b)/(1.0-b);
return(s);
*/
