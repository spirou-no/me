/*							ellpjf.c
 *
 *	Jacobian Elliptic Functions
 *
 *
 *
 * SYNOPSIS:
 *
 * float u, m, sn, cn, dn, phi;
 * int ellpj();
 *
 * ellpj( u, m, _&sn, _&cn, _&dn, _&phi );
 *
 *
 *
 * DESCRIPTION:
 *
 *
 * Evaluates the Jacobian elliptic functions sn(u|m), cn(u|m),
 * and dn(u|m) of parameter m between 0 and 1, and real
 * argument u.
 *
 * These functions are periodic, with quarter-period on the
 * real axis equal to the complete elliptic integral
 * ellpk(1.0-m).
 *
 * Relation to incomplete elliptic integral:
 * If u = ellik(phi,m), then sn(u|m) = sin(phi),
 * and cn(u|m) = cos(phi).  Phi is called the amplitude of u.
 *
 * Computation is by means of the arithmetic-geometric mean
 * algorithm, except when m is within 1e-9 of 0 or 1.  In the
 * latter case with m close to 1, the approximation applies
 * only for phi < pi/2.
 *
 * ACCURACY:
 *
 * Tested at random points with u between 0 and 10, m between
 * 0 and 1.
 *
 *            Absolute error (* = relative error):
 * arithmetic   function   # trials      peak         rms
 *    IEEE      sn          10000       1.7e-6      2.2e-7
 *    IEEE      cn          10000       1.6e-6      2.2e-7
 *    IEEE      dn          10000       1.4e-3      1.9e-5
 *    IEEE      phi         10000       3.9e-7*     6.7e-8*
 *
 *  Peak error observed in consistency check using addition
 * theorem for sn(u+v) was 4e-16 (absolute).  Also tested by
 * the above relation to the incomplete elliptic integral.
 * Accuracy deteriorates when u is large.
 *
 */

/*							ellpj.c		*/


/*
Cephes Math Library Release 2.2:  July, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/math/constants.h>
#include <cmath>


namespace oyrke { namespace algorithm { namespace sse { namespace math {

    void 
    ellpj(
        const float4v& u, 
        float m,
        float4v& sn, 
        float4v& cn, 
        float4v& dn, 
        float4v& ph
    ) {
        /* Check for special cases */

        if (m < 0.0f || m > 1.0f) {
            sn = cn = dn = ph = float4v::NaNs();
            return;
        }

        if (m < 1.0e-5) {
	        float4v t, b;
            sincos(u, t, b);
	        float4v ai = 0.25f * m * (u - t*b);
	        sn = t - ai*b;
	        cn = b + ai*t;
	        ph = u - ai;
	        dn = 1.0f - 0.5f*m*t*t;
	    }
        else if (m >= 0.99999f) {
	        float4v ai = 0.25f * (1.0f-m);
	        float4v b = cosh(u);
	        float4v t = tanh(u);
	        float4v phi = 1.0f/b;
	        float4v twon = b * sinh(u);
	        sn = t + ai * (twon - u)/(b*b);
	        ph = 2.0*atan(exp(u)) - PI_HALF + ai*(twon - u)/b;
	        ai *= t * phi;
	        cn = phi - ai * (twon - u);
	        dn = phi + ai * (twon + u);
	    }
        else {
            /*	A. G. M. scale		*/
            float a[10], c[10];
            int i = 0;
            float twon = 1.0;
            
            float b = std::sqrt(1.0f - m);

            c[0] = std::sqrt(m);
            a[0] = 1.0;

            while(std::abs((c[i]/a[i])) > MACHEPF) {
	            if( i > 8 ) {
                    /*	mtherr( "ellpjf", OVERFLOW );*/
		            break;
		        }
	            float ai = a[i];
	            ++i;
	            c[i] = 0.5f * (ai - b);
	            float t = std::sqrt(ai * b);
	            a[i] = 0.5f * (ai + b);
	            b = t;
	            twon += twon;
	        }

            /* backward recurrence */
            float4v phi = twon * a[i] * u;
            float4v t;
            float4v bphi;
            do {
	            t = c[i] * sin(phi) / a[i];
	            bphi = phi;
	            phi = 0.5 * (asin(t) + phi);
	        } while (--i);

            sincos(phi, sn, cn);
            dn = t/cos(phi-b);
            ph = phi;
        }
    }
 }}}}
