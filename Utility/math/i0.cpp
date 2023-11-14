/*							i0f.c
 *
 *	Modified Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, i0();
 *
 * y = i0f( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns modified Bessel function of order zero of the
 * argument.
 *
 * The function is defined as i0(x) = j0( ix ).
 *
 * The range is partitioned into the two intervals [0,8] and
 * (8, infinity).  Chebyshev polynomial expansions are employed
 * in each interval.
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        100000      4.0e-7      7.9e-8
 *
 */
/*							i0ef.c
 *
 *	Modified Bessel function of order zero,
 *	exponentially scaled
 *
 *
 *
 * SYNOPSIS:
 *
 * float x, y, i0ef();
 *
 * y = i0ef( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns exponentially scaled modified Bessel function
 * of order zero of the argument.
 *
 * The function is defined as i0e(x) = exp(-|x|) j0( ix ).
 *
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE      0,30        100000      3.7e-7      7.0e-8
 * See i0f().
 *
 */

/*							i0.c		*/


/*
Cephes Math Library Release 2.2:  June, 1992
Copyright 1984, 1987, 1992 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#include <utility/math.h>
#include <utility/sse.h>
//#include <utility/math/constants.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	/* Chebyshev coefficients for exp(-x) I0(x)
	 * in the interval [0,8].
	 *
	 * lim(x->0){ exp(-x) I0(x) } = 1.
	 */

	static const float A[] = {
 		 6.76795274409476084995E-1f,
		-3.04682672343198398683E-1f,
		 1.71620901522208775349E-1f,
		-9.49010970480476444210E-2f,
		 4.93052842396707084878E-2f,
		-2.37374148058994688156E-2f,
		 1.05464603945949983183E-2f,
		-4.32430999505057594430E-3f,
		 1.63947561694133579842E-3f,
		-5.76375574538582365885E-4f,
		 1.88502885095841655729E-4f,
		-5.75419501008210370398E-5f,
		 1.64484480707288970893E-5f,
		-4.41673835845875056359E-6f,
		 1.11738753912010371815E-6f,
		-2.67079385394061173391E-7f,
		 6.04699502254191894932E-8f,
		-1.30002500998624804212E-8f,
	};


	/* Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
	 * in the inverted interval [8,infinity].
	 *
	 * lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
	 */

	static const float B[] = {
	 8.04490411014108831608E-1f,
	 3.36911647825569408990E-3f,
	 6.88975834691682398426E-5f,
	 2.89137052083475648297E-6f,
	 2.04891858946906374183E-7f,
	 2.26666899049817806459E-8f,
	 3.39623202570838634515E-9f
	};

 
	float4v i0(const float4v& xin) {
		return exp(xin) * i0e(xin);
	}


	float4v
	i0e(const float4v& xin) {
		
		float4v x = abs(xin);

		mask4v smallarg = x <= 8.0;
		mask4v bigarg = x > 8.0;
		float4v ysmall;
		float4v ybig;
		if (smallarg.any_true()) {
			float4v y = 0.5f*x - 2.0f;
			ysmall = eval_chebychev(y, A);
		}

		if (bigarg.any_true()) {
			float4v y = 32.0f/x - 2.0f;
			ybig = eval_chebychev(y, B) / sqrt(x);
		}

		float4v result = blend(ysmall, ybig);
		result = select(isnan(x), xin, result);

		return result;
	}

/*							chbevlf.c
 *
 *	Evaluate Chebyshev series
 *
 *
 *
 * SYNOPSIS:
 *
 * int N;
 * float x, y, coef[N], chebevlf();
 *
 * y = chbevlf( x, coef, N );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the series
 *
 *        N-1
 *         - '
 *  y  =   >   coef[i] T (x/2)
 *         -            i
 *        i=0
 *
 * of Chebyshev polynomials Ti at argument x/2.
 *
 * Coefficients are stored in reverse order, i.e. the zero
 * order term is last in the array.  Note N is the number of
 * coefficients, not the order.
 *
 * If coefficients are for the interval a to b, x must
 * have been transformed to x -> 2(2x - b - a)/(b-a) before
 * entering the routine.  This maps x from (a, b) to (-1, 1),
 * over which the Chebyshev polynomials are defined.
 *
 * If the coefficients are for the inverted interval, in
 * which (a, b) is mapped to (1/b, 1/a), the transformation
 * required is x -> 2(2ab/x - b - a)/(b-a).  If b is infinity,
 * this becomes x -> 4a/x - 1.
 *
 *
 *
 * SPEED:
 *
 * Taking advantage of the recurrence properties of the
 * Chebyshev polynomials, the routine requires one more
 * addition per loop than evaluating a nested polynomial of
 * the same degree.
 *
 */
/*							chbevl.c	*/

/*
Cephes Math Library Release 2.0:  April, 1987
Copyright 1985, 1987 by Stephen L. Moshier
Direct inquiries to 30 Frost Street, Cambridge, MA 02140
*/

#if 0
    template <typename U, typename T, size_t N>
    __forceinline T 
    eval_chebychev(const T& x, const U (&c)[N]) {

        typedef std::pair<T, T> recurse_t;

        recurse_t b = chev(x, c+1);

        recurse_t chev(const T& x, const U (&c)[N]) {
            recurse_t b = chev(x, c+1);
            T b0 = x*b.first - b.second + c[0];
            return std::make_pair(b0, b.first);
        }

        recurse_t chev<0>(const T& x, const U (&c)[N]) {
            return std::make_pair(c[0], T(0));
        }

        // special case:
        // N
        // 0    c[0]
        chev(x, c+1, b1, b2);
        return c[0] -b2 + b1*x;

        void chev(const T& x, const U (&c)[N], T& b0, T& b1) {
            T d1, d2;
            chev(x, c+1, d1, d2);
            b0 = x*d1 - d2 + c[0]; // why not 2*x*d1? // b0 = (x+x)*d1 - d2 + c[0];
            b1 = d1;
        }

        void chev<0>(const T& x, const U (&c)[N], T& b0, T& b1) {
            b0 = c[0];
            b1 = T(0);
        }

        /*
        function p = ClenshawChebyshev( t, c )
        % Clenshaw algorithm for evaluation of the series
        m = length(c) - 1;
        if m<0,
            % empty c ==> zero polynomial
            p = zeros(size(t));
        elseif m==0,
            % degree 0 ==> constant polynomial
            p = c(1)*ones(size(t));
        else
            % degree 1 or more, use Clenshaw algorithm
            y2 = zeros(size(t));  % y_{m+1}
            y1 = c(m+1)*ones(size(t));  % y_{m}
            % 
            % p = sum( c(j)*T(j-1,t), j=1..m-1) + c(m)*T(m-1,t) 
            %        + y1*T(m,t) - y2*T(m-1,t)  Loop invariant initialization
            for k=m:-1:2,
                % p = sum( c(j)*T(j-1,t), j=1..k-1) + c(k)*T(k-1,t) + y1*T(k,t) 
                %  - y2*T(k-1,t)                Loop invariant
                y0 = c(k) +2*t.*y1 - y2;   % y_{k-1} 
                % p = sum( c(j)*T(j-1,t), j=1..k-1) + y0*T(k-1,t) +
                % y1*(T(k,t)-2tT(k-1,t)) [= -y1*T(k-2,t) ]
                y2 = y1; 
                y1 = y0; 
                % p = sum( c(j)*T(j-1,t), j=1..k-2) + c(k-1)*T(k-2,t) + 
                % y1*T(k-1,t) - y2*T(k-2,t)     Proving correctness
            end
            % p = c(1)*T(0,t) + y1*T(1,t) - y2*T(0,t)
            p = c(1)-y2 + y1.*t;
        end
        */

        // check if 0.5*(b0-b2) is needed in recursive formula.
        // 
	U b0 = c[0];
	U b1 = U(0);
	U b2;
	size_t index = 1;
	int i = N-1;

	do {
		b2 = b1;
		b1 = b0;
		b0 = x * b1  - b2 + c[index];
		++index;
	} while (--i);

	return 0.5f*(b0-b2);
  }
  

#ifdef ANSIC
float chbevlf( float x, float *array, int n )
#else
float chbevlf( x, array, n )
float x;
float *array;
int n;
#endif
{
float b0, b1, b2, *p;
int i;

p = array;
b0 = *p++;
b1 = 0.0;
i = n - 1;

do
	{
	b2 = b1;
	b1 = b0;
	b0 = x * b1  -  b2  + *p++;
	}
while( --i );

return( 0.5*(b0-b2) );
}
#endif
}}}}