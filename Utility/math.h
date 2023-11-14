#pragma once

#include <utility/sse.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {
	
	float sin_cephes(float xx );
	float4v cbrt_halleys(const float4v& x);																											// arc sine, returned angle in range -pi/2..pi/2


	// triginometric functions
	float4v acos(const float4v& x);																											// arc sine, returned angle in range -pi/2..pi/2
	float4v asin(const float4v& x);																											// arc sine, returned angle in range -pi/2..pi/2
	float4v atan(const float4v& x);
	float4v atan2(const float4v& y, const float4v& x);
	float4v cos(const float4v& x);
	float4v cosm1(const float4v& x);																										// cos(x)-1
	float4v cosdg(const float4v& x);
	float4v sin(const float4v& x);
	void    sincos(const float4v& x, float4v& sin, float4v& cos);
	float4v sindg(const float4v& x);
	float4v tan(const float4v& x);
	float4v tandg(const float4v& x);
	float4v cot(const float4v& x);
	float4v cotdg(const float4v& x);

	// hyperbolic functions
	float4v cosh(const float4v& x);
	float4v sinh(const float4v& x);
	float4v tanh(const float4v& x);
	float4v acosh(const float4v& x);																										// arc hyperbolic cosine
	float4v asinh(const float4v& x);																										// arc hyperbolic sine
	float4v atanh(const float4v& x);																										// arc hyperbolic tangent

	float4v sinc(const float4v& x);
	float4v sinhc(const float4v& x);

	// exponential and logarithms
	float4v exp(const float4v& x);
	float4v exp10(const float4v& x);
	float4v exp2(const float4v& x);
	float4v expm1(const float4v& x);
	float4v expx2(const float4v& x);
	float4v log(const float4v& x);
	float4v log2(const float4v& x);
	float4v log10(const float4v& x);
	float4v log1p(const float4v& x);
	float4v logb(const float4v& x, float base);

	float4v pow(const float4v& x, const float4v& y);
	float4v pow(const float4v& x, int N);
	float4v pow(const float4v& x, const int4v& y);
	float4v pow_v2(const float4v& x, const float4v& y);
	int4v   pow(const int4v& x, int n);
	int4v   pow(const int4v& x, const int4v& n);
	float4v powm1(const float4v& x, const float4v& y);                                  // x^y - 1

	float4v cbrt(const float4v&);                                                       // cube root
	float4v nthrt(const float4v&, int n);                                               // nth root
	__forceinline float4v sqrt(const float4v& x) { return _mm_sqrt_ps(x.value()); }
	float4v sqrt1pm1(const float4v& x);                                                 // sqrt(1+x) - 1
	float4v hypot(const float4v& x, const float4v& y);

	mask4v isfinite(const float4v& x);
    
	// polynomial evaluations
	float4v chbevl(const float4v& x, const float coeff[], int N);			// Chebyshev polynom evaluation
	float4v polevl(const float4v& x, const float coeff[], int N);			// poly evaluation, 
	float4v p1evl(const float4v& x, const float coeff[], int N);			// poly evaluation, assume coeff[0] is 1.0 and omitted
	float4v legendre_p(int n, int m, const float4v& x);
	float4v legendre_q(int n, int m, const float4v& x);


	// Special functions
	void airy(const float4v& x, float4v& ai, float4v& aip, float4v& bi, float4v& bip);	// Airy function
	float4v beta(const float4v& a, const float4v& b);									// beta function
	float4v dawsn(const float4v& x);													// Dawson integral
	float4v fac(const int4v& i);														// factorial
	float4v gamma(const float4v& x);													// Gamma function
	float4v lgam(const float4v& x);														// log of Gamma function
	float4v psi(const float4v& x);														// psi (digamma) function, i.e. log derivative of log-gamma function
	float4v rgamma(const float4v& x);													// reciprocal Gamma function
	float4v spence(const float4v& x);													// dilogarithm function
	float4v zeta(const float4v& x, const float4v& q);									// Riemann Zeta function of 2 args
	float4v zetac(const float4v& x);													// Riemann Zeta function

	// Integrals
	float4v ellie(const float4v& phi, float m);											// incomplete elliptiptic integral of the second kind
	float4v ellik(const float4v& phi, float m);											// incomplete elliptiptic integral of the first kind
	float4v ellpe(const float4v& m1);													// complete elliptiptic integral of the second kind
	float4v ellpj(const float4v& u, const float4v& m, float4v& sn, float4v& cn, float4v& dn, float4v& phi); // jacobian elliptic functions
	float4v ellpk(const float4v& m1);													// complete elliptiptic integral of the first kind
	float4v ei(const float4v& x);
	float4v expn(int n, const float4v& x);												// exponential integral
	void    fresnl(const float4v& x, float4v& S, float4v& C);							// Fresnel integral
	// fresnel_integral
	float4v igam(const float4v& a, const float4v& x);									// incomplete gamma integral
	// gamma_integral_incomplete
	float4v igamc(const float4v& a, const float4v& x);									// complemented incomplete gamma integral
	float4v igami(const float4v& a, const float4v& y);									// inverse of complemented incomplete gamma integral
	float4v incbet(const float4v& a, const float4v& b, const float4v& x);               // incomplete beta integral
	float4v incbi(const float4v& a, const float4v& b, const float4v& y);                // inverse of incomplete beta integral
	// beta_integral_incomplete_inverse
	void    shichi(const float4v& x, float4v& Shi, float4v& Chi);						// Hyperbolic sine and cosine integrals
	// sinh_cosh_integral
	void    sici(const float4v& x, float4v& Si, float4v& Ci);						    // sine and cosine integrals
	// sin_cos_integral
	float4v struve(const float4v& v, const float4v& x);								    // Struve function

	// statistics
	float4v bdtr(const int4v& k, const int4v& n, const float4v& p);                     // binomial distribution
	float4v bdtrc(const int4v& k, const int4v& n, const float4v& p);                    // complemented binomial distribution
	float4v bdtri(const int4v& k, const int4v& n, const float4v& p);                    // inverse binomial distribution
	float4v chdtr(float df, const float4v& x);											// chi-square distribution
	float4v chdtrc(float v, const float4v& x);											// complemented chi-square distribution
	float4v chdtic(float df, const float4v& x);											// inverse of complemented chi-square distribution
	float4v fdtr(const int4v& df1, const int4v& df2, const float4v& x);                 // F distribution
	float4v fdtrc(const int4v& df1, const int4v& df2, const float4v& x);                // complemented F distribution
	float4v fdtri(const int4v& df1, const int4v& df2, const float4v& x);                // inverse of complemented F distribution
	float4v gdtr(const float4v& a, const float4v& b, const float4v& x);                 // Gamma distribution
	float4v gdtrc(const float4v& a, const float4v& b, const float4v& x);                // complemented Gamma distribution
	float4v hyp2f1(const float4v& a, const float4v& b, const float4v& c, const float4v& x);  // Gauss hypergeometric function
	float4v hyperg(const float4v& a, const float4v& b, const float4v& x);               // confluent hypergeometric expansion
	float4v nbdtr(const int4v& k, int n, const float4v& p);								// negative binomial distribution
	float4v nbdtrc(const int4v& k, int n, const float4v& p);							// complemented negative binomial distribution
	float4v ndtr(const float4v& x);								                        // normal distribution function
	float4v ndtri(const float4v& x);								                    // inverse of normal distribution function
	float4v pdtr(const int4v& k, const float4v& m);				                        // Poisson distribution
	float4v pdtrc(const int4v& k, const float4v& m);				                    // complemented Poisson distribution
	float4v pdtri(const int4v& k, const float4v& y);				                    // inverse Poisson distribution
	float4v stdtr(int k, const float4v& t);				                                // Student's distribution


	// error functions
	float4v erf(const float4v& x);								                        // error function
	//float4v erfi(const float4v& x);								                    //  inverse error function
	float4v erfc(const float4v& x);								                        // complimentary error function
	//float4v erfci(const float4v& x);							                        // inverse complimentary error function

	// Bessel functions
	float4v i0(const float4v& x);													    // Modified Bessel function of order zero
	float4v i0e(const float4v& x);														// Modified Bessel function of order zero, exponentially scaled
	float4v i1(const float4v& x);														// Modified Bessel function of order one
	float4v i1e(const float4v& x);														// Modified Bessel function of order one, exponentially scaled
	float4v iv(const float4v& v, const float4v& x);										// Modified Bessel function of non-integer order
	float4v j0(const float4v& x);														// Bessel function order zero
	float4v j1(const float4v& x);														// Bessel function order one
	float4v jn(int n, const float4v& x);												// Bessel function of integer order
	float4v jv(float v, const float4v& x);												// Bessel function of non-integer order
	float4v k0(const float4v& x);														// Bessel function 3rd kind, order zero
	float4v k0e(const float4v& x);														// Bessel function 3rd kind, order zero, exponentially scaled
	float4v k1(const float4v& x);														// Bessel function 3rd kind, order one
	float4v k1e(const float4v& x);														// Bessel function 3rd kind, order one, exp scaled
	float4v kn(int n, const float4v& x);												// Bessel function 3rd kind, integer order
	// kv?
	float4v y0(const float4v& x);														// Bessel function of 2nd kind, order zero
	float4v y1(const float4v& x);														// Bessel function of 2nd kind, order one
	float4v yn(int n, const float4v& x);												// Bessel function 2nd kind, integer order
	float4v yv(float v, const float4v& x);												// Bessel function of 2nd kind, non-integer order
	
    void cartesian2polar(const float4v& x, const float4v& y, float4v& r, float4v& phi);
    void polar2cartesian(const float4v& r, const float4v& phi, float4v& x, float4v& y);
    void cartesian2spherical(const float4v& x,const float4v& y,const float4v& z, float4v& r, float4v& phi, float4v& theta);
    void spherical2cartesian(const float4v& r,const float4v& phi,const float4v& theta, float4v& x, float4v& y, float4v& z);

    //! Calculate cross product of 2 vectors. Vectors described with elements 0..2, #3 is undefined.
    __forceinline
    float4v cross(const float4v& u, const float4v& v) {
        /*
        w.x = u.y*v.z - u.z*v.y;
        w.y = u.x*v.z - u.z*v.x;
        w.z = u.x*v.y - u.y*v.x;
        */
        float4v w = shuffle<1,0,0,0>(u,u)*shuffle<2,2,1,1>(v,v) 
                  - shuffle<2,2,1,1>(u,u)*shuffle<1,0,0,0>(v,v);
        return w;
    }


    template <typename T, size_t N>
    inline size_t array_sizeof(const T (&)[N]) { return  N; }


	void
	quadratic_roots(const float4v& a, const float4v& b, const float4v& c,
					float4v& r0, float4v& r1, int4v& what);
	void
	quadratic_real_roots(const float4v& a, const float4v& b, const float4v& c,
						 float4v& r0, float4v& r1);


	__forceinline float4v reciprocal(const float4v& x)                  { return _mm_rcp_ps(x.value()); }
	__forceinline float4v sqrt_reciprocal(const float4v& x)             { return _mm_rsqrt_ps(x.value()); }
	float4v reciprocal_nr(const float4v& x);
	float4v sqrt_reciprocal_nr(const float4v& x);




	// ordering
	__forceinline float4v min(const float4v& x, const float4v& y) { return _mm_min_ps(x.value(), y.value()); }
	__forceinline float4v min(const float4v& x, const float4v& y, const float4v& z) { return min(x, min(y, z)); }
	__forceinline float4v min(const float4v& x, const float4v& y, const float4v& z, const float4v& w) { return min(min(x, y), min(z, w)); }
	__forceinline float4v max(const float4v& x, const float4v& y) { return _mm_max_ps(x.value(), y.value()); }
	__forceinline float4v max(const float4v& x, const float4v& y, const float4v& z) { return max(x, max(y, z)); }
	__forceinline float4v max(const float4v& x, const float4v& y, const float4v& z, const float4v& w) { return max(max(x, y), max(z, w)); }
	__forceinline float4v clip(const float4v& x, const float4v& a, const float4v& b) { return min(max(x, a), b); }


	float4v	min_vector(const float4v& r);
	float4v	min_vector(const float4v& r0, float4v& r1);
	float4v	min_vector(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3);
	float4v	max_vector(const float4v& r);
	float4v	max_vector(const float4v& r0, float4v& r1);
	float4v	max_vector(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3);

    float4v sort_2_asc(const float4v& x);
    float4v sort_3_asc(const float4v& x);
    float4v sort_4_asc(const float4v& x);
    float4v sort_2_desc(const float4v& x);
    float4v sort_3_desc(const float4v& x);
    float4v sort_4_desc(const float4v& x);

	//! Compute y = { sum(r), x, x, x }
	__forceinline float4v
	sum_vector(const float4v& r) {
		return sum_row(r);
	}

	//! Compute y = { sum(r0), sum(r1), x, x }
	__forceinline float4v
	sum_vector(const float4v& r0, const float4v& r1) {
		return sum_row(r0, r1);
	}


#if 0
	template <typename T>
	struct minimum : public std::binary_function<T, T, T> {
		result_type operator(const T& lhs, const T& rhs) const { return min(lhs, rhs); }
	};

	template <typename T>
	struct maximum : public std::binary_function<T, T, T> {
		result_type operator(const T& lhs, const T& rhs) const { return max(lhs, rhs); }
	};

	template <typename T>
	struct plus : public std::binary_function<T, T, T> {
		result_type operator(const T& lhs, const T& rhs) const { return lhs + rhs; }
	};


	//! Compute y = { min(r), x, x, x }
	template <typename Op, typename T>
	__forceinline T
	acc_vector(const float4v& r, const Op& operation) {
	//! Compute y = { min(r), x, x, x }
    T tmp = operation(r, _mm_movehl_ps(r.value(), r.value()));    // { r_02, r_13, x, x} ignore upper 2 floats
    T tmp2 = shuffle<1, 0, 0, 1>(tmp, tmp); //_mm_set_ss(tmp.m128_f32[1]);
    return operation(tmp.value(), tmp2.value());
	}

	//! Compute y = { min(r0), min(r1), x, x }
	template <typename Op, typename T>
	__forceinline T
	acc_vector(const T& r0, const T& r1, const Op& operation) {
    T ab = operation(shuffle<0,1,0,1>(r0, r1), shuffle<2,3,2,3>(r0, r1));  // {r0_02, r0_13, r1_02, r1_13}
    T vmin = operation(shuffle<0,2,0,2>(ab, ab), shuffle<1,3,1,3>(ab, ab));
    return vmin;
	}

	//! Compute y = { min(r0), min(r1), min(r2), min(r3) }
	template <typename Op, typename T>
	__forceinline T
	acc_vector(const T& r0, const T& r1, const T& r2, const T& r3, const Op& operation) {
    T ab = operation(shuffle<0,1,0,1>(r0, r1), shuffle<2,3,2,3>(r0, r1));  // {r0_02, r0_13, r1_02, r1_13}
    T cd = operation(shuffle<0,1,0,1>(r2, r3), shuffle<2,3,2,3>(r2, r3));  // {r2_02, r2_13, r3_02, r3_13}

    T vmin = operation(shuffle<0,2,0,2>(ab, cd), shuffle<1,3,1,3>(ab, cd));
    return vmin;
	}


		//! Compute y = { min(r), x, x, x }
	__forceinline float4v
	min_vector(const float4v& r) {
		return acc_vector(r, minimum<float4v>());
	}

	//! Compute y = { min(r0), min(r1), x, x }
	__forceinline float4v
	min_vector(const float4v& r0, float4v& r1) {
		return acc_vector(r0, r1, minimum<float4v>());
	}

	//! Compute y = { min(r0), min(r1), min(r2), min(r3) }
	__forceinline float4v
	min_vector(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3) {
		return acc_vector(r0, r1, r2, r3, minimum<float4v>());
	}


	//! Compute y = { max(r), x, x, x }
	__forceinline float4v
	max_vector(const float4v& r) {
		return acc_vector(r, maximum<float4v>());
	}

	//! Compute y = { max(r0), max(r1), x, x }
	__forceinline float4v
	max_vector(const float4v& r0, const float4v& r1) {
		return acc_vector(r0, r1, maximum<float4v>());
	}

	//! Compute y = { max(r0), max(r1), max(r2), max(r3) }
	__forceinline float4v
	max_vector(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3) {
		return acc_vector(r0, r1, r2, r3, maximum<float4v>());
	}
	
	//! Compute y = { sum(r0), sum(r1), sum(r2), sum(r3) }
	//! Compute y = { max(r), x, x, x }
	__forceinline float4v
	sum_vector(const float4v& r) {
		return acc_vector(r, plus<float4v>());
	}

	//! Compute y = { max(r0), max(r1), x, x }
	__forceinline float4v
	sum_vector(const float4v& r0, const float4v& r1) {
		return acc_vector(r0, r1, plus<float4v>());
	}

	//! Compute y = { max(r0), max(r1), max(r2), max(r3) }
	__forceinline float4v
	sum_vector(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3) {
		return acc_vector(r0, r1, r2, r3, plus<float4v>());
	}
#endif



	// normalization
	//float4v ceil(const float4v& x);
	//float4v floor(const float4v& x);
	//float4v truncate(const float4v& x);
	//float4v nearest(const float4v& x);
		__forceinline float4v abs(const float4v& x) {
		// return x & float4v::valuemask()
		return float4v::create_reinterpreted(x.reinterpret_bits() & float4v::valuemask());
	}

	__forceinline int4v   signum(const float4v& x) { 
		int4v  result = (int4v(x < 0) & (-1)) | (int4v(x > 0) & 1);
		return result;
	}
	
	__forceinline int4v   signum(const float4v& x, const float4v& lo, const float4v& hi) {
		int4v  result = (int4v(x < lo) & (-1)) | (int4v(x > hi) & 1);
		return result;
	}


  float4v fmod(const float4v& x, const float4v& y);
  float4v modf(const float4v& x, float4v& int_part);

	__forceinline float4v floor(const float4v& x)							{ return _mm_floor_ps(x.value()); }
	__forceinline float4v ceil(const float4v& x)							{ return _mm_ceil_ps(x.value()); }
	__forceinline float4v round(const float4v& x)							{ return _mm_round_ps(x.value(), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_RAISE_EXC); }
	__forceinline float4v trunc(const float4v& x)							{ return _mm_round_ps(x.value(), _MM_FROUND_TO_ZERO | _MM_FROUND_RAISE_EXC); }
	__forceinline int4v   round_int(const float4v& x)						{ return _mm_cvtps_epi32(x.value()); }
	__forceinline int4v   trunc_int(const float4v& x)						{ return _mm_cvttps_epi32(x.value()); }
	template <size_t N> __forceinline int trunc_int(const float4v& x) 		{ return _mm_cvttss_si32(float4v::copy_to_0<N>(x).value()); }
	template <size_t N> __forceinline int round_int(const float4v& x) 		{ return _mm_cvtss_si32(float4v::copy_to_0<N>(x).value()); }
	template <size_t N> __forceinline float ceil(const float4v& x) 			{ return _mm_ceil_ss(x.value(), float4v::copy_to_0<N>(x).value()); }
	template <size_t N> __forceinline float floor(const float4v& x)			{ return _mm_floor_ss(x.value(), float4v::copy_to_0<N>(x).value()); }
	template <size_t N> __forceinline float round(const float4v& x)			{ return _mm_round_ss(x.value(), float4v::copy_to_0<N>(x).value(), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_RAISE_EXC); }
	template <size_t N> __forceinline float trunc(const float4v& x)			{ return _mm_round_ss(x.value(), float4v::copy_to_0<N>(x).value(), _MM_FROUND_TO_ZERO | _MM_FROUND_RAISE_EXC); }


	float4v frexp(const float4v& x, int4v& expnt);
	float4v ldexp(const float4v& x, const int4v& expnt);
	mask4v signbit(const float4v& x);
	__forceinline mask4v isnan(const float4v& x) {
			return int4v(abs(x).reinterpret_bits()) > float4v::infinites().reinterpret_bits();
	}

	__forceinline mask4v isinf(const float4v& x) {
		return abs(x).reinterpret_bits() == float4v::infinites().reinterpret_bits();
	}

	__forceinline mask4v isfinite(const float4v& x) {
		return int4v(abs(x).reinterpret_bits()) < float4v::infinites().reinterpret_bits();
	}

	__forceinline mask4v isnormal(const float4v& x) {
		int4v xbits = abs(x).reinterpret_bits();
		return xbits > 0x007fffff  && xbits < 0x7f800000;  // 0x007fffff is the largest subnormal number 
	}

	__forceinline bool any_nan(const float4v& x)				{ return isnan(x).any_true(); }
	__forceinline bool any_infinite(const float4v& x)			{ return isinf(x).any_true(); }
	__forceinline bool any_finite(const float4v& x)				{ return isfinite(x).any_true(); }

	__forceinline mask4v signbits(const float4v& x) {
		return _mm_and_ps(x.value(), float4v::signmask().value());  // return the signbits of x
	}

	__forceinline float4v copysign_unsafe(const float4v& x, const mask4v& signbits) { 
		// return abs(x) | signbits
		// return blend((abs(x), signbits)
		return float4v::create_reinterpreted(abs(x).reinterpret_bits() | signbits); 
	}
	
	__forceinline float4v copysign(const float4v& x, const mask4v& signbits) { 
		return copysign_unsafe(x, signbits & float4v::signmask()); 
	}
	
	__forceinline float4v copysign(const float4v& x, const int4v& signbits) { 
		return copysign(x, mask4v(signbits.value())); 
	}
	
	__forceinline float4v copysign(const float4v& x, const float4v& signbits) { 
		return copysign(x, signbits.reinterpret_bits());
	}
		
	__forceinline float4v changesign_unsafe(const float4v& x, const mask4v& signbits) { 
		// return x ^ signbits
		return float4v::create_reinterpreted(x.reinterpret_bits() ^ signbits); 
	}
	
	__forceinline float4v changesign(const float4v& x, const mask4v& signbits) { 
		return changesign_unsafe(x, signbits & float4v::signmask()); 
	}

	__forceinline float4v changesign(const float4v& x, const int4v& signbits) { 
		return changesign(x, mask4v(signbits.value())); 
	}


	__forceinline float4v changesign(const float4v& x, const float4v& signbits) { 
    return changesign(x, signbits.reinterpret_bits());
	}

	// geometry
	float4v row_matrix_multiply(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3, const float4v& x);
	float4v col_matrix_multiply(const float4v& c0, const float4v& c1, const float4v& c2, const float4v& c3, const float4v& x);


  /* NewtonRaphson Reciprocal
	[2 * rcpps(x) - (x * rcpps(x) * rcpps(x))] */
	__forceinline float4v reciprocal_nr(const float4v& x) {
		float4v Ra0 = reciprocal(x);
		return Ra0+Ra0 - (x*Ra0*Ra0);
		// is this really faster than the exact 1.0/x?  YES
		// i.e. return float4v(1.0f)/x?
	}

	/*	NewtonRaphson Reciprocal Square Root
	  0.5 * rsqrtps * (3 - x * rsqrtps(x) * rsqrtps(x)) */
	__forceinline float4v sqrt_reciprocal_nr(const float4v& x) {
		float4v Ra0 = sqrt_reciprocal(x);
		return (float4v(0.5f) * Ra0) * (float4v(3.0f) - (x * Ra0) * Ra0);
		// is this really faster than the exact 1.0/sqrt(x)?  YES
		// i.e. return float4v(1.0f)/x.sqrt()?
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	//
	// Implementations
	//
	//////////////////////////////////////////////////////////////////////////////////////////////////



	//! Compute y = { min(r), x, x, x }
	__forceinline float4v
	min_vector(const float4v& r) {
		float4v tmp = min(r, _mm_movehl_ps(r.value(), r.value()));    // { r_02, r_13, x, x} ignore upper 2 floats
		float4v tmp2 = shuffle<1, 0, 0, 1>(tmp, tmp); //_mm_set_ss(tmp.m128_f32[1]);
		return min(tmp.value(), tmp2.value());
	}

	//! Compute y = { min(r0), min(r1), x, x }
	__forceinline float4v
	min_vector(const float4v& r0, float4v& r1) {
		float4v ab = min(shuffle<0,1,0,1>(r0, r1), shuffle<2,3,2,3>(r0, r1));  // {r0_02, r0_13, r1_02, r1_13}
		float4v vmin = min(shuffle<0,2,0,2>(ab, ab), shuffle<1,3,1,3>(ab, ab));
		return vmin;
	}

	//! Compute y = { min(r0), min(r1), min(r2), min(r3) }
	__forceinline float4v
	min_vector(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3) {
		float4v ab = min(shuffle<0,1,0,1>(r0, r1), shuffle<2,3,2,3>(r0, r1));  // {r0_02, r0_13, r1_02, r1_13}
		float4v cd = min(shuffle<0,1,0,1>(r2, r3), shuffle<2,3,2,3>(r2, r3));  // {r2_02, r2_13, r3_02, r3_13}

		float4v vmin = min(shuffle<0,2,0,2>(ab, cd), shuffle<1,3,1,3>(ab, cd));
		return vmin;
	}


	//! Compute y = { max(r), x, x, x }
	__forceinline float4v
	max_vector(const float4v& r) {
		float4v tmp = max(r, _mm_movehl_ps(r.value(), r.value()));    // { r_02, r_13, x, x} ignore upper 2 floats
		float4v tmp2 = shuffle<1, 0, 0, 1>(tmp, tmp); //_mm_set_ss(tmp.m128_f32[1]);
		return max(tmp.value(), tmp2.value());
	}

	//! Compute y = { max(r0), max(r1), x, x }
	__forceinline float4v
	max_vector(const float4v& r0, const float4v& r1) {
		float4v ab = max(shuffle<0,1,0,1>(r0, r1), shuffle<2,3,2,3>(r0, r1));  // {r0_02, r0_13, r1_02, r1_13}
		float4v vmax = max(shuffle<0,2,0,2>(ab, ab), shuffle<1,3,1,3>(ab, ab));
		return vmax;
	}

	//! Compute y = { max(r0), max(r1), max(r2), max(r3) }
	__forceinline float4v
	max_vector(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3) {
		float4v ab = max(shuffle<0,1,0,1>(r0, r1), shuffle<2,3,2,3>(r0, r1));  // {r0_02, r0_13, r1_02, r1_13}
		float4v cd = max(shuffle<0,1,0,1>(r2, r3), shuffle<2,3,2,3>(r2, r3));  // {r2_02, r2_13, r3_02, r3_13}

		float4v vmax = max(shuffle<0,2,0,2>(ab, cd), shuffle<1,3,1,3>(ab, cd));
		return vmax;
	}

	// sum_vector
	// and_vector
	// or_vector


	namespace intern {

		// y = c0 + x*(c1 + x*( ... + x*(c{N-1} + c{N}*x)...)
		template <size_t N>
		struct eval_poly_impl {
			template <typename U, typename T>
			static __forceinline auto 
			eval(const T& x, const U *c) -> decltype(*c+x) {
				return x * eval_poly_impl<N-1>::eval(x, c+1) + c[0];
			}
		};

		template <>
		struct eval_poly_impl<0> {
			template <typename U, typename T>
#pragma warning(push)
#pragma warning(disable:4100)  // x unreferenced formal parameter
			static __forceinline auto eval(const T& x, const U *c) -> decltype(*c+x) { return c[0];	}
#pragma warning(pop)
		};

		template <size_t N>
		struct eval_1poly_impl {
			template <typename U, typename T>
			static __forceinline auto 
			eval(const T& x, const U *c) -> decltype(*c+x) {
				return x * eval_1poly_impl<N-1>::eval(x, c+1) + c[0];
			}
		};

#if 0
		// y = (((c0*x + c1)*x + ... -c{N-1})*x + c{N}
		template <size_t N>
		struct eval_poly_reverse_impl {
			template <typename U, typename T>
			static __forceinline auto 
			eval(const U *c, const T& x, acc) -> decltype(*c+x) {
				acc = (acc + c[0])*x;
				return eval_poly_reverse_impl<N-1>::eval(c+1, x, acc);
			}
		};
#endif

		template <>
		struct eval_1poly_impl<1> {
			template <typename U, typename T>
			static __forceinline auto eval(const T& x, const U *c) -> decltype(c[0] + x) { x + c[0]; }
		};

		template <>
		struct eval_1poly_impl<0> {
			template <typename U, typename T>
			static __forceinline auto eval(const T& x, const U *) -> decltype(declval<U>() +x) { 1;	}
		};


    template <size_t N>
    struct pow {
      template <typename T>
      static __forceinline T eval(const T& x) {
        return (N & 1) != 0 ? x*pow<N/2>::eval(x*x) : pow<N/2>::eval(x*x);
      }
    };

    template <>
    struct pow<1> {
      template <typename T>
      static __forceinline T eval(const T& x) {
        return x; 
      }
    };

    template <>
    struct pow<0> {
      template <typename T>
      static __forceinline T eval(const T& ) {
        return 1; 
      }
    };

	}
	// polynomials

  /*! Polynomial evaluation by Horners rule.  Polynomial order is N-1.
				y = c0 + x*(c1 + x*(c1 + ...x*(c{N-1} + x*c{N})...)
	*/
  template <typename U, typename T, size_t N>
  __forceinline auto 
	eval_polynomial(const T& x, const U (&c)[N]) -> decltype(intern::eval_poly_impl<N>::eval(x, c)) {
		return intern::eval_poly_impl<N-1>::eval(x, c);
	}

	/*! Polynomial evaluation by Horners rule.  Polynomial order is N.
	This is an optimization of eval_polynomial where the coefficient of x^(N-1) is 1.0.
	Thus this constant is omitted from c.  i.e. 
		y = c0 + x*(c1 + x*(c1 + ...x*(c{N-1} + x)...)


	N=2:  y = c0 + x*(c1 + x)
	N=1:  y = c0 + x
	N=0:  y = 1
	*/
  template <typename U, typename T, size_t N>
  __forceinline auto 
	eval_1polynomial(const T& x, const U (&c)[N]) -> decltype(intern::eval_1poly_impl(x, c)) {
		return intern::eval_1poly_impl(x, c);
	}

  /*! Chebychev polynomial evaluation by Horners rule.  Polynomial order is N-1.
				y = c0 + x*(c1 + x*(c1 + ...x*(c{N-1} + x*c{N})...)
	*/
  template <typename U, typename T, size_t N>
  __forceinline T 
	  eval_chebychev(const T& x, const U (&c)[N]) {
          return x + c[0];  // TODO implement this
  }

  template <size_t N, typename T>
  __forceinline T pow(const T& x) {
    return intern::pow<N>::eval(x);
  }

}}}}

	
/*
----  void airyf( float, float *, float *, float *, float * );  LEVEL = 13
c---  float acosf( float );
c---  float acoshf( float );
c---  float asinf( float );
c---  float asinhf( float );
c---  float atanf( float );
c---  float atanhf( float );
c---  float atan2f( float, float );
c---  float bdtrcf( int, int, float );                          LEVEL = 2
c---  float bdtrf( int, int, float );                           LEVEL = 3
c---  float bdtrif( int, int, float );                          LEVEL = 3
c---  float betaf( float, float );                              LEVEL = 5
c---  float cbrtf( float );
----  float chbevlf( float, float *, int );                     --> eval_chebychev
c---  float chdtrcf(float, float);                              done LEVEL = 1
c---  float chdtrf(float, float);                               done LEVEL = 1
c---  float chdtrif( float, float );                            done LEVEL = 1
ctp-  float ceilf( float );
----  float cosdgf( float );                                    LEVEL = 5
ct--  float cosf(float);
c---  float coshf(float);
c---  float cot( float );
----  float cotdgf( float );
c---  float dawsnf(float);                                      LEVEL = 5
x---- float eif(const float x);
c---  float ellief( float, float );
c---  float ellikf( float, float );
c---  float ellpef(float);
c---  float ellpjf( float, float, float *, float *, float *, float * );
c---  float ellpkf(float);
c---  float erff(float);
c---  float erfcf(float);
c---  float erfif(float);
c---  float erfcif(float);
c---  float expf(float);
c---  float expf2(float);
c---  float exp10f(float);
----  float expnf( int, float );                                LEVEL = 13
c---  float facf( int );                                        done LEVEL = 1
c---  float fdtrcf( int, int, float );                          done LEVEL = 2
c---  float fdtrf( int, int, int );                             done LEVEL = 2
c---  float fdtrif( int, int, int );                            done LEVEL = 2
ctp-  float floorf(float);
c---  void fresnlf( float, float *, float * );                  done LEVEL=3
ct--  float frexpf(float, int *);
----  float gammaf(float);                                      LEVEL = 21 WIP
c---  float gdtrf( float, float, float );                       done LEVEL = 1
c---  float gdtrcf( float, float, float );                      done LEVEL = 1
----  NEXT float hyp2f1f( float, float, float, float );         LEVEL = 8
----  NEXT float hyp2f0f(float, float, float, int, float *);    LEVEL = 5
----  NEXT float hypergf( float, float, float );                LEVEL = 3
c---  float i0f( float );
c---  float i0ef( float );
c---  float i1f( float );
c---  float i1ef( float );
c---  float igamcf(float, float);                               LEVEL = 8
c---  float igamf(float, float);                                LEVEL = 5
c---  float igamif(float, float);                               LEVEL = 8
----  float incbetf(float, float, float);                       LEVEL = 34
----  float incbcff(float, float, float);                       LEVEL = 8
----  float incbdf(float, float, float);                        LEVEL = 8
----  float incbif( float, float, float );                      LEVEL = 21
----  float ivf( float, float );                                LEVEL = 5
c---  float j0f( float );
c---  float j1f( float );                                       LEVEL = 5
----  float jnf( int, float );                                  LEVEL = 5
----  float jvf( float, float );                                LEVEL = BIG
c---  float k0f( float );                                       done LEVEL = 1
c---  float k0ef( float );                                      done LEVEL = 1
c---  float k1f( float );                                       done LEVEL = 1
c---  float k1ef( float );                                      done LEVEL = 2
----  float knf( int, float );                                  LEVEL = 13
ct--  float ldexpf(float, int);
----  float lgamf(float);                                       LEVEL = 21
c---  float logf( float );
c---  float log2f( float );
c---  float log10f( float );
----	float4v log1p(const float4v& x);                            ??? need to find power series
c---	float4v logb(const float4v& x, float base);
c---  float nbdtrcf( int, int, float );                         done LEVEL = 1
c---  float nbdtrf( int, int, float );                          done LEVEL = 1
c---  float ndtrf( float );
c---  float ndtrif( float );                                    LEVEL = 5  similar to ndtrf
----  float nthrt( float, int n );                              see cbrt
----  float onef2f( float, float, float, float, float * );      LEVEL = 8 struve
c---  float pdtrcf( int, float );                               done LEVEL = 1
c---  float pdtrf( int, float );                                done LEVEL = 1
c---  float pdtrif( int, float );                               done LEVEL = 1
c---  float polevlf( float, float *, int );                     --> eval_polynomial
x---  float p1evlf( float, float *, int );                      --> eval_polynomial (special case omitted)
c-p-  float powf(float, float);
c-p-  float powif(float, int);
----  float psif( float );                                      LEVEL = 13
c---  float rgammaf( float );                                   done LEVEL = 8
----  int shichif( float, float *, float * );                   LEVEL = 8
c---  int sicif( float, float *, float * );                     LEVEL = 8
----  float sindgf( float );                                    LEVEL = 5
ctp-  float sinf( float );
ctp-  float sincos( float );
c---  float sinhf( float );
----  float spencef( float );                                   LEVEL = 8
ct--  float sqrtf( float );
----  float stdtrf( int, float );                               LEVEL = 13
----  float struvef( float, float );                            LEVEL = 5
----  float tandgf( float );
c---  float tanf( float );
c---  float tanhf( float );
----  float threef0f( float, float, float, float, float * );    LEVEL = 8 struve
c---  float y1f( float );                                       done LEVEL = 5
----  float ynf( int, float );                                  LEVEL = 8
c---  float yvf( float, float );                                done LEVEL = 1 struve
c---  float zetacf( float );
c---  float zetaf( float, float );

http://code.google.com/p/project-eneida/source/browse/util/util_amath_c.inc?r=95da92bcfed3070df60b0912c245fee21c81b074

*/