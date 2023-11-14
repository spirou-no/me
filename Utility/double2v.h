#pragma once

#include <utility/float4v.h>
#include <utility/int4v.h>
#include <utility/mask2v.h>
#include <utility/sse_internal.h>
#include <xmmintrin.h>   // SSE
#include <emmintrin.h>   // SSE2
#include <pmmintrin.h>   // SSE3
#include <smmintrin.h>   // SSE4



namespace oyrke { namespace algorithm { namespace sse {

  class double2v {
  private:
    __m128d value_;

  public:
    static __forceinline double2v set_all(double x)                      { return double2v(x); }
    static __forceinline double2v create_aligned(const double data[2])   { return double2v(data, aligned); }
    static __forceinline double2v create_unaligned(const double data[2]) { return double2v(data, unaligned); }
    static __forceinline double2v zeros()                                { return _mm_setzero_pd(); }
    static __forceinline double2v compose_01(const double2v& a, const double2v& b) { return _mm_move_sd(a.value(), b.value()); }
    static __forceinline double2v compose_00(const double2v& a, const double2v& b) { return _mm_unpacklo_pd(a.value(), b.value()); }
    static __forceinline double2v compose_11(const double2v& a, const double2v& b) { return _mm_unpackhi_pd(a.value(), b.value()); }

    __forceinline double2v()                                             { value_ = zeros().value(); }
    __forceinline double2v(__m128d x)                                    { value_ = x; }
    __forceinline explicit double2v(double x)                            { value_ = _mm_set1_pd(x); }
    __forceinline          double2v(double x, double y)                  { value_ = _mm_set_pd(y,x); }
    __forceinline          double2v(const double data[2], UNALIGNED)     { value_ = _mm_loadu_pd(data); }
    __forceinline          double2v(const double data[2], ALIGNED)       { value_ = _mm_load_pd(data); }
    __forceinline explicit double2v(const int4v& x)                      { value_ = _mm_cvtepi32_pd(x.value()); }
    __forceinline explicit double2v(const float4v& x)                    { value_ = _mm_cvtps_pd(x.value()); }

    __forceinline __m128d value() const                                  { return value_; }
		template <size_t N> __forceinline __m128d copy_to_0() const          { return _mm_castsi128_pd(_mm_srli_si128(_mm_castpd_si128(value_), N*8)); } //_mm_shuffle_ps(value_, value_, N); }
		template <>         __forceinline __m128d copy_to_0<0>() const       { return value_; }


		template <size_t N> double at() const																 { return _mm_cvtsd_f64(copy_to_0<N>()); }
		template <size_t N> __forceinline void set_at(double x)							 { value_ = _mm_blend_pd(value_, double2v(x).value_, 1<<N); }
    
		__forceinline int4v trunc_int() const																 { return _mm_cvttpd_epi32(value_); }
		__forceinline int4v round_int() const																 { return _mm_cvtpd_epi32(value_); }
		__forceinline float4v convert_float() const                          { return _mm_cvtpd_ps(value_); }
		__forceinline double2v ceil() const																	 { return _mm_ceil_pd(value_); }
		__forceinline double2v floor() const															   { return _mm_floor_pd(value_); }
		__forceinline double2v round() const																 { return _mm_round_pd(value_, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_RAISE_EXC); }
		__forceinline double2v truncate() const															 { return _mm_round_pd(value_, _MM_FROUND_TO_ZERO | _MM_FROUND_RAISE_EXC); }

    template <size_t N> __forceinline int truncate_int() const           { return _mm_cvttsd_si32(copy_to_0<N>()); }
		template <size_t N> __forceinline int round_int() const							 { return _mm_cvtsd_epi32(copy_to_0<N>()); }
		template <size_t N> __forceinline double ceil() const								 { return _mm_cvtsd_f64(_mm_ceil_sd(value_, copy_to_0<0>())); }
		template <size_t N> __forceinline double floor() const							 { return _mm_cvtsd_f64(_mm_floor_sd(value_, copy_to_0<0>())); }
		template <size_t N> __forceinline double round() const							 { return _mm_cvtsd_f64(_mm_round_sd(value_, copy_to_0<0>(), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_RAISE_EXC)); }
		template <size_t N> __forceinline double trunc() const							 { return _mm_cvtsd_f64(_mm_round_sd(value_, copy_to_0<0>(), _MM_FROUND_TO_ZERO | _MM_FROUND_RAISE_EXC)); }

    __forceinline void store(double *ptr, ALIGNED)                       { _mm_store_pd(ptr, value_);}
    __forceinline void store(double *ptr, UNALIGNED)                     { _mm_storeu_pd(ptr, value_);}
		template <size_t N> __forceinline void store(double& dest)           { _mm_store_sd(&dest, copy_to_0<N>()); }

    __forceinline double  operator[](size_t i) const                     { return value_.m128d_f64[i]; }
    __forceinline double& operator[](size_t i)                           { return value_.m128d_f64[i]; }

    __forceinline const double2v& operator+=(const double2v& rhs)        { value_ = (*this + rhs).value_; return *this; }
    __forceinline const double2v& operator*=(const double2v& rhs)        { value_ = (*this * rhs).value_; return *this; }
    __forceinline const double2v& operator-=(const double2v& rhs)        { value_ = (*this - rhs).value_; return *this; }
    __forceinline const double2v& operator/=(const double2v& rhs)        { value_ = (*this / rhs).value_; return *this; }

    __forceinline double2v operator+(const double2v& rhs) const          { return _mm_add_pd(value_, rhs.value_); }
    __forceinline double2v operator+(double rhs) const                   { return *this + double2v(rhs); }

    __forceinline double2v operator-(const double2v& rhs) const          { return _mm_sub_pd(value_, rhs.value_); }
    __forceinline double2v operator-(double rhs) const                   { return *this - double2v(rhs); }

    __forceinline double2v operator*(const double2v& rhs) const          { return _mm_mul_pd(value_, rhs.value_); }
    __forceinline double2v operator*(double rhs) const                   { return *this * double2v(rhs); }

    __forceinline double2v operator/(const double2v& rhs) const          { return _mm_div_pd(value_, rhs.value_); }
    __forceinline double2v operator/(double rhs) const                   { return *this / double2v(rhs); }

    __forceinline mask2v   operator==(const double2v& rhs) const         { return _mm_cmpeq_pd(value_, rhs.value_); }
    __forceinline mask2v   operator!=(const double2v& rhs) const         { return _mm_cmpneq_pd(value_,rhs.value_); }
    __forceinline mask2v   operator< (const double2v& rhs) const         { return _mm_cmplt_pd(value_, rhs.value_); }
    __forceinline mask2v   operator<=(const double2v& rhs) const         { return _mm_cmple_pd(value_, rhs.value_); }
    __forceinline mask2v   operator> (const double2v& rhs) const         { return _mm_cmpgt_pd(value_, rhs.value_); }
    __forceinline mask2v   operator>=(const double2v& rhs) const         { return _mm_cmpge_pd(value_, rhs.value_); }

    __forceinline double2v sqrt() const                                  { return _mm_sqrt_pd(value_); }
    //__forceinline double2v reciprocal() const                            { return _mm_rcp_pd(value_); }
    //__forceinline double2v sqrt_reciprocal() const                       { return _mm_rsqrt_pd(value_); }
    __forceinline double2v min(const double2v& rhs) const                { return _mm_min_pd(value_, rhs.value_); }
    __forceinline double2v max(const double2v& rhs) const                { return _mm_max_pd(value_, rhs.value_); }
  };


  __forceinline double2v operator+(double lhs, const double2v& rhs)  { return rhs + lhs; }
  __forceinline double2v operator-(double lhs, const double2v& rhs)  { return double2v(lhs) - rhs; }
  __forceinline double2v operator*(double lhs, const double2v& rhs)  { return rhs * lhs; }
  __forceinline double2v operator/(double lhs, const double2v& rhs)  { return double2v(lhs) / rhs; }
  __forceinline double2v operator+(__m128d lhs, const double2v& rhs) { return rhs + lhs; }
  __forceinline double2v operator-(__m128d lhs, const double2v& rhs) { return double2v(lhs) - rhs; }
  __forceinline double2v operator*(__m128d lhs, const double2v& rhs) { return rhs * lhs; }
  __forceinline double2v operator/(__m128d lhs, const double2v& rhs) { return double2v(lhs) / rhs; }



  

  //! Merges x and y, selecting 2 floats from x and 2 from y.
  // x[x0], x[x1], y[y0]m y[y1]
  template <size_t x0, size_t y0>
  __forceinline double2v
  shuffle(const double2v& x, const double2v& y) {
    return _mm_shuffle_pd(x.value(), y.value(), _MM_SHUFFLE2(y0, x0));
  }

	//! Return 2-component dot product
	__forceinline double2v dot2(const double2v& lhs, const double2v& rhs) {
		return double2v(_mm_dp_pd(lhs.value(), rhs.value(), 0x31));
	}

	__forceinline double2v dot2(const double2v& x, const double2v& y0, const double2v& y1) {
		__m128d p = _mm_dp_pd(x.value(), y0.value(), 0x31);
		__m128d q = _mm_dp_pd(x.value(), y1.value(), 0x32);

		__m128d result = _mm_or_pd(p, q);
		return result;
	}



  //! Compute the sum of the 4 elements in r0 and r1.  Returned result is [sum(r0), sum(r1)].
  __forceinline double2v
  sum_row(const double2v &r0, const double2v& r1) {
#if 1
		return double2v(_mm_hadd_pd(r0.value(), r1.value()));
#else
		return double2v(_mm_unpackhi_pd(r0.value(), r1.value())) + _mm_unpacklo_pd(r0.value(), r1.value());
#endif
  }

  //! Compute the sum of the 4 elements in r.  Returned result is [sum(r), x], x i undefined
  __forceinline double2v
  sum_row(const double2v &r) {
#if 1
		return double2v(_mm_hadd_pd(r.value(), r.value()));
#else
		return r + _mm_unpackhi_pd(r.value(), r.value());
#endif
  }

  
}}}