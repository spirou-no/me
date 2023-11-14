#pragma once

#include <utility/int4v.h>
#include <utility/mask4v.h>
#include <utility/sse_internal.h>
#include <xmmintrin.h>   // SSE
#include <emmintrin.h>   // SSE2
#include <pmmintrin.h>   // SSE3
#include <smmintrin.h>   // SSE4
#include <immintrin.h>   // AVX


namespace oyrke { namespace algorithm { namespace sse {
  
  class int8v {
    __m256i value_;
  public:
    int8v(__m256i i)  : value_(i) {}
    __m256i value() const { return value_; }
  };

  class mask8v {
    __m256 value_;
  public:
    mask8v(__m256 m)  : value_(m) {}
    mask8v(unsigned int m)  : value_(_mm256_castsi256_ps(_mm256_set1_epi32(m))) {}
    __m256 value() const { return value_; }
  };

	class float8v {
  private:
    __m256 value_;

  public:
		__forceinline float8v()                                             { value_ = zeros().value(); }
    __forceinline float8v(__m256 x)                                     { value_ = x; }
    __forceinline float8v(float x)                                      { value_ = _mm256_set1_ps(x); }
    //__forceinline float8v(int x)                                      { value_ = _mm_set_ps1(float(x)); }
    __forceinline float8v(float a, float b, float c, float d,
                          float e, float f, float g, float h)           { value_ = _mm256_setr_ps(a, b, c, d, e, f, g, h); }
    __forceinline float8v(const float data[8], UNALIGNED)               { value_ = _mm256_loadu_ps(data); }
    __forceinline float8v(const float data[8], ALIGNED)                 { value_ = _mm256_load_ps(data); }
		__forceinline explicit float8v(const int8v& data)										{ value_ = _mm256_cvtepi32_ps(data.value()); }
		__forceinline float8v& operator=(const float8v& rhs)                { value_ = rhs.value_; return *this; }
		__forceinline float8v& operator=(const int8v& rhs)                  { value_ = float8v(rhs).value_; return *this; }
		__forceinline float8v& operator=(float rhs)													{ value_ = _mm_set_ps1(rhs); return *this; }

		static __forceinline float8v create0(float x)                       { return _mm256_set_ss(x); }
    static __forceinline float8v create_all(float x)                    { return float8v(x); }
    static __forceinline float8v create_aligned(const float data[4])    { return float8v(data, aligned); }
    static __forceinline float8v create_unaligned(const float data[4])  { return float8v(data, unaligned); }
		static __forceinline float8v create_reinterpreted(const int8v&  bitmask) { return _mm256_castsi256_ps(bitmask.value()); }
		static __forceinline float8v create_reinterpreted(const mask8v& bitmask) { return bitmask.value(); }
		static __forceinline float8v zeros()                                { return _mm256_setzero_ps(); }
		static __forceinline float8v ones()																	{ return float8v(1.0f); }
		static __forceinline float8v infinites()														{ __m256i m = _mm256_set1_epi32(0x7fe00000); return _mm256_castsi256_ps(m); }
		static __forceinline float8v NaNs()																	{ __m256i m = _mm256_set1_epi32(0x7fefffff); return _mm256_castsi256_ps(m); }
		static __forceinline mask8v  signmask()                             { return mask8v(0x80000000); }
		static __forceinline mask8v  valuemask()		                        { return mask8v(0x7fffffff); }

		template <size_t N> static __forceinline float8v copy_to_0(const float8v& x)    { return _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(x.value_), N*4)); } //_mm_shuffle_ps(value_, value_, N); }
		template <>         static __forceinline float8v copy_to_0<0>(const float8v& x) { return x.value_; }

    __forceinline void store(float *p, ALIGNED) const                   { _mm256_store_ps(p, value_); }
    __forceinline void store(float *p, UNALIGNED) const                 { _mm256_storeu_ps(p, value_); }
		template <size_t N> __forceinline void store(float& dest) const     { _mm_store_ss(&dest, copy_to_0<N>(*this).value_); } //_mm_store_ss(p, _mm_shuffle_ps(value_, value_, _MM_SHUFFLE(N,N,N,N))); }

    __forceinline __m256 value() const                                  { return value_; }
    __forceinline float  operator[](size_t i) const                     { return value_.m256_f32[i]; }
    __forceinline float& operator[](size_t i)                           { return value_.m256_f32[i]; }
		template <size_t N> __forceinline float at() const                  { return _mm256_cvtss_f32(copy_to_0<N>(*this).value_); } 
    
		template <size_t N> __forceinline void set_at(float x)              { value_ = _mm256_insert_ps(value_, float8v(x).value_, _MM_MK_INSERTPS_NDX(0,N,0)); } 
		template <size_t N> __forceinline void clear()                      { value_ = _mm256_insert_ps(value_, zeros().value_, N); }

		__forceinline mask8v reinterpret_bits() const                       { return mask8v(value_); }

		__forceinline const float8v& operator+=(const float8v& rhs)         { value_ = _mm256_add_ps(value_, rhs.value_); return *this; }
    __forceinline const float8v& operator*=(const float8v& rhs)         { value_ = _mm256_mul_ps(value_, rhs.value_); return *this; }
    __forceinline const float8v& operator-=(const float8v& rhs)         { value_ = _mm256_sub_ps(value_, rhs.value_); return *this; }
    __forceinline const float8v& operator/=(const float8v& rhs)         { value_ = _mm256_div_ps(value_, rhs.value_); return *this; }

		__forceinline bool equals(const float8v& rhs) const									{	return intern::is_all_set(_mm_cmpeq_ps(value_, rhs.value_)); }
  };

  __forceinline float8v operator+(const float8v& lhs, const float8v& rhs) { return _mm256_add_ps(lhs.value(), rhs.value()); }
  __forceinline float8v operator+(const float8v& lhs, float rhs)					{ return lhs + float8v(rhs); }
  __forceinline float8v operator+(float lhs, const float8v& rhs)          { return rhs + lhs; }

  __forceinline float8v operator-(const float8v& lhs, const float8v& rhs) { return _mm256_sub_ps(lhs.value(), rhs.value()); }
  __forceinline float8v operator-(const float8v& lhs, float rhs)					{ return lhs - float8v(rhs); }
  __forceinline float8v operator-(float lhs, const float8v& rhs)          { return float8v(lhs) - rhs; }
	__forceinline float8v operator-(const float8v& x)												{ return _mm256_xor_ps(x.value(), float8v::signmask().value()); }

  __forceinline float8v operator*(const float8v& lhs, const float8v& rhs) { return _mm256_mul_ps(lhs.value(), rhs.value()); }
  __forceinline float8v operator*(const float8v& lhs, float rhs)					{ return lhs * float8v(rhs); }
  __forceinline float8v operator*(float lhs, const float8v& rhs)          { return rhs * lhs; }

  __forceinline float8v operator/(const float8v& lhs, const float8v& rhs) { return _mm256_div_ps(lhs.value(), rhs.value()); }
  __forceinline float8v operator/(const float8v& lhs, float rhs)					{ return lhs / float8v(rhs); }
  __forceinline float8v operator/(float lhs, const float8v& rhs)          { return float8v(lhs) / rhs; }

  __forceinline mask8v  operator==(const float8v& lhs, const float8v& rhs) { return _mm256_cmp_ps(lhs.value(), rhs.value(), _CMP_EQ_OQ); }
  __forceinline mask8v  operator==(const float8v& lhs, float rhs)					 { return lhs == float8v(rhs); }
  __forceinline mask8v  operator==(float lhs, const float8v& rhs)          { return float8v(lhs) == rhs; }

  __forceinline mask8v  operator!=(const float8v& lhs, const float8v& rhs) { return _mm256_cmp_ps(lhs.value(), rhs.value(), _CMP_NEQ_OQ); }
  __forceinline mask8v  operator!=(const float8v& lhs, float rhs)					 { return lhs != float8v(rhs); }
  __forceinline mask8v  operator!=(float lhs, const float8v& rhs)          { return float8v(lhs) != rhs; }
		
  __forceinline mask8v  operator<(const float8v& lhs, const float8v& rhs)	 { return _mm256_cmp_ps(lhs.value(), rhs.value(), _CMP_LT_OQ); }
  __forceinline mask8v  operator<(const float8v& lhs, float rhs)					 { return lhs < float8v(rhs); }
  __forceinline mask8v  operator<(float lhs, const float8v& rhs)           { return float8v(lhs) < rhs; }

  __forceinline mask8v  operator<=(const float8v& lhs, const float8v& rhs) { return _mm256_cmp_ps(lhs.value(), rhs.value(), _CMP_LE_OQ); }
  __forceinline mask8v  operator<=(const float8v& lhs, float rhs)					 { return lhs <= float8v(rhs); }
  __forceinline mask8v  operator<=(float lhs, const float8v& rhs)          { return float8v(lhs) <= rhs; }

  __forceinline mask8v  operator>(const float8v& lhs, const float8v& rhs)  { return _mm256_cmp_ps(lhs.value(), rhs.value(), _CMP_GT_OQ); }
  __forceinline mask8v  operator>(const float8v& lhs, float rhs)					 { return lhs > float8v(rhs); }
  __forceinline mask8v  operator>(float lhs, const float8v& rhs)           { return float8v(lhs) > rhs; }

  __forceinline mask8v  operator>=(const float8v& lhs, const float8v& rhs) { return _mm256_cmp_ps(lhs.value(), rhs.value(), _CMP_GE_OQ); }
  __forceinline mask8v  operator>=(const float8v& lhs, float rhs)					 { return lhs >= float8v(rhs); }
  __forceinline mask8v  operator>=(float lhs, const float8v& rhs)          { return float8v(lhs) >= rhs; }

	
	__forceinline float8v blend(const float8v& x0, const float8v& x1) { 
		return _mm256_or_ps(x0.value(), x1.value());
	}

	__forceinline float8v blend(const float8v& x0, const float8v& x1, const float8v& x2) { 
		return blend(blend(x0, x1), x2);
	}

	__forceinline float8v blend(const float8v& x0, const float8v& x1, const float8v& x2, const float8v& x3) { 
		return blend(blend(x0, x1), blend(x2, x3));
	}

		/* Shift elements in lhs up N times, e.g. out[i]=x[i-N]; out[i]=0 if i<N */
		template <size_t N>
		__forceinline float8v shift_up(const float8v& x) {
			__m256i xi = _mm256_castps_si256(x.value());
			xi = _mm256_slli_si256(xi, N*sizeof(float)*8);
			return _mm256_castsi256_ps(xi);
		}

		/* Shift elements in lhs up N times, e.g. out[i]=x[i+N]; out[i]=0 if i>3-N */
		template <size_t N>
		__forceinline float8v shift_down(const float8v& x) {
			__m128i xi = _mm_castps_si128(x.value());
			xi = _mm_srli_si128(xi, N*sizeof(float)*8);
			return _mm_castsi128_ps(xi);
		}

		
		__forceinline float8v select(const mask8v& mask, const float8v& a, const float8v& b) {
			//			return _mm_or_ps(_mm_and_ps(mask.value(), a.value()), _mm_andnot_ps(mask.value(), b.value()));
			//return _mm_blendv_ps(b.value(), a.value(), mask.value());
			return intern::blend(mask.value(), a.value(), b.value());
		}

		__forceinline float8v select(const mask8v& mask, const float8v& a) {
			//			return _mm_or_ps(_mm_and_ps(mask.value(), a.value()), _mm_andnot_ps(mask.value(), b.value()));
			//return _mm_blendv_ps(b.value(), a.value(), mask.value());
			return _mm256_and_ps(mask.value(), a.value());
		}


		//! Return 4-component dot product
	__forceinline float8v dot4(const float8v& x, const float8v& y0, const float8v& y1, const float8v& y2, const float8v& y3) {
		__m256 p = _mm256_dp_ps(x.value(), y0.value(), 0xf1);
		__m256 q = _mm256_dp_ps(x.value(), y1.value(), 0xf2);
		__m256 r = _mm256_dp_ps(x.value(), y2.value(), 0xf4);
		__m256 s = _mm256_dp_ps(x.value(), y3.value(), 0xf8);

		__m256 result = _mm_or_ps(_mm_or_ps(p, q), _mm_or_ps(r, s));
		return result;
	}

	__forceinline float8v dot4(const float8v& x, const float8v& y0, const float8v& y1) {
		__m256 p = _mm256_dp_ps(x.value(), y0.value(), 0xf1);
		__m256 q = _mm256_dp_ps(x.value(), y1.value(), 0xf2);

		__m256 result = _mm_or_ps(p, q);
		return result;
	}

	namespace internal {
		template <size_t N>
		struct dotproduct_mask {
			enum mask { components = (0x0f >> (4-N)) << (4 + (4-N)); };
		};
	}

	//! Return 4-component dot product
	template <size_t N>
	__forceinline float8v dot(const float8v& x, const float8v& y0, const float8v& y1, const float8v& y2, const float8v& y3) {
		__m256 p = _mm256_dp_ps(x.value(), y0.value(), internal::dotproduct_mask<N>::components | 0x01);
		__m256 q = _mm256_dp_ps(x.value(), y1.value(), internal::dotproduct_mask<N>::components | 0x02);
		__m256 r = _mm256_dp_ps(x.value(), y2.value(), internal::dotproduct_mask<N>::components | 0x04);
		__m256 s = _mm256_dp_ps(x.value(), y3.value(), internal::dotproduct_mask<N>::components | 0x08);

		__m256 result = _mm256_or_ps(_mm_or_ps(p, q), _mm_or_ps(r, s));
		return result;
	}

	template <size_t N>
	__forceinline float8v dot(const float8v& x, const float8v& y0, const float8v& y1) {
		__m256 p = _mm256_dp_ps(x.value(), y0.value(), internal::dotproduct_mask<N>::components | 0x01);
		__m256 q = _mm256_dp_ps(x.value(), y1.value(), internal::dotproduct_mask<N>::components | 0x02);

		__m256 result = _mm256_or_ps(p, q);
		return result;
	}

	template <size_t N>
	//! Return N-component dot product
	__forceinline float8v dot(const float8v& lhs, const float8v& rhs) {
		return _mm256_dp_ps(lhs.value(), rhs.value(), internal::dotproduct_mask<N>::components | 0x1);
	}



	//! Return 4-component dot product
	__forceinline float8v dot4(const float8v& lhs, const float8v& rhs) {
		return _mm256_dp_ps(lhs.value(), rhs.value(), 0xf1);
	}

	//! Return 3-component dot product, using first 3 elements, ignoring the 4th
	__forceinline float8v dot3(const float8v& lhs, const float8v& rhs) {
		return _mm256_dp_ps(lhs.value(), rhs.value(), 0x71);
	}

	//! Return 2-component dot product, using first 2 elements, ignoring the 2 last
	__forceinline float8v dot2(const float8v& lhs, const float8v& rhs) {
		return _mm256_dp_ps(lhs.value(), rhs.value(), 0x31);
	}


  //! Merges x and y, selecting 2 floats from x and 2 from y.
  // x[x0], x[x1], y[y0], y[y1]
  template <size_t x0, size_t x1, size_t y0, size_t y1>
  __forceinline float8v
  shuffle(const float8v& x, const float8v& y) {
    return _mm256_shuffle_ps(x.value(), y.value(), _MM_SHUFFLE(y1, y0, x1, x0));
  }




  //! Compute the sum of r0..3. Returned result is [sum(r0), sum(r1), sum(r2), sum(r3)]
  __forceinline float8v
  sum_row(const float8v& r0, const float8v& r1, const float8v& r2, const float8v& r3) {
		float8v r01 = _mm256_hadd_ps(r0.value(), r1.value());
		float8v r23 = _mm256_hadd_ps(r2.value(), r3.value());
		float8v sum = _mm256_hadd_ps(r01.value(), r23.value());
		return sum;
#if 0
    float8v ab = shuffle<0,1,0,1>(r0, r1) + shuffle<2,3,2,3>(r0, r1);  // {r0_02, r0_13, r1_02, r1_13}
    float8v cd = shuffle<0,1,0,1>(r2, r3) + shuffle<2,3,2,3>(r2, r3);  // {r2_02, r2_13, r3_02, r3_13}

    float8v sum = shuffle<0,2,0,2>(ab, cd) + shuffle<1,3,1,3>(ab, cd);
    return sum;
#endif
  }




  //! Compute the sum of the 4 elements in r0 and r1.  Returned result is [sum(r0), sum(r1), x, x], x i undefined
  __forceinline float8v
  sum_row(const float8v &r0, const float8v& r1) {
		float8v r01 = _mm256_hadd_ps(r0.value(), r1.value());  // r0_01, r0_23, r1_01, r1_23
		float8v sum = _mm256_hadd_ps(r01.value(), r01.value()); // r0_0123, r1_0123, r0, r1
		return sum;
#if 0
    float8v tmp = shuffle<2, 3, 0, 1>(r0, r1);   // {r0_2, r0_3, r1_0, r1_1}
    float8v tmp0= r0 + tmp;  //  {r0_02, r0_13, x    , x    }
    float8v tmp1= r1 + tmp;  //  {x    , x    , r1_02, r1_13}
    tmp = shuffle<1, 1, 3, 3>(tmp0, tmp1);  // {r0_13, x, r1_13, x}
    tmp0 += tmp;  // {sum(r0), x, x, x}
    tmp1 += tmp;  // {x, x, sum(r1), x}

    tmp1 = _mm_movehl_ps(tmp1.value(), tmp1.value());  // {sum(r1), x, sum(r1), x}
    return _mm_unpacklo_ps(tmp0.value(), tmp1.value());  // {sum(r0), sum(r1), x, x}
#endif
  }



  //! Compute the sum of the 4 elements in r.  Returned result is [sum(r), x, x, x], x i undefined
  __forceinline float8v
  sum_row(const float8v &r) {
		float8v t = _mm256_hadd_ps(r.value(), r.value());
		t = _mm_hadd_ps(t.value(), t.value());
		return t;  // sum in all 4 elements
#if 0
    float8v tmp = r.value() + _mm_movehl_ps(r.value(), r.value());    // { r_02, r_13, x, x} ignore upper 2 floats
    float8v tmp2 = shuffle<1, 0, 0, 1>(tmp, tmp); //_mm_set_ss(tmp.m128_f32[1]);
    return _mm_add_ss(tmp.value(), tmp2.value());
#endif
  }


}}}