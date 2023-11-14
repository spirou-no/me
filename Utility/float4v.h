#pragma once

#include <utility/int4v.h>
#include <utility/mask4v.h>
#include <utility/sse_internal.h>
#include <xmmintrin.h>   // SSE
#include <emmintrin.h>   // SSE2
#include <pmmintrin.h>   // SSE3
#include <smmintrin.h>   // SSE4



namespace oyrke { namespace algorithm { namespace sse {
  
#if 0
    class float4v_traits {

        static float4v infinites();
        static float4v quiet_NaNs();
        static float4v signalling_NaNs();
        static float4v zeros();
        static float4v ones();
    };
#endif

	class float4v {
	private:
		__m128 value_;

	public:
		typedef float value_type;
		typedef __m128 impl_type;
		enum { size = 4 } size_type;

		__forceinline float4v()                                             { value_ = zeros().value(); }
		__forceinline float4v(__m128 x)                                     { value_ = x; }
		__forceinline float4v(float x)                                      { value_ = _mm_set_ps1(x); }
		//__forceinline float4v(int x)                                      { value_ = _mm_set_ps1(float(x)); }
		__forceinline float4v(float x, float y, float z, float w)           { value_ = _mm_setr_ps(x, y, z, w); }
		__forceinline float4v(const float data[4], UNALIGNED)               { value_ = _mm_loadu_ps(data); }
		__forceinline float4v(const float data[4], ALIGNED)                 { value_ = _mm_load_ps(data); }
		__forceinline explicit float4v(const int4v& data)					{ value_ = _mm_cvtepi32_ps(data.value()); }
		__forceinline float4v& operator=(const float4v& rhs)                { value_ = rhs.value_; return *this; }
		__forceinline float4v& operator=(const int4v& rhs)                  { value_ = float4v(rhs).value_; return *this; }
		__forceinline float4v& operator=(float rhs)							{ value_ = _mm_set_ps1(rhs); return *this; }

		static __forceinline float4v create0(float x)                       { return _mm_set_ss(x); }
		static __forceinline float4v create_all(float x)                    { return float4v(x); }
		static __forceinline float4v create_aligned(const float data[4])    { return float4v(data, aligned); }
		static __forceinline float4v create_unaligned(const float data[4])  { return float4v(data, unaligned); }
		template <size_t N> static __forceinline float4v create_replicated(const float4v& x) {
			if constexpr (N > 0) {
				return _mm_shuffle_ps(x.value(), x.value(), _MM_SHUFFLE(N, N, N, N));
			}
			else {
				return x.value_;
			}
		}
		//template <>  /*static*/ __forceinline          float4v create_replicated<0>(const float4v& x) { return x.value_; }
    static __forceinline float4v create_indexed(const float data[], const int4v& ix) {
			return float4v(data[ix.at<0>()], data[ix.at<1>()], data[ix.at<2>()], data[ix.at<3>()]);
		}


		static __forceinline float4v create_reinterpreted(const int4v&  bitmask) { return _mm_castsi128_ps(bitmask.value()); }
		static __forceinline float4v create_reinterpreted(const mask4v& bitmask) { return bitmask.value(); }
		static __forceinline float4v zeros()                                { return _mm_setzero_ps(); }
		static __forceinline float4v ones()									{ return float4v(1.0f); }
		static __forceinline float4v infinites()							{ __m128i m = _mm_set1_epi32(0x7fe00000); return _mm_castsi128_ps(m); }
		static __forceinline float4v NaNs()									{ __m128i m = _mm_set1_epi32(0x7fefffff); return _mm_castsi128_ps(m); }
		static __forceinline mask4v  signmask()                             { return mask4v(0x80000000); }
		static __forceinline mask4v  valuemask()		                    { return mask4v(0x7fffffff); }
		static __forceinline mask4v  specialmask()		                    { return mask4v(0x7fe00000); }

		//template <size_t N> static __forceinline float4v copy_to_0(const float4v& x)    { rturn _mm_castsi128_ps(_mm_srli_si128(_mm_castps_si128(x.value_), N*4)); } //_mm_shuffle_ps(value_, value_, N); }
    //template <size_t N> static __forceinline float4v copy_to_0(const float4v& x)    { return _mm_shuffle_ps(x.value(), x.value(), _MM_SHUFFLE(N, N, N, N)); } 
    template <size_t N> static __forceinline float4v copy_to_0(const float4v& x)    { return create_replicated<N>(x); } 
		//template <>         /*static*/ __forceinline float4v copy_to_0<0>(const float4v& x) { return x.value_; }
		__forceinline void store(float *p, ALIGNED) const                   { _mm_store_ps(p, value_); }
		__forceinline void store(float *p, UNALIGNED) const                 { _mm_storeu_ps(p, value_); }
		template <size_t N> __forceinline void store(float& dest) const     { _mm_store_ss(&dest, copy_to_0<N>(*this).value_); } //_mm_store_ss(p, _mm_shuffle_ps(value_, value_, _MM_SHUFFLE(N,N,N,N))); }

		__forceinline __m128 value() const                                  { return value_; }
		__forceinline float  operator[](size_t i) const                     { return value_.m128_f32[i]; }
		__forceinline float& operator[](size_t i)                           { return value_.m128_f32[i]; }
		template <size_t N> __forceinline float at() const                  { return _mm_cvtss_f32(copy_to_0<N>(*this).value_); } 
    
		template <size_t N> __forceinline void set_at(float x)              { value_ = _mm_insert_ps(value_, float4v(x).value_, _MM_MK_INSERTPS_NDX(0,N,0)); } 
		template <size_t N> __forceinline void clear()                      { value_ = _mm_insert_ps(value_, zeros().value_, N); }

		__forceinline int4v reinterpret_bits() const                        { return int4v(value_); }

		__forceinline const float4v& operator+=(const float4v& rhs)         { value_ = _mm_add_ps(value_, rhs.value_); return *this; }
		__forceinline const float4v& operator*=(const float4v& rhs)         { value_ = _mm_mul_ps(value_, rhs.value_); return *this; }
		__forceinline const float4v& operator-=(const float4v& rhs)         { value_ = _mm_sub_ps(value_, rhs.value_); return *this; }
		__forceinline const float4v& operator/=(const float4v& rhs)         { value_ = _mm_div_ps(value_, rhs.value_); return *this; }

		__forceinline bool equals(const float4v& rhs) const					{ return mask4v(_mm_cmpeq_ps(value_, rhs.value_)).all_true(); }
	};

    __forceinline void streaming_store(const float4v& x, float *p) { 
        // assert((reinterpret_cast<std::ptrdiff_t>(p) & 0xf)  == 0);
        _mm_stream_ps(p, x.value()); 
    }
    
    __forceinline float4v streaming_load(const float* p) { 
        // assert((reinterpret_cast<std::ptrdiff_t>(p) & 0xf)  == 0);
        int4v tmp(_mm_stream_load_si128(reinterpret_cast<__m128i *>(const_cast<float *>(p))));
        return float4v::create_reinterpreted(tmp);
    }

    __forceinline void prefetch(const void *p) { _mm_prefetch(reinterpret_cast<const char*>(p), 0); }
    // move these to generic sse
    __forceinline void store_fence()  { _mm_sfence(); } 
    __forceinline void load_fence()   { _mm_lfence(); }
    __forceinline void memory_fence() { _mm_mfence(); }

	__forceinline float4v operator+(const float4v& lhs, const float4v& rhs) { return _mm_add_ps(lhs.value(), rhs.value()); }
	__forceinline float4v operator+(const float4v& lhs, float rhs)			{ return lhs + float4v(rhs); }
	__forceinline float4v operator+(float lhs, const float4v& rhs)          { return rhs + lhs; }

	__forceinline float4v operator-(const float4v& lhs, const float4v& rhs) { return _mm_sub_ps(lhs.value(), rhs.value()); }
	__forceinline float4v operator-(const float4v& lhs, float rhs)			{ return lhs - float4v(rhs); }
	__forceinline float4v operator-(float lhs, const float4v& rhs)          { return float4v(lhs) - rhs; }
	__forceinline float4v operator-(const float4v& x)						{ return _mm_xor_ps(x.value(), float4v::signmask().value()); }

	__forceinline float4v operator*(const float4v& lhs, const float4v& rhs) { return _mm_mul_ps(lhs.value(), rhs.value()); }
	__forceinline float4v operator*(const float4v& lhs, float rhs)			{ return lhs * float4v(rhs); }
	__forceinline float4v operator*(float lhs, const float4v& rhs)          { return rhs * lhs; }

	__forceinline float4v operator/(const float4v& lhs, const float4v& rhs) { return _mm_div_ps(lhs.value(), rhs.value()); }
	__forceinline float4v operator/(const float4v& lhs, float rhs)			{ return lhs / float4v(rhs); }
	__forceinline float4v operator/(float lhs, const float4v& rhs)          { return float4v(lhs) / rhs; }

	__forceinline mask4v  operator==(const float4v& lhs, const float4v& rhs) { return _mm_cmpeq_ps(lhs.value(), rhs.value()); }
	__forceinline mask4v  operator==(const float4v& lhs, float rhs)			 { return lhs == float4v(rhs); }
	__forceinline mask4v  operator==(float lhs, const float4v& rhs)          { return float4v(lhs) == rhs; }

	__forceinline mask4v  operator!=(const float4v& lhs, const float4v& rhs) { return _mm_cmpneq_ps(lhs.value(), rhs.value()); }
	__forceinline mask4v  operator!=(const float4v& lhs, float rhs)			 { return lhs != float4v(rhs); }
	__forceinline mask4v  operator!=(float lhs, const float4v& rhs)          { return float4v(lhs) != rhs; }
		
	__forceinline mask4v  operator<(const float4v& lhs, const float4v& rhs)	 { return _mm_cmplt_ps(lhs.value(), rhs.value()); }
	__forceinline mask4v  operator<(const float4v& lhs, float rhs)			 { return lhs < float4v(rhs); }
	__forceinline mask4v  operator<(float lhs, const float4v& rhs)           { return float4v(lhs) < rhs; }

	__forceinline mask4v  operator<=(const float4v& lhs, const float4v& rhs) { return _mm_cmple_ps(lhs.value(), rhs.value()); }
	__forceinline mask4v  operator<=(const float4v& lhs, float rhs)			 { return lhs <= float4v(rhs); }
	__forceinline mask4v  operator<=(float lhs, const float4v& rhs)          { return float4v(lhs) <= rhs; }

	__forceinline mask4v  operator>(const float4v& lhs, const float4v& rhs)  { return _mm_cmpgt_ps(lhs.value(), rhs.value()); }
	__forceinline mask4v  operator>(const float4v& lhs, float rhs)			 { return lhs > float4v(rhs); }
	__forceinline mask4v  operator>(float lhs, const float4v& rhs)           { return float4v(lhs) > rhs; }

	__forceinline mask4v  operator>=(const float4v& lhs, const float4v& rhs) { return _mm_cmpge_ps(lhs.value(), rhs.value()); }
	__forceinline mask4v  operator>=(const float4v& lhs, float rhs)			 { return lhs >= float4v(rhs); }
	__forceinline mask4v  operator>=(float lhs, const float4v& rhs)          { return float4v(lhs) >= rhs; }

	
	__forceinline float4v blend(const float4v& x0, const float4v& x1) { 
		return _mm_or_ps(x0.value(), x1.value());
	}

	__forceinline float4v blend(const float4v& x0, const float4v& x1, const float4v& x2) { 
		return blend(blend(x0, x1), x2);
	}

	__forceinline float4v blend(const float4v& x0, const float4v& x1, const float4v& x2, const float4v& x3) { 
		return blend(blend(x0, x1), blend(x2, x3));
	}

	/* Shift elements in lhs up N times, e.g. out[i]=x[i-N]; out[i]=0 if i<N */
	template <size_t N>
	__forceinline float4v shift_up(const float4v& x) {
		__m128i xi = _mm_castps_si128(x.value());
		xi = _mm_slli_si128(xi, N*sizeof(float)*8);
		return _mm_castsi128_ps(xi);
	}

	/* Shift elements in lhs up N times, e.g. out[i]=x[i+N]; out[i]=0 if i>3-N */
	template <size_t N>
	__forceinline float4v shift_down(const float4v& x) {
		__m128i xi = _mm_castps_si128(x.value());
		xi = _mm_srli_si128(xi, N*sizeof(float)*8);
		return _mm_castsi128_ps(xi);
	}

		
	__forceinline float4v select(const mask4v& mask, const float4v& a, const float4v& b) {
		return _mm_or_ps(_mm_and_ps(mask.value(), a.value()), _mm_andnot_ps(mask.value(), b.value()));
		//return _mm_blendv_ps(b.value(), a.value(), mask.value());
		//return intern::blend(mask.value(), a.value(), b.value());
	}

	__forceinline float4v select(const mask4v& mask, const float4v& a) {
		//			return _mm_or_ps(_mm_and_ps(mask.value(), a.value()), _mm_andnot_ps(mask.value(), b.value()));
		//return _mm_blendv_ps(b.value(), a.value(), mask.value());
		return _mm_and_ps(mask.value(), a.value());
	}

	__forceinline float4v select_not(const mask4v& mask, const float4v& a) {
		//			return _mm_or_ps(_mm_and_ps(mask.value(), a.value()), _mm_andnot_ps(mask.value(), b.value()));
		//return _mm_blendv_ps(b.value(), a.value(), mask.value());
		return _mm_andnot_ps(mask.value(), a.value());
	}

	__forceinline float4v if_then_else(const mask4v& mask, const float4v& a, const float4v& b) {
		return select(mask, a, b);
	}

	__forceinline float4v if_then(const mask4v& mask, const float4v& a) {
		return select(mask, a);
	}

  //! Merges x and y, selecting 2 floats from x and 2 from y.
  // x[x0], x[x1], y[y0], y[y1]
  template <size_t x0, size_t x1, size_t y0, size_t y1>
  __forceinline float4v
  shuffle(const float4v& x, const float4v& y) {
    return _mm_shuffle_ps(x.value(), y.value(), _MM_SHUFFLE(y1, y0, x1, x0));
  }




  //! Compute the sum of r0..3. Returned result is [sum(r0), sum(r1), sum(r2), sum(r3)]
  __forceinline float4v
  sum_row(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3) {
		float4v r01 = _mm_hadd_ps(r0.value(), r1.value());
		float4v r23 = _mm_hadd_ps(r2.value(), r3.value());
		float4v sum = _mm_hadd_ps(r01.value(), r23.value());
		return sum;
#if 0
    float4v ab = shuffle<0,1,0,1>(r0, r1) + shuffle<2,3,2,3>(r0, r1);  // {r0_02, r0_13, r1_02, r1_13}
    float4v cd = shuffle<0,1,0,1>(r2, r3) + shuffle<2,3,2,3>(r2, r3);  // {r2_02, r2_13, r3_02, r3_13}

    float4v sum = shuffle<0,2,0,2>(ab, cd) + shuffle<1,3,1,3>(ab, cd);
    return sum;
#endif
  }




  //! Compute the sum of the 4 elements in r0 and r1.  Returned result is [sum(r0), sum(r1), x, x], x i undefined
  __forceinline float4v
  sum_row(const float4v &r0, const float4v& r1) {
		float4v r01 = _mm_hadd_ps(r0.value(), r1.value());  // r0_01, r0_23, r1_01, r1_23
		float4v sum = _mm_hadd_ps(r01.value(), r01.value()); // r0_0123, r1_0123, r0, r1
		return sum;
#if 0
    float4v tmp = shuffle<2, 3, 0, 1>(r0, r1);   // {r0_2, r0_3, r1_0, r1_1}
    float4v tmp0= r0 + tmp;  //  {r0_02, r0_13, x    , x    }
    float4v tmp1= r1 + tmp;  //  {x    , x    , r1_02, r1_13}
    tmp = shuffle<1, 1, 3, 3>(tmp0, tmp1);  // {r0_13, x, r1_13, x}
    tmp0 += tmp;  // {sum(r0), x, x, x}
    tmp1 += tmp;  // {x, x, sum(r1), x}

    tmp1 = _mm_movehl_ps(tmp1.value(), tmp1.value());  // {sum(r1), x, sum(r1), x}
    return _mm_unpacklo_ps(tmp0.value(), tmp1.value());  // {sum(r0), sum(r1), x, x}
#endif
  }



  //! Compute the sum of the 4 elements in r.  Returned result is [sum(r), x, x, x], x i undefined
  __forceinline float4v
  sum_row(const float4v &r) {
		float4v t = _mm_hadd_ps(r.value(), r.value());
		t = _mm_hadd_ps(t.value(), t.value());
		return t;  // sum in all 4 elements
#if 0
    float4v tmp = r.value() + _mm_movehl_ps(r.value(), r.value());    // { r_02, r_13, x, x} ignore upper 2 floats
    float4v tmp2 = shuffle<1, 0, 0, 1>(tmp, tmp); //_mm_set_ss(tmp.m128_f32[1]);
    return _mm_add_ss(tmp.value(), tmp2.value());
#endif
  }


  		//! Return 4-component dot product
	__forceinline float4v dot4(const float4v& x, const float4v& y0, const float4v& y1, const float4v& y2, const float4v& y3) {
		__m128 p = _mm_dp_ps(x.value(), y0.value(), 0xf1);
		__m128 q = _mm_dp_ps(x.value(), y1.value(), 0xf2);
		__m128 r = _mm_dp_ps(x.value(), y2.value(), 0xf4);
		__m128 s = _mm_dp_ps(x.value(), y3.value(), 0xf8);

		__m128 result = _mm_or_ps(_mm_or_ps(p, q), _mm_or_ps(r, s));
		return result;
	}

	__forceinline float4v dot4(const float4v& x, const float4v& y0, const float4v& y1) {
		__m128 p = _mm_dp_ps(x.value(), y0.value(), 0xf1);
		__m128 q = _mm_dp_ps(x.value(), y1.value(), 0xf2);

		__m128 result = _mm_or_ps(p, q);
		return result;
	}

	namespace internal {
		template <size_t N>
		struct dotproduct_mask {
			enum mask { components = (0x0f >> (4 - N)) << (4 + (4 - N)) };
		};
	}

	//! Return 4-component dot product
	template <size_t N>
	__forceinline float4v dot(const float4v& x, const float4v& y0, const float4v& y1, const float4v& y2, const float4v& y3) {
		__m128 p = _mm_dp_ps(x.value(), y0.value(), internal::dotproduct_mask<N>::components | 0x01);
		__m128 q = _mm_dp_ps(x.value(), y1.value(), internal::dotproduct_mask<N>::components | 0x02);
		__m128 r = _mm_dp_ps(x.value(), y2.value(), internal::dotproduct_mask<N>::components | 0x04);
		__m128 s = _mm_dp_ps(x.value(), y3.value(), internal::dotproduct_mask<N>::components | 0x08);

		__m128 result = _mm_or_ps(_mm_or_ps(p, q), _mm_or_ps(r, s));
		return result;
	}

	template <size_t N>
	__forceinline float4v dot(const float4v& x, const float4v& y0, const float4v& y1) {
		__m128 p = _mm_dp_ps(x.value(), y0.value(), internal::dotproduct_mask<N>::components | 0x01);
		__m128 q = _mm_dp_ps(x.value(), y1.value(), internal::dotproduct_mask<N>::components | 0x02);

		__m128 result = _mm_or_ps(p, q);
		return result;
	}

	template <size_t N>
	//! Return N-component dot product
	__forceinline float4v dot(const float4v& lhs, const float4v& rhs) {
		return _mm_dp_ps(lhs.value(), rhs.value(), internal::dotproduct_mask<N>::components | 0x1);
	}



	//! Return 4-component dot product
	__forceinline float4v dot4(const float4v& lhs, const float4v& rhs) {
		return _mm_dp_ps(lhs.value(), rhs.value(), 0xf1);
	}

	//! Return 3-component dot product, using first 3 elements, ignoring the 4th
	__forceinline float4v dot3(const float4v& lhs, const float4v& rhs) {
		return _mm_dp_ps(lhs.value(), rhs.value(), 0x71);
	}

	//! Return 2-component dot product, using first 2 elements, ignoring the 2 last
	__forceinline float4v dot2(const float4v& lhs, const float4v& rhs) {
		return _mm_dp_ps(lhs.value(), rhs.value(), 0x31);
	}

}}}