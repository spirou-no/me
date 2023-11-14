#pragma once

#include <utility/sse_internal.h>
#include <xmmintrin.h>   // SSE
#include <emmintrin.h>   // SSE2
#include <pmmintrin.h>   // SSE3
#include <smmintrin.h>   // SSE4


namespace oyrke { namespace algorithm { namespace sse {

	/* 	TODO
	bit algo to create mask for element 0
	struct bitmask
	*/
	template <size_t N>
	struct maskbits {
		enum { value = 0x1 << N | maskbits<N-1>::value } bits;
	};

	template <>
	struct maskbits<0> {
		enum { value = 0x1  } bits;
	};



	struct m128_ps_bit_traits {
		typedef __m128 impl_t;
		__forceinline static impl_t and(const impl_t& lhs, const impl_t& rhs) { return _mm_and_ps(lhs, rhs); }
		__forceinline static impl_t andnot(const impl_t& lhs, const impl_t& rhs) { return _mm_andnot_ps(lhs, rhs); }
		__forceinline static impl_t or(const impl_t& lhs, const impl_t& rhs) { return _mm_or_ps(lhs, rhs); }
		__forceinline static impl_t xor(const impl_t& lhs, const impl_t& rhs) { return _mm_xor_ps(lhs, rhs); }
	};

	struct m128_ps_const_traits {
		typedef __m128 impl_t;
		typedef __m128i int_impl_t;

		__forceinline static impl_t cast(const int_impl_t& x) { return _mm_castsi128_ps(x); }
		__forceinline static int_impl_t cast(const impl_t& x) { return _mm_castps_si128(x); }
		__forceinline static impl_t zeros() { return _mm_setzero_ps(); }
		__forceinline static impl_t ones() { return _mm_set1_ps(1.0f); }
		__forceinline static impl_t xFFFs() { impl_t x; return _mm_cmpeq_ps(x, x); }
		
		template <size_t N> 
		static __forceinline impl_t mask() { 
			unsigned int x = ~0u;
			__m128i m = _mm_castps_si128(_mm_set_ss(reinterpret_cast<const float&>(x)));
			__m128i res = _mm_slli_si128(m, 4*N); 
			return res;
		}
		__forceinline static impl_t signmask() { return cast(_mm_set1_epi32(0x80000000)); }
		__forceinline static impl_t valuemask() { return cast(_mm_set1_epi32(0x80000000));  }

		__forceinline impl_t nan() { return cast(_mm_set1_epi32(0x7f8fffff)); }
		__forceinline impl_t inf() { return cast(_mm_set1_epi32(0x7f800000)); }
	};

	struct m128_ps_compare_traits {
		typedef __m128 impl_t;
		__forceinline static impl_t eq(const impl_t& lhs, const impl_t& rhs) { return _mm_cmpeq_ps(lhs, rhs); }
		__forceinline static impl_t neq(const impl_t& lhs, const impl_t& rhs) { return _mm_cmpneq_ps(lhs, rhs); }
		__forceinline static impl_t lt(const impl_t& lhs, const impl_t& rhs) { return _mm_cmplt_ps(lhs, rhs); }
		__forceinline static impl_t le(const impl_t& lhs, const impl_t& rhs) { return _mm_cmple_ps(lhs, rhs); }
		__forceinline static impl_t gt(const impl_t& lhs, const impl_t& rhs) { return _mm_cmpgt_ps(lhs, rhs); }
		__forceinline static impl_t ge(const impl_t& lhs, const impl_t& rhs) { return _mm_cmpge_ps(lhs, rhs); }
		__forceinline static impl_t min(const impl_t& lhs, const impl_t& rhs) { return _mm_min_ps(lhs, rhs); }
		__forceinline static impl_t max(const impl_t& lhs, const impl_t& rhs) { return _mm_max_ps(lhs, rhs); }
	};

	struct m128_ps_arithmetic_traits {
		typedef __m128 impl_t;
		__forceinline static impl_t plus(const impl_t& lhs, const impl_t& rhs) { return _mm_add_ps(lhs, rhs); }
		__forceinline static impl_t minus(const impl_t& lhs, const impl_t& rhs) { return _mm_sub_ps(lhs, rhs); }
		__forceinline static impl_t multiply(const impl_t& lhs, const impl_t& rhs) { return _mm_mul_ps(lhs, rhs); }
		__forceinline static impl_t divide(const impl_t& lhs, const impl_t& rhs) { return _mm_div_ps(lhs, rhs); }
	};


	struct m128_ps_function_traits {
		typedef __m128 impl_t;
		__forceinline static impl_t sqrt(const impl_t& x) { return _mm_sqrt_ps(x); }
		__forceinline static impl_t rcp(const impl_t& x) { return _mm_rcp_ps(x); }
		__forceinline static impl_t sqrt_rcp(const impl_t& x) { return _mm_rsqrt_ps(x); }
	};

	/*
	min, max
	abs
	signbits
	insert
	extract
	shuffle

	*/


	// Mask array for float4v and int32_4v
	class mask_vec {
	public:
		typedef __m128 impl_type;
		typedef unsigned int value_type;
		const int size = 4;
		const size_t max_movemask = maskbits<size-1>::value;

		typedef m128_ps_bit_traits bit_traits;
		typedef m128_ps_const_traits const_traits;
		typedef m128_ps_compare_traits compare_traits;

	private:
		impl_type value_;

		static __forceinline __m128  cast(__m128i m) { return _mm_castsi128_ps(m); }
		static __forceinline __m128i cast(__m128 m)  { return _mm_castps_si128(m); }
	public:

		__forceinline mask_vec()						 : value_(const_traits::zeros()) {}
		__forceinline mask_vec(impl_type mask)  : value_(mask) {}
		// gone __forceinline explicit mask_vec(unsigned int mask) : value_(cast(_mm_set1_epi32(int(mask)))) {}

		__forceinline __m128 value() const { return value_; }

		static __forceinline mask_vec zeros() { return mask_vec(); }
		static __forceinline mask_vec falses() { return mask_vec(); }
		static __forceinline mask_vec trues() { return xFFFs(); }
		static __forceinline mask_vec ones() { return xFFFs(); }
		static __forceinline mask_vec xFFFs()  { return const_traits::xFFFs(); }  // equality also for uninited x

		template <size_t N> __forceinline void set();
		template <size_t N> __forceinline void clear();
		template <size_t N> __forceinline __m128 at() const;
		template <size_t N> static __forceinline mask_vec mask() { return const_traits<N>::mask(); }

		__forceinline mask_vec& operator=(const mask_vec& rhs)	{ value_ = rhs.value_; return *this; }
		__forceinline mask_vec& operator&=(const mask_vec& rhs) { value_ = bit_traits::and(value_, rhs.value_) ; return *this; }
		__forceinline mask_vec& operator|=(const mask_vec& rhs) { value_ = bit_traits::or(value_, rhs.value_) ; return *this; }
		__forceinline mask_vec& operator^=(const mask_vec& rhs) { value_ = bit_traits::xor(value_, rhs.value_) ; return *this; }

		__forceinline int movemask() const									{ return _mm_movemask_ps(value_); }
		__forceinline bool all_false() const								{ return movemask() == 0; }
		__forceinline bool all_true() const									{ return movemask() == max_movemask; }
		__forceinline bool any_false() const								{ return movemask() != max_movemask; }
		__forceinline bool any_true() const									{ return movemask() != 0; }
	};


	__forceinline mask_vec select(const mask_vec& mask, const mask_vec& a, const mask_vec& b) {
		// return _mm_or_si128(_mm_and_si128(mask.value_i(), a.value_i()), _mm_andnot_si128(mask.value_i(), b.value_i()));
		return intern::blend(mask.value(), a.value(), b.value());
	}


	__forceinline mask_vec operator& (const mask_vec& lhs, const mask_vec& rhs) { return mask_vec::bit_traits::and(lhs.value(), rhs.value()); }
	__forceinline mask_vec operator&&(const mask_vec& lhs, const mask_vec& rhs) { return lhs & rhs; }
	__forceinline mask_vec operator| (const mask_vec& lhs, const mask_vec& rhs) { return mask_vec::bit_traits::or(lhs.value(), rhs.value()); }
	__forceinline mask_vec operator||(const mask_vec& lhs, const mask_vec& rhs) { return lhs | rhs; }
	__forceinline mask_vec operator^ (const mask_vec& lhs, const mask_vec& rhs) { return mask_vec::bit_traits::xor(lhs.value(), rhs.value()); }
//	__forceinline	mask_vec operator~ (const mask_vec& x) { return _mm_xor_ps(x.value(), mask_vec::ones().value()); }
	__forceinline	mask_vec operator~ (const mask_vec& x)											{ return _mm_xor_ps(x.value(), mask_vec::ones().value()); }

	__forceinline	mask_vec operator! (const mask_vec& x) { return ~x; }
	__forceinline mask_vec operator==(const mask_vec& lhs, const mask_vec& rhs) { return mask_vec::compare_traits::eq(lhs.value(), rhs.value()); }
	__forceinline mask_vec operator!=(const mask_vec& lhs, const mask_vec& rhs) { return ~(lhs == rhs); }  
	
}}}
