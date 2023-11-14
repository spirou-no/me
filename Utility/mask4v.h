#pragma once

#include <utility/sse_internal.h>
#include <xmmintrin.h>   // SSE
#include <emmintrin.h>   // SSE2
#include <pmmintrin.h>   // SSE3
#include <smmintrin.h>   // SSE4


namespace oyrke { namespace algorithm { namespace sse {
	// Mask array for float4v and int32_4v
	class mask4v {
		__m128 value_;

		static __forceinline __m128  cast(__m128i m) { return _mm_castsi128_ps(m); }
		static __forceinline __m128i cast(__m128 m)  { return _mm_castps_si128(m); }
	public:

		__forceinline mask4v()						 : value_(cast(_mm_setzero_si128())) {}
		__forceinline mask4v(__m128 mask)  : value_(mask) {}
		__forceinline mask4v(__m128i mask) : value_(cast(mask)) {}
		__forceinline explicit mask4v(unsigned int mask) : value_(cast(_mm_set1_epi32(int(mask)))) {}

		__forceinline __m128 value() const { return value_; }
		__forceinline __m128i value_i() const { return cast(value_); }

		static __forceinline mask4v zeros() { return mask4v(); }
		static __forceinline mask4v falses() { return mask4v(); }
		static __forceinline mask4v trues() { return xFFFs(); }
		static __forceinline mask4v ones() { return xFFFs(); }
		static __forceinline mask4v xFFFs()  { __m128i x = _mm_setzero_si128(); return _mm_cmpeq_epi32(x, x); }  // equality also for uninited x

		template <size_t N> __forceinline void set();
		template <size_t N> __forceinline void clear();
		template <size_t N> __forceinline __m128 at() const;
		template <size_t N> static __forceinline mask4v mask() { 
			//__m128i m = { 0xffffffff, 0x0, 0x0, 0x0 }; 
			// m = _mm_srli_si128(ones.value_i(), 3*4); // alternative way to set r0:=0xffffffff, r1:=r2:=r3:=0
			// or
			unsigned int x = ~0u;
			__m128i m = cast(_mm_set_ss(reinterpret_cast<const float&>(x)));
			__m128i res = _mm_slli_si128(m, 4*N); 
			return res;
		}

		__forceinline mask4v& operator=(const mask4v& rhs)	{ value_ = rhs.value_; return *this; }
		__forceinline mask4v& operator&=(const mask4v& rhs) { value_ = _mm_and_ps(value_, rhs.value_) ; return *this; }
		__forceinline mask4v& operator|=(const mask4v& rhs) { value_ = _mm_or_ps(value_, rhs.value_) ; return *this; }
		__forceinline mask4v& operator^=(const mask4v& rhs) { value_ = _mm_xor_ps(value_, rhs.value_) ; return *this; }

		__forceinline int movemask() const									{ return _mm_movemask_ps(value_); }
		__forceinline bool all_false() const								{ return movemask() == 0; }
		__forceinline bool all_true() const									{ return movemask() == 0xf; }
		__forceinline bool any_false() const								{ return movemask() != 0xf; }
		__forceinline bool any_true() const									{ return movemask() != 0; }
	};


	__forceinline mask4v select(const mask4v& mask, const mask4v& a, const mask4v& b) {
		// return _mm_or_si128(_mm_and_si128(mask.value_i(), a.value_i()), _mm_andnot_si128(mask.value_i(), b.value_i()));
		return intern::blend(mask.value(), a.value(), b.value());
	}


	__forceinline mask4v operator& (const mask4v& lhs, const mask4v& rhs) { return _mm_and_ps(lhs.value(), rhs.value()); }
	__forceinline mask4v operator&&(const mask4v& lhs, const mask4v& rhs) { return lhs & rhs; }
	__forceinline mask4v operator| (const mask4v& lhs, const mask4v& rhs) { return _mm_or_ps(lhs.value(), rhs.value()); }
	__forceinline mask4v operator||(const mask4v& lhs, const mask4v& rhs) { return lhs | rhs; }
	__forceinline mask4v operator^ (const mask4v& lhs, const mask4v& rhs) { return _mm_xor_ps(lhs.value(), rhs.value()); }
//	__forceinline	mask4v operator~ (const mask4v& x) { return _mm_xor_ps(x.value(), mask4v::ones().value()); }
	__forceinline	mask4v operator~ (const mask4v& x)											{ return _mm_xor_ps(x.value(), mask4v::ones().value()); }

	__forceinline	mask4v operator! (const mask4v& x) { return ~x; }
	__forceinline mask4v operator==(const mask4v& lhs, const mask4v& rhs) { return _mm_cmpeq_epi32(lhs.value_i(), rhs.value_i()); }
	__forceinline mask4v operator!=(const mask4v& lhs, const mask4v& rhs) { return ~(lhs == rhs); }  
	
}}}
