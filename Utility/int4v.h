#pragma once

#include <utility/mask4v.h>
#include <utility/sse_internal.h>
#include <xmmintrin.h>   // SSE
#include <emmintrin.h>   // SSE2
#include <pmmintrin.h>   // SSE3
#include <smmintrin.h>   // SSE4



namespace oyrke { namespace algorithm { namespace sse {

  
	class int4v {
		__m128i value_;

	public:
		__forceinline int4v()												{ value_ = _mm_setzero_si128(); }
		__forceinline int4v(__m128i x)										{ value_ = x; }
		__forceinline int4v(int i)											{ value_ = _mm_set1_epi32(i); }
		__forceinline int4v(int i0, int i1, int i2, int i3)                 { value_ = _mm_setr_epi32(i0, i1, i2, i3); }
		__forceinline explicit int4v(const mask4v& mask)                    { value_ = mask.value_i(); }

		static __forceinline int4v zeros()									{ return int4v(); }
		static __forceinline int4v ones()									{ return int4v(1); }
		static __forceinline int4v xFFFs()								    { __m128i m; return _mm_cmpeq_epi32(m, m); /* return int4v(-1); */}
		static __forceinline int4v create0(int i0)                          { return _mm_castps_si128(_mm_set_ss(reinterpret_cast<float&>(i0))); }


		template <size_t N> __forceinline __m128i copy_to_0() const         { return _mm_srli_si128(value_, N*4); } //_mm_shuffle_ps(value_, value_, N); }
		template <>         __forceinline __m128i copy_to_0<0>() const      { return value_; }

		__forceinline __m128i value() const                                 { return value_; }
		__forceinline int operator[](int i)                                 { return value_.m128i_i32[i]; }
		template <size_t N> __forceinline int at()    const                 { return _mm_cvtsi128_si32(copy_to_0<N>()); }
		template <>         __forceinline int at<0>() const                 { return _mm_cvtsi128_si32(value_); }

		// SSE4
		//__forceinline int4v operator*(const int4v& rhs) const          { return _mm_mul_epi32(value_, rhs.value_); }
		//__forceinline int4v operator*(int rhs) const                   { return *this * int4v(rhs); }

		__forceinline const int4v& operator+=(const int4v& rhs)             { value_ = _mm_add_epi32(value_, rhs.value_); return *this; }
		__forceinline const int4v& operator+=(int rhs)                      { return *this += int4v(rhs); }
		__forceinline const int4v& operator-=(const int4v& rhs)             { value_ = _mm_sub_epi32(value_, rhs.value_); return *this; }
		__forceinline const int4v& operator-=(int rhs)                      { return *this -= int4v(rhs); }
		__forceinline const int4v& operator*=(const int4v& rhs)             { value_ = _mm_mullo_epi32(value_, rhs.value_); return *this; }
		__forceinline const int4v& operator*=(int rhs)                      { return *this *= int4v(rhs); }
		__forceinline const int4v& operator&=(const int4v& rhs)             { value_ = _mm_and_si128(value_, rhs.value_); return *this; }
		__forceinline const int4v& operator&=(const mask4v& rhs)            { value_ = _mm_and_si128(value_, rhs.value_i()); return *this; }
		__forceinline const int4v& operator&=(int rhs)                      { return *this &= int4v(rhs); }
		__forceinline const int4v& operator|=(const int4v& rhs)             { value_ = _mm_or_si128(value_, rhs.value_); return *this; }
		__forceinline const int4v& operator|=(const mask4v& rhs)            { value_ = _mm_or_si128(value_, rhs.value_i()); return *this; }
		__forceinline const int4v& operator|=(int rhs)                      { return *this |= int4v(rhs); }
		__forceinline const int4v& operator^=(const int4v& rhs)             { value_ = _mm_xor_si128(value_, rhs.value_); return *this; }
		__forceinline const int4v& operator^=(const mask4v& rhs)            { value_ = _mm_xor_si128(value_, rhs.value_i()); return *this; }
		__forceinline const int4v& operator^=(int rhs)                      { return *this ^= int4v(rhs); }
		__forceinline const int4v& operator<<=(const int4v& rhs)            { value_ = _mm_sll_epi32(value_, rhs.value_); return *this; }
		__forceinline const int4v& operator<<=(const int rhs)               { value_ = _mm_slli_epi32(value_, rhs); return *this; }
		__forceinline const int4v& operator>>=(const int4v& rhs)            { value_ = _mm_srl_epi32(value_, rhs.value_); return *this; }
		__forceinline const int4v& operator>>=(const int rhs)               { value_ = _mm_srli_epi32(value_, rhs); return *this; }
	};

	__forceinline int4v operator&(const int4v& lhs, const int4v& rhs)	{ return _mm_and_si128(lhs.value(), rhs.value()); }
	__forceinline int4v operator&(const mask4v& lhs, const int4v& rhs)	{ return int4v(lhs) & rhs; }
	__forceinline int4v operator&(const int4v& lhs, const mask4v& rhs)	{ return lhs & int4v(rhs); }
	__forceinline int4v operator|(const int4v& lhs, const int4v& rhs)	{ return _mm_or_si128(lhs.value(), rhs.value()); }
	__forceinline int4v operator|(const mask4v& lhs, const int4v& rhs)	{ return int4v(lhs) | rhs; }
	__forceinline int4v operator|(const int4v& lhs, const mask4v& rhs)	{ return lhs | int4v(rhs); }
	__forceinline int4v operator^(const int4v& lhs, const int4v& rhs)	{ return _mm_xor_si128(lhs.value(), rhs.value()); }
	__forceinline int4v operator^(const mask4v& lhs, const int4v& rhs)	{ return int4v(lhs) ^ rhs; }
	__forceinline int4v operator^(const int4v& lhs, const mask4v& rhs)	{ return lhs ^ int4v(rhs); }
	__forceinline int4v operator~(const int4v& x)						{ return x ^ int4v::xFFFs(); }

	__forceinline int4v operator+(const int4v& lhs, const int4v& rhs)	{ return _mm_add_epi32(lhs.value(), rhs.value()); }
	__forceinline int4v operator-(const int4v& lhs, const int4v& rhs)   { return _mm_sub_epi32(lhs.value(), rhs.value()); }
	__forceinline int4v operator-(const int4v& x)						{ return int4v::zeros() - x; }
	__forceinline int4v operator*(const int4v& lhs, const int4v& rhs)   { return _mm_mullo_epi32(lhs.value(), rhs.value()); }
	__forceinline mask4v operator==(const int4v& lhs, const int4v& rhs) { return _mm_cmpeq_epi32(lhs.value(), rhs.value()); }
	__forceinline mask4v operator!=(const int4v& lhs, const int4v& rhs) { return ~(lhs == rhs); }
	__forceinline mask4v operator<(const int4v& lhs, const int4v& rhs)	{ return _mm_cmplt_epi32(lhs.value(), rhs.value()); }
	__forceinline mask4v operator>(const int4v& lhs, const int4v& rhs)	{ return _mm_cmpgt_epi32(lhs.value(), rhs.value()); }
	__forceinline mask4v operator<=(const int4v& lhs, const int4v& rhs) { return ~(lhs > rhs); } // no _mm_cmple_epi32(lhs.value(), rhs.value()); }
	__forceinline mask4v operator>=(const int4v& lhs, const int4v& rhs) { return ~(lhs < rhs); } // no _mm_cmpge_epi32(lhs.value(), rhs.value()); }

	__forceinline int4v operator<<(const int4v& lhs, const int4v& rhs)  { int4v r = lhs; r <<= rhs; return r; }
	__forceinline int4v operator<<(const int4v& lhs, const int rhs)		{ int4v r = lhs; r <<= rhs; return r; }
	__forceinline int4v operator>>(const int4v& lhs, const int4v& rhs)  { int4v r = lhs; r >>= rhs; return r; }
	__forceinline int4v operator>>(const int4v& lhs, const int rhs)     { int4v r = lhs; r >>= rhs; return r; }

	//! Merges x and y, selecting 2 floats from x and 2 from y.
	// x[x0], x[x1], y[y0], y[y1]
	template <size_t x0, size_t x1, size_t y0, size_t y1>
	__forceinline int4v
	shuffle(const int4v& x, const int4v& y) {
		return _mm_shuffle_epi32(x.value(), y.value(), _MM_SHUFFLE(y1, y0, x1, x0));
	}


	// if_then_else
	__forceinline int4v select(const mask4v& mask, const int4v& a, const int4v& b) {
		return intern::blend_epi32(mask.value_i(), a.value(), b.value());
	}

	// if_then_else_zero
	__forceinline int4v select(const mask4v& mask, const int4v& a) {
		return _mm_and_si128(mask.value_i(), a.value());
	}

    // if_not_then_else_zero
	__forceinline int4v select_not(const mask4v& mask, const int4v& a) {
		return _mm_andnot_si128(mask.value_i(), a.value());
	}


	__forceinline int4v abs(const int4v& a) { return _mm_abs_epi32(a.value()); }

	__forceinline int4v changesign(const int4v& a, const int4v& sign) {
		mask4v negate = sign < 0;
		// negate in 2-complement is same as subtract 1 then toggle all bits
		int4v t = select(negate, a-1, a) ^ negate;  
		return t;
	}

	__forceinline int4v copysign(const int4v& a, const int4v& sign) {
		int4v pos = abs(a);
		mask4v negate = sign < 0;
		return select(negate, -pos, pos);
	}


	__forceinline int4v min(const int4v& a, const int4v& b) {
		return _mm_min_epi32(a.value(), b.value()); // SSE4 select(a < b, a, b); 
	}

	__forceinline int4v min(const int4v& a, const int4v& b, const int4v& c) {
		return min(a, min(b, c));
	}

	__forceinline int4v min(const int4v& a, const int4v& b, const int4v& c, const int4v& d) {
		return min(min(a, b), min(c, d));
	}

	__forceinline int4v max(const int4v& a, const int4v& b) {
		return _mm_max_epi32(a.value(), b.value()); // SSE4 select(a > b, a, b);
	}

 	__forceinline int4v max(const int4v& a, const int4v& b, const int4v& c) {
		return max(a, max(b, c));
	}

	__forceinline int4v max(const int4v& a, const int4v& b, const int4v& c, const int4v& d) {
		return max(max(a, b), max(c, d));
	}

	__forceinline int4v clip(const int4v& x, const int4v& lo, const int4v& hi) {
		return min(max(x, lo), hi);
	}


}}}