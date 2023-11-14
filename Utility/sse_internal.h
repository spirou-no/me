#pragma once

#include <xmmintrin.h>   // SSE
#include <emmintrin.h>   // SSE2
#include <pmmintrin.h>   // SSE3
#include <smmintrin.h>   // SSE4

#undef UNALIGNED

#define USE_SSE4

namespace oyrke { namespace algorithm { namespace sse {

	enum ALIGNED   { aligned };
	enum UNALIGNED { unaligned };


	namespace intern {
		// TODO C++20 Give proper name
		__forceinline __m128i or20(const __m128i& a) {
			__m128i t2= _mm_srli_epi64(a, 32);
			__m128i ored = _mm_or_si128(a, t2);  // return here on x64
			t2 = _mm_srli_si128(ored, 2*4);
			return _mm_or_si128(ored, t2);
		}

		__forceinline __m128i and20(const __m128i& a) {
			__m128i t2= _mm_srli_epi64(a, 32);
			__m128i anded = _mm_and_si128(a, t2); // return here on x64
			t2 = _mm_srli_si128(anded, 2*4);
			return _mm_and_si128(anded, t2);
		}

		__forceinline bool     is_all_zero0(const __m128i& a)	{ return _mm_cvtsi128_si32(a) == 0; }
		__forceinline bool     is_zero0(const __m128d& a)	    { return (_mm_movemask_pd(a) & 0x1) == 0; }
		__forceinline bool     is_zero0(const __m128& a)	    { return (_mm_movemask_ps(a) & 0x1) == 0; }
		__forceinline bool     is_all_set0(const __m128i& a)	{ return _mm_cvtsi128_si32(a) == ~int(0); }
		__forceinline bool     is_set0(const __m128d& a)		{ return (_mm_movemask_pd(a) & 0x1) != 0;; }
		__forceinline bool     is_set0(const __m128& a)			{ return (_mm_movemask_ps(a) & 0x1) == 0;; }

		__forceinline bool     is_all_zero(const __m128i& a)	{ return is_all_zero0(or20(a)); }
		__forceinline bool     is_all_zero(const __m128d& a)	{ return _mm_movemask_pd(a) == 0; }
		__forceinline bool     is_all_zero(const __m128& a)		{ return _mm_movemask_ps(a) == 0; }
		__forceinline bool     is_all_set(const __m128i& a)		{ return is_all_set0(and20(a)); }
		__forceinline bool     is_all_set(const __m128d& a)		{ return _mm_movemask_pd(a) == 0x3; }
		__forceinline bool     is_all_set(const __m128& a)		{ return _mm_movemask_ps(a) == 0xf; }
		__forceinline bool     is_any_zero(const __m128i& a)	{ return is_all_zero0(and20(a));	}
		__forceinline bool     is_any_zero(const __m128d& a)	{ return _mm_movemask_pd(a) != 0x3;	}
		__forceinline bool     is_any_zero(const __m128& a)		{ return _mm_movemask_ps(a) != 0xf; }
		__forceinline bool     is_any_set(const __m128i& a)		{ return is_all_set0(or20(a));	}
		__forceinline bool     is_any_set(const __m128d& a)		{ return _mm_movemask_pd(a) != 0;	}
		__forceinline bool     is_any_set(const __m128& a)		{ return _mm_movemask_ps(a) != 0; }

		__forceinline __m128   blend(const __m128& mask, const __m128& a, const __m128& b) {
		#ifdef USE_SSE4
			return _mm_blendv_ps(b, a, mask);
		#else
			return _mm_or_ps(_mm_and_ps(mask, a), _mm_andnot_ps(mask, b));
		#endif
		}

		__forceinline __m128i blend_epi32(const __m128i& mask, const __m128i& a, const __m128i& b) {
		#ifdef USE_SSE4
			return _mm_castps_si128(blend(_mm_castsi128_ps(mask), _mm_castsi128_ps(a), _mm_castsi128_ps(b)));
		#else
			return _mm_or_si128(_mm_and_si128(mask, a), _mm_andnot_si128(mask, b));
		#endif
		}

		__forceinline __m128d compose_01(const __m128d& a, const __m128d& b) { return _mm_move_sd(b, a); }
		__forceinline __m128d compose_00(const __m128d& a, const __m128d& b) { return _mm_unpacklo_pd(b, a); }
		__forceinline __m128d compose_11(const __m128d& a, const __m128d& b) { return _mm_unpackhi_pd(b, a); }

  
		template <size_t ai, size_t bi> __forceinline __m128d
		shuffle(const __m128d& a, const __m128d& b)						{ return _mm_shuffle_pd(a, b, _MM_SHUFFLE2(ai, bi)); }

		struct m128_float {
			static __forceinline __m128 set0(float x)					{ return _mm_set_ss(x); }
			//static __forceinline __m128 { }
		};

		struct m128_int32 {
			//static __forceinline __m128 { }
		};

		struct m128_uint32 {
			//static __forceinline __m128 { }
		};

		struct m128_int16 {
			//static __forceinline __m128 { }
		};

		struct m128_uint16 {
			//static __forceinline __m128 { }
		};

		struct m128_int8 {
			//static __forceinline __m128 { }
		};

		struct m128_uint8 {
			//static __forceinline __m128 { }
		};


		struct m128i_int128 {

			//static __forceinline __m128 { }
		};

		struct m128i_int64 {
			static __forceinline __m128i add(const __m128i& x, const __m128i& y)			{ return _mm_add_epi64(x, y); }
			static __forceinline __m128i sub(const __m128i& x, const __m128i& y)			{ return _mm_sub_epi64(x, y); }

			//static __forceinline __m128 { }
		};


		struct m128i_int32 {
			static __forceinline __m128i add(const __m128i& x, const __m128i& y)			{ return _mm_add_epi32(x, y); }
			static __forceinline __m128i sub(const __m128i& x, const __m128i& y)			{ return _mm_sub_epi32(x, y); }
			static __forceinline __m128i max(const __m128i& x, const __m128i& y)			{ return _mm_max_epi32(x, y); }
			static __forceinline __m128i multiply(const __m128i& x, const __m128i& y)		{ return _mm_mul_epi32(x, y); }
			//static __forceinline __m128 { }
		};

		struct m128i_uint32 {
			//static __forceinline __m128 { }

			static __forceinline __m128i max(const __m128i& x, const __m128i& y)			{ return _mm_max_epu32(x, y); }
			static __forceinline __m128i multiply(const __m128i& x, const __m128i& y)		{ return _mm_mul_epu32(x, y); }
		};

		struct m128i_int16 {
			static __forceinline __m128i add(const __m128i& x, const __m128i& y)			{ return _mm_add_epi16(x, y); }
			static __forceinline __m128i add_saturate(const __m128i& x, const __m128i& y)	{ return _mm_adds_epi16(x, y); }
			static __forceinline __m128i sub(const __m128i& x, const __m128i& y)			{ return _mm_sub_epi16(x, y); }
			static __forceinline __m128i sub_saturate(const __m128i& x, const __m128i& y)	{ return _mm_subs_epi16(x, y); }
			static __forceinline __m128i multiply_add(const __m128i& x, const __m128i& y)	{ return _mm_madd_epi16(x, y); }
			static __forceinline __m128i multiply_high(const __m128i& x, const __m128i& y)	{ return _mm_mulhi_epi16(x, y); }
			static __forceinline __m128i multiply_low(const __m128i& x, const __m128i& y)	{ return _mm_mullo_epi16(x, y); }

			static __forceinline __m128i max(const __m128i& x, const __m128i& y)			{ return _mm_max_epi16(x, y); }
			static __forceinline __m128i min(const __m128i& x, const __m128i& y)			{ return _mm_min_epi16(x, y); }

			//static __forceinline __m128i { }
		};

		struct m128i_uint16 {
			static __forceinline __m128i add_saturate(const __m128i& x, const __m128i& y)	{ return _mm_adds_epu16(x, y); }
			static __forceinline __m128i sub_saturate(const __m128i& x, const __m128i& y)	{ return _mm_subs_epu16(x, y); }
			static __forceinline __m128i average(const __m128i& x, const __m128i& y)		{ return _mm_avg_epu16(x, y); }
			static __forceinline __m128i multiply_high(const __m128i& x, const __m128i& y)	{ return _mm_mulhi_epu16(x, y); }
			static __forceinline __m128i max(const __m128i& x, const __m128i& y)			{ return _mm_max_epu16(x, y); }
			static __forceinline __m128i min(const __m128i& x, const __m128i& y)			{ return _mm_min_epu16(x, y); }

			//static __forceinline __m128i { }
		};
		
#define MM_FUNC(MM_FUNC, TYPE, X, Y) _mm##MM_FUNC##_##YPE(X, Y)
#define MM_EPI8(mm_func, x, y) MM_FUNC(add, epi8, x, y)
		struct m128i_int8 {
			static __forceinline __m128i add(const __m128i& x, const __m128i& y)			{ return _mm_add_epi8(x, y); }
			static __forceinline __m128i add_saturate(const __m128i& x, const __m128i& y)	{ return _mm_adds_epi8(x, y); }
			static __forceinline __m128i sub(const __m128i& x, const __m128i& y)			{ return _mm_sub_epi8(x, y); }
			static __forceinline __m128i sub_saturate(const __m128i& x, const __m128i& y)	{ return _mm_subs_epi8(x, y); }

			static __forceinline __m128i max(const __m128i& x, const __m128i& y)			{ return _mm_max_epi8(x, y); }
			static __forceinline __m128i min(const __m128i& x, const __m128i& y)			{ return _mm_min_epi8(x, y); }

			//static __forceinline __m128i { }
		};

		struct m128i_uint8 {
			static __forceinline __m128i add_saturate(const __m128i& x, const __m128i& y)	{ return _mm_adds_epu8(x, y); }
			static __forceinline __m128i average(const __m128i& x, const __m128i& y)		{ return _mm_avg_epu8(x, y); }
			static __forceinline __m128i sub_saturate(const __m128i& x, const __m128i& y)	{ return _mm_subs_epu8(x, y); }
			static __forceinline __m128i max(const __m128i& x, const __m128i& y)			{ return _mm_max_epu8(x, y); }
			static __forceinline __m128i min(const __m128i& x, const __m128i& y)			{ return _mm_min_epu8(x, y); }

			//static __forceinline __m128i { }
		};
	}

}}}