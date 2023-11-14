#pragma once

#include <xmmintrin.h>   // SSE
#include <emmintrin.h>   // SSE2
#include <pmmintrin.h>   // SSE3
#include <smmintrin.h>   // SSE4


#include <Utility/mask4v.h>
#include <Utility/int4v.h>
#include <Utility/float4v.h>
#include <Utility/double2v.h>

//#define USE_SSE4

namespace oyrke { namespace algorithm { namespace sse {

  /*

  ==== SSE 128 bit ====
  char16v     16 chars                       __m128i  int8x16_t       imask8x16_t
  uchar16v    16 unsigned chars              __m128i  uint8x16_t      imask8x16_t
  short8v      8 short                       __m128i  int16x8_t       imask16x8_t
  ushort8v     8 unsigned short              __m128i  uint16x8_t      imask16x8_t
  int4v        4 int                         __m128i  int32x4_t       imask32x4_t
  uint4v       4 unsigned int                __m128i  uint32x4_t      imask32x4_t
  llint2v      2 long long int               __m128i  int64x2_t       imask64x2_t
  ullint2v     2 unsigned long long int      __m128i  uint64x2_t      imask64x2_t
  128 bit signed int                         __m128i  int128x1_t      imask128x1_t
  128 bit unsigned signed int                __m128i  int128x1_t      imask128x1_t

  float32x4                                  __m128   float32x4       fmask32x4_t
  float64x2                                  __m128d  float64x2       fmask64x2_t


  ==== AVX 256 bit ====
  char16v     16 chars                       __m256i  int8x32_t       imask8x32_t
  uchar16v    16 unsigned chars              __m256i  uint8x32_t      imask8x32_t
  short8v      8 short                       __m256i  int16x16_t      imask16x16_t
  ushort8v     8 unsigned short              __m256i  uint16x16_t     imask16x16_t
  int4v        4 int                         __m256i  int32x8_t       imask32x8_t
  uint4v       4 unsigned int                __m256i  uint32x8_t      imask32x8_t
  llint2v      2 long long int               __m256i  int64x4_t       imask64x4_t
  ullint2v     2 unsigned long long int      __m256i  uint64x4_t      imask64x4_t
  128 bit signed int                         __m256i  int128x2_t      imask128x2_t
  128 bit unsigned signed int                __m256i  int128x2_t      imask128x2_t
  256 bit signed int                         __m256i  int256x1_t      imask256x1_t
  256 bit unsigned signed int                __m256i  uint256x1_t     imask256x1_t

  float32x8                                  __m256   float32x8       fmask32x8_t
  float64x4                                  __m256d  float64x4       fmask64x4_t


	float32v4
	float32x4v
	float32_v4
	float32v4_t

	fmask_vec
	imask_vec
	int_vec
	float_vec

	float_traits
	int_traits
	logic_traits

  */
#if 0

	class int32_4v {
	public:
		typedef signed int value_type;
		typedef __m128d impl_type;

		int32_4v();
		int32_4v(value_type a);
		int32_4v(value_type a0, value_type a1, value_type a2, value_type a3);
		int32_4v(const value_type* data, ALIGNED);
		int32_4v(const value_type* data, UNALIGNED);

		int32_4v& operator=(const int32_4v& rhs);
		int32_4v& operator=(value_type rhs);

		impl_type value() const;

		template <size_t N> at() const;
		value_type operator[](int ix) const;
		value_type& operator[](int ix);

		void store(value_type*, ALIGNED) const;
		void store(value_type*, UNALIGNED) const;

		template <size_t N> void store(value_type& dest);

		int32_4v& operator+=(const int32_4v& rhs);
		int32_4v& operator+=(value_type rhs);
		int32_4v& operator-=(const int32_4v& rhs);
		int32_4v& operator-=(value_type rhs);
		int32_4v& operator*=(const int32_4v& rhs);
		int32_4v& operator*=(value_type rhs);
		int32_4v& operator/=(const int32_4v& rhs);
		int32_4v& operator/=(value_type rhs);
		int32_4v& operator&=(const int32_4v& rhs);
		int32_4v& operator|=(const int32_4v& rhs);
		int32_4v& operator^=(const int32_4v& rhs);
		int32_4v& operator>>=(const int32_4v& rhs) { _mm_srl_epi64
		int32_4v& operator>>=(int count);
		int32_4v& operator<<=(const int32_4v& rhs);
		int32_4v& operator<<=(int count);


	};

  class int8_16v;  // __m128i xxx_epi8
  class uint8_16v;
  class int16_8v; // __m128i xxx_epi16
  class uint16_8v;
  class int32_4v; // __m128i xxx_epi32
  class uint32_4v;
  class int64_2v; // __m128i xxx_epi32

  // traints for 8, 16, 32 and 64 bit, modulo or saturated arithmetic
  template <typename T, typename Traits>
  class int_Nv {
    __m128i value_;

    typedef int_Nv self_t;
    typedef int_traits<T> traits;

  public:
  
    self_t operator+(const self_t& rhs) const { return traits::add(value_, rhs.value_); }
    self_t operator-(const self_t& rhs) const { return traits::subtract(value_, rhs.value_); }

    self_t operator&(const self_t& rhs) const { return traits::subtract(value_, rhs.value_); }
    self_t operator|(const self_t& rhs) const { return traits::subtract(value_, rhs.value_); }
    self_t operator^(const self_t& rhs) const { return traits::subtract(value_, rhs.value_); }
    self_t and_not(const self_t& rhs) const;

    self_t operator<<(int count) const;
    self_t operator>>(int count) const;

    // available for only a few types, free function?
    self_t average(const self_t& rhs) const;  // 8, 16 only
    Traits::add_mul_t add_multiply(const self_t& rhs) const { return traits::add_multiply(value_, rhs.value_); }

  };

#endif

}}}