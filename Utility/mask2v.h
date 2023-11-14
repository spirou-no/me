#pragma once

#include <utility/sse_internal.h>
#include <xmmintrin.h>   // SSE
#include <emmintrin.h>   // SSE2
#include <pmmintrin.h>   // SSE3
#include <smmintrin.h>   // SSE4


namespace oyrke { namespace algorithm { namespace sse {
  
  
  class mask2v {
    __m128d value_;

  public:
    __forceinline mask2v(__m128d x)                      { value_ = x; }

    __forceinline __m128d value() const { return value_; }
    // TODO C++20  proper names
    __forceinline mask2v   and20(const mask2v& rhs) const         { return _mm_and_pd(value_, rhs.value_); }
    __forceinline mask2v   or20(const mask2v& rhs) const          { return _mm_or_pd(value_,rhs.value_); }
    __forceinline mask2v   xor20(const mask2v& rhs) const         { return _mm_xor_pd(value_, rhs.value_); }
    __forceinline mask2v   and_not(const mask2v& rhs) const     { return _mm_andnot_pd(value_, rhs.value_); }

		__forceinline bool     all_true()  const { return intern::is_all_set(value_); }
    __forceinline bool     any_true()  const { return intern::is_any_set(value_); }
    __forceinline bool     all_false() const { return intern::is_all_zero(value_); }
    __forceinline bool     any_false() const { return intern::is_any_zero(value_); }
  };



}}}