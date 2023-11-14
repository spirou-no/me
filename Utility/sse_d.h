#pragma once

#include <xmmintrin.h>

#undef UNALIGNED

namespace oyrke { namespace algorithm { namespace sse {

  enum ALIGNED   { aligned };

    enum UNALIGNED { unaligned };


  class float4v {
  private:
    __m128 value_;

  public:
    __forceinline float4v(__m128 x)                                     { value_ = x; }
    __forceinline explicit float4v(float x)                             { value_ = _mm_set_ps1(x); }
    __forceinline explicit float4v(const float data[4], UNALIGNED)      { value_ = _mm_loadu_ps(data); }
    __forceinline explicit float4v(const float data[4], ALIGNED)        { value_ = _mm_load_ps(data); }

    static __forceinline float4v set_all(float x)                       { return float4v(x); }
    static __forceinline float4v create_aligned(const float data[4])    { return float4v(data, aligned); }
    static __forceinline float4v create_unaligned(const float data[4])  { return float4v(data, unaligned); }

    __forceinline __m128 value() const                                  { return value_; }
    float  operator[](size_t i) const                                   { return value_.m128_f32[i]; }
    float& operator[](size_t i)                                         { return value_.m128_f32[i]; }

    __forceinline const float4v& operator+=(const float4v& rhs)         { value_ = (*this + rhs).value_; return *this; }
    __forceinline const float4v& operator*=(const float4v& rhs)         { value_ = (*this * rhs).value_; return *this; }
    __forceinline const float4v& operator-=(const float4v& rhs)         { value_ = (*this - rhs).value_; return *this; }
    __forceinline const float4v& operator/=(const float4v& rhs)         { value_ = (*this / rhs).value_; return *this; }

    __forceinline float4v operator+(const float4v& rhs) const           { return _mm_add_ps(value_, rhs.value_); }
    __forceinline float4v operator+(float rhs) const                    { return *this + float4v(rhs); }

    __forceinline float4v operator-(const float4v& rhs) const           { return _mm_sub_ps(value_, rhs.value_); }
    __forceinline float4v operator-(float rhs) const                    { return *this - float4v(rhs); }

    __forceinline float4v operator*(const float4v& rhs) const           { return _mm_mul_ps(value_, rhs.value_); }
    __forceinline float4v operator*(float rhs) const                    { return *this * float4v(rhs); }

    __forceinline float4v operator/(const float4v& rhs) const           { return _mm_div_ps(value_, rhs.value_); }
    __forceinline float4v operator/(float rhs) const                    { return *this / float4v(rhs); }

    __forceinline float4v sqrt() const                                  { return _mm_sqrt_ps(value_); }
    __forceinline float4v reciprocal() const                            { return _mm_rcp_ps(value_); }
    __forceinline float4v sqrt_reciprocal() const                       { return _mm_rsqrt_ps(value_); }
    __forceinline float4v min(const float4v& rhs) const                 { return _mm_min_ps(value_, rhs.value_); }
    __forceinline float4v max(const float4v& rhs) const                 { return _mm_max_ps(value_, rhs.value_); }
  };

  __forceinline float4v operator+(__m128 lhs, const float4v& rhs) { return rhs + lhs; }
  __forceinline float4v operator-(__m128 lhs, const float4v& rhs) { return float4v(lhs) - rhs; }
  __forceinline float4v operator*(__m128 lhs, const float4v& rhs) { return rhs * lhs; }
  __forceinline float4v operator/(__m128 lhs, const float4v& rhs) { return float4v(lhs) / rhs; }



  //! Merges x and y, selecting 2 floats from x and 2 from y.
  // x[x0], x[x1], y[y0]m y[y1]
  template <size_t x0, size_t x1, size_t y0, size_t y1>
  __forceinline float4v
  shuffle(const float4v& x, const float4v& y) {
    return _mm_shuffle_ps(x.value(), y.value(), _MM_SHUFFLE(y1, y0, x1, x0));
  }




  //! Compute the sum of r0..3. Returned result is [sum(r0), sum(r1, sum(r2), sum(r3)]
  __forceinline float4v
  sum_row(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3) {
    float4v ab = shuffle<0,1,0,1>(r0, r1) + shuffle<2,3,2,3>(r0, r1);  // {r0_02, r0_13, r1_02, r1_13}
    float4v cd = shuffle<0,1,0,1>(r2, r3) + shuffle<2,3,2,3>(r2, r3);  // {r2_02, r2_13, r3_02, r3_13}

    float4v sum = shuffle<0,2,0,2>(ab, cd) + shuffle<1,3,1,3>(ab, cd);
    return sum;
  }



  //! Compute the sum of the 4 elements in r0 and r1.  Returned result is [sum(r0), sum(r1), x, x], x i undefined
  __forceinline float4v
  sum_row(const float4v &r0, const float4v& r1) {
    float4v tmp = shuffle<2, 3, 0, 1>(r0, r1);   // {r0_2, r0_3, r1_0, r1_1}
    float4v tmp0= r0 + tmp;  //  {r0_02, r0_13, x    , x    }
    float4v tmp1= r1 + tmp;  //  {x    , x    , r1_02, r1_13}
    tmp = shuffle<1, 1, 3, 3>(tmp0, tmp1);  // {r0_13, x, r1_13, x}
    tmp0 += tmp;  // {sum(r0), x, x, x}
    tmp1 += tmp;  // {x, x, sum(r1), x}

    tmp1 = _mm_movehl_ps(tmp1.value(), tmp1.value());  // {sum(r1), x, sum(r1), x}
    return _mm_unpacklo_ps(tmp0.value(), tmp1.value());  // {sum(r0), sum(r1), x, x}
  }

  //! Compute the sum of the 4 elements in r.  Returned result is [sum(r), x, x, x], x i undefined
  __forceinline float4v
  sum_row(const float4v &r) {
    float4v tmp = r.value() + _mm_movehl_ps(r.value(), r.value());    // { r_02, r_13, x, x} ignore upper 2 floats
    float4v tmp2 = shuffle<1, 0, 0, 1>(tmp, tmp); //_mm_set_ss(tmp.m128_f32[1]);
    return _mm_add_ss(tmp.value(), tmp2.value());
  }

}}}