/*


15.24   bicubic_spline_unrolled_sse_v2(float& xin, float& yin, const float xknots[16], const float yknots[16])

*/
#include <utility/performance_testing.h>
#include <utility/sse.h>
#include <utility/tiny_algo.h>
#include <time.h>
#include <functional>
#include <numeric>

#if 0
struct coord2d {
  double x, y;

  coord2d() : x(0.0), y(0.0) {}
  coord2d(double xx, double yy) : x(xx), y(yy) {}
};

__forceinline coord2d
operator+(const coord2d& lhs, const coord2d& rhs) {
  return coord2d(lhs.x + rhs.x, lhs.y + rhs.y);
}

__forceinline coord2d
operator*(const double lhs, const coord2d& rhs) {
  return coord2d(lhs * rhs.x, lhs * rhs.y);
}

__forceinline coord2d
operator*(const coord2d& lhs, const double rhs) {
  return rhs * lhs;

}

bicubic spline (coord2d)             :      2038.891 ms    13518195      150.83 ns
bicubic spline (coord2d, recompute)  :      1985.403 ms    12047300      164.80 ns
bicubic spline (coord2d opt)         :      1995.724 ms    21671767       92.09 ns
bicubic spline (coord2d, recompute)  :      2012.198 ms    16797403      119.79 ns

#else
//bicubic spline (coord2d)             :      1986.281 ms    14265874      139.23 ns
//bicubic spline (coord2d, recompute)  :      2031.920 ms    14455556      140.56 ns
//bicubic spline (coord2d opt)         :      1999.608 ms    21126385       94.65 ns
//bicubic spline (coord2d, recompute)  :      2001.522 ms    16070179      124.55 ns

struct coord2d {
  typedef oyrke::algorithm::sse::double2v impl_t;
  impl_t xy;

  coord2d() {}
  coord2d(const impl_t& values) : xy(values) {}
  coord2d(double xx, double yy) : xy(xx, yy) {}

  __forceinline double x() const { return xy[0]; }
};

__forceinline coord2d
operator+(const coord2d& lhs, const coord2d& rhs) {
  return coord2d(lhs.xy + rhs.xy);
}

__forceinline coord2d
operator*(const coord2d& lhs, const coord2d& rhs) {
  return rhs.xy * lhs.xy;
}

__forceinline coord2d
operator*(const double lhs, const coord2d& rhs) {
  return coord2d(lhs * rhs.xy);
}

__forceinline coord2d
operator*(const coord2d& lhs, const double rhs) {
  return rhs * lhs.xy;
}

#endif

template <typename T>
__forceinline 
T weighted_sum(T x0, T w0, T x1, T w1) {
  return x0 * w0  +  x1 * w1;
}

template <typename T>
inline 
T linear_interpolate(T a, T b, T x) {
  return weighted_sum(a, 1-x, b, x);
}

template <typename T>
inline 
T linear_interpolate2(T a, T b, T x) {
  return a*(1-x) + b*x;
}

template <typename T>
inline 
T linear_interpolate_generic(const T *data, const T *x) {
  return data[0]*(1-x[0]) + data[1]*x[0];
}


template <typename T>
inline 
T bilinear_interpolate_plain(T f00, T f10, T f01, T f11, T x, T y) {
  return f00*(1-x)*(1-y) + f10*x*(1-y) + f01*(1-x)*y + f11*x*y;
}

template <typename T>
inline 
T bilinear_interpolate_opt(T f00, T f10, T f01, T f11, T x, T y) {
  T s1 = weighted_sum(f00, 1-x, f10, x);
  T s2 = weighted_sum(f01, 1-x, f11, x);
  
  return linear_interpolate(s1, s2, y);
}

template <typename T, int N>
inline 
T linear_interpolate(const T data[], const T x[N]) {
  
  T wsum[1<<N];
  std::copy(data, data+(1<<N), wsum);
  for (int dim = 0; dim < N; ++dim) {
    T w1 =  1-x[dim];
    int nof_samples = 1<<(N-dim-1);
    for (int i = 0; i < nof_samples; ++i) {
      wsum[i] = weighted_sum(wsum[2*i], w1, wsum[2*i+1], x[dim]);
    }
  }

  return wsum[0];
}

static const double one_sixth = 1.0/6.0;

template <typename T>
__forceinline T spline_0(const T&x) { 
  return static_cast<T>(one_sixth*(((-x + 3.0)*x -3.0)*x + 1.0));
}

template <typename T>
__forceinline T spline_1(const T&x) { 
  return static_cast<T>(one_sixth * ((3.0*x -6.0)*x*x + 4.0)); 
}

template <typename T>
__forceinline T spline_2(const T&x) { 
  return static_cast<T>(one_sixth * (((-x +1.0)*x + 1.0)*x*3.0 + 1.0)); 
  
}
template <typename T>
__forceinline T spline_3(const T&x) { 
  return static_cast<T>(one_sixth) * x*x*x; 
}

template <typename T>
__forceinline T spline_opt_0(const T&x) { 
  //return ((3-x)*x -3)*x + 1;
  T x3 = x+x+x;
  T a = x3 - x*x;
  return a*x - x3 + 1;
  //return ((3-x)*x -3)*x + 1;
}

template <typename T>
__forceinline T spline_opt_1(const T&x) { 
  //return 3 * ((x -2)*x + 1)*x;  
  //return (x-2)*x*x*3 + 4;
  return (x*x-(x+x))*x*3 + 4;
}

template <typename T>
__forceinline T spline_opt_2(const T&x) { 
  //return 3 * (-x*x + 1)*x; 
  //return 3*((1-x)*x + 1)*x + 1;
  
  return 3*((x-x*x)*x + x) + 1; 
}

template <typename T>
__forceinline T spline_opt_3(const T&x) {
  return x*x*x; 
}


template <typename T, typename V>
__forceinline 
V cubic_spline(const T& x, const V knots[4]) {
  return knots[0] * spline_0(x)
       + knots[1] * spline_1(x)
       + knots[2] * spline_2(x)
       + knots[3] * spline_3(x);
}


#include <xmmintrin.h>
__forceinline
float 
cubic_spline_sse(float xin, const float knots[4]) {
  /*
      -1  3  -3  1
       3 -6   3  0
      -3  0   3  0
      -1  4   1  0
  */

  //((a*x + b)*x + c)*x + d
  __declspec(align(16)) static float a_coeffs[] = {-1.0f/6.0f,  3.0f/6.0f,-3.0f/6.0f,  1.0f/6.0f};
  __declspec(align(16)) static float b_coeffs[] = { 3.0f/6.0f, -1.0f     ,      0.0f,  4.0f/6.0f};
  __declspec(align(16)) static float c_coeffs[] = {-3.0f/6.0f,  3.0f/6.0f, 3.0f/6.0f,  1.0f/6.0f};
  __declspec(align(16)) static float d_coeffs[] = { 1.0f/6.0f, 0.0f,       0.0f,       0.0f};

  __m128 x = _mm_set_ps1(xin);  // { xin, xin, xin, xin}
  __m128 a = _mm_load_ps(a_coeffs);
  __m128 sum = _mm_mul_ps(a, x);  // a*x
  __m128 b = _mm_load_ps(b_coeffs);
  sum = _mm_add_ps(sum, b);       // a*x + b
  sum = _mm_mul_ps(sum, x);       // (a*x + b)*x
  __m128 c = _mm_load_ps(c_coeffs);
  sum = _mm_add_ps(sum, c);          // (a*x + b)*x + c
  sum = _mm_mul_ps(sum, x);       // ((a*x + b)*x + c)*x
  __m128 d = _mm_load_ps(d_coeffs);
  sum = _mm_add_ps(sum, d);          // ((a*x + b)*x + c)*x + d

  // sum now holds the 4 cubic B-spline coeffisients to weigh the knots with
  __m128 k = _mm_load_ps(knots);
  sum = _mm_mul_ps(sum, k);   // each knot weigthed by B-spline factor

  __m128 sum2 = _mm_movehl_ps(sum, sum); // {a3, a2, a3, a2}
  sum = _mm_add_ps(sum, sum2);    // ignore upper 2 floats
  sum2= _mm_set_ss(sum.m128_f32[1]);
  sum = _mm_add_ss(sum, sum2);
  return sum.m128_f32[0];
}

__forceinline
float 
cubic_spline_sse_unaligned(float xin, const float knots[4]) {
  /*
      -1  3  -3  1
       3 -6   3  0
      -3  0   3  0
      -1  4   1  0
  */

  //((a*x + b)*x + c)*x + d
  static float a_coeffs[] = {-1.0f/6.0f,  3.0f/6.0f,-3.0f/6.0f,  1.0f/6.0f};
  static float b_coeffs[] = { 3.0f/6.0f, -1.0f     ,      0.0f,  4.0f/6.0f};
  static float c_coeffs[] = {-3.0f/6.0f,  3.0f/6.0f, 3.0f/6.0f,  1.0f/6.0f};
  static float d_coeffs[] = { 1.0f/6.0f,  0.0f,           0.0f,  0.0f     };

  __m128 x = _mm_set_ps1(xin);  // { xin, xin, xin, xin}
  __m128 a = _mm_loadu_ps(a_coeffs);
  __m128 sum = _mm_mul_ps(a, x);  // a*x
  __m128 b = _mm_loadu_ps(b_coeffs);
  sum = _mm_add_ps(sum, b);       // a*x + b
  sum = _mm_mul_ps(sum, x);       // (a*x + b)*x
  __m128 c = _mm_loadu_ps(c_coeffs);
  sum = _mm_add_ps(sum, c);          // (a*x + b)*x + c
  sum = _mm_mul_ps(sum, x);       // ((a*x + b)*x + c)*x
  __m128 d = _mm_loadu_ps(d_coeffs);
  sum = _mm_add_ps(sum, d);          // ((a*x + b)*x + c)*x + d

  // sum now holds the 4 cubic B-spline coeffisients to weigh the knots with
  __m128 k = _mm_loadu_ps(knots);
  sum = _mm_mul_ps(sum, k);   // each knot weigthed by B-spline factor

  __m128 sum2 = _mm_movehl_ps(sum, sum); // {a3, a2, a3, a2}
  sum = _mm_add_ps(sum, sum2);    // ignore upper 2 floats
  sum2= _mm_set_ss(sum.m128_f32[1]);
  sum = _mm_add_ss(sum, sum2);
  return sum.m128_f32[0];
}

__forceinline
float 
cubic_spline_sse_v2(float xin, const float knots[4]) {
  /*
      -1  3  -3  1
       3 -6   3  0
      -3  0   3  0
      -1  4   1  0
  */

  //((a*x + b)*x + c)*x + d
  static const __m128 a = {-1.0f/6.0f,  3.0f/6.0f,-3.0f/6.0f,  1.0f/6.0f};
  static const __m128 b = { 3.0f/6.0f, -1.0f     ,      0.0f,  4.0f/6.0f};
  static const __m128 c = {-3.0f/6.0f,  3.0f/6.0f, 3.0f/6.0f,  1.0f/6.0f};
  static const __m128 d = { 1.0f/6.0f,  0.0f,           0.0f,  0.0f};

  __m128 x = _mm_set_ps1(xin);  // { xin, xin, xin, xin}
  __m128 sum = _mm_mul_ps(a, x);  // a*x
  sum = _mm_add_ps(sum, b);       // a*x + b
  sum = _mm_mul_ps(sum, x);       // (a*x + b)*x
  sum = _mm_add_ps(sum, c);          // (a*x + b)*x + c
  sum = _mm_mul_ps(sum, x);       // ((a*x + b)*x + c)*x
  sum = _mm_add_ps(sum, d);          // ((a*x + b)*x + c)*x + d

  __m128 k = _mm_load_ps(knots);
  sum = _mm_mul_ps(sum, k);   // each knot weigthed by B-spline factor

  __m128 sum2 = _mm_movehl_ps(sum, sum); // {a3, a2, a3, a2}
  sum = _mm_add_ps(sum, sum2);    // ignore upper 2 floats
  return sum.m128_f32[0] + sum.m128_f32[1];
}


__forceinline float 
cubic_spline_sse_v4(float xin, const float knots[4]) {
  /*
      -1  3  -3  1
       3 -6   3  0
      -3  0   3  0
      -1  4   1  0
  */

  //((a*x + b)*x + c)*x + d
  static const __m128 a = {-1.0f/6.0f,  3.0f/6.0f,-3.0f/6.0f, -1.0f/6.0f};
  static const __m128 b = { 3.0f/6.0f, -1.0f     ,      0.0f,  4.0f/6.0f};
  static const __m128 c = {-3.0f/6.0f,  3.0f/6.0f, 3.0f/6.0f,  1.0f/6.0f};
  static const __m128 d = { 1.0f/6.0f, 0.0f,       0.0f, 0.0f};

  namespace sse = oyrke::algorithm::sse;
  sse::float4v x(xin);
  sse::float4v result = ((x*a + b)*x + c)*x + d;
  return sse::sum_row(result)[0];
}


__forceinline
oyrke::algorithm::sse::float4v
compute_spline_weigths(float xin) {
  static const __m128 a = {-1.0f/6.0f,  3.0f/6.0f,-3.0f/6.0f, -1.0f/6.0f};
  static const __m128 b = { 3.0f/6.0f, -1.0f     ,      0.0f,  4.0f/6.0f};
  static const __m128 c = {-3.0f/6.0f,  3.0f/6.0f, 3.0f/6.0f,  1.0f/6.0f};
  static const __m128 d = { 1.0f/6.0f, 0.0f,       0.0f, 0.0f};

  namespace sse = oyrke::algorithm::sse;
  sse::float4v x(xin);
  sse::float4v result = ((x*a + b)*x + c)*x + d;
  //__m128 x = _mm_set_ps1(xin);  // { xin, xin, xin, xin}
  //__m128 sum = _mm_mul_ps(a, x);  // a*x
  //sum = _mm_add_ps(sum, b);       // a*x + b
  //sum = _mm_mul_ps(sum, x);       // (a*x + b)*x
  //sum = _mm_add_ps(sum, c);          // (a*x + b)*x + c
  //sum = _mm_mul_ps(sum, x);       // ((a*x + b)*x + c)*x
  //sum = _mm_add_ps(sum, d);          // ((a*x + b)*x + c)*x + d
  
  return result;
}


float
bicubic_spline_sse(float x, float y, const float knots[16]) {
  namespace sse = oyrke::algorithm::sse;
  sse::float4v xw = compute_spline_weigths(x);

  sse::float4v y0 = xw * sse::float4v(&knots[ 0], sse::unaligned);
  sse::float4v y1 = xw * sse::float4v(&knots[ 4], sse::unaligned);
  sse::float4v y2 = xw * sse::float4v(&knots[ 8], sse::unaligned);
  sse::float4v y3 = xw * sse::float4v(&knots[12], sse::unaligned);

  sse::float4v ysum = compute_spline_weigths(y) * sse::sum_row(y0, y1, y2, y3);
  //sse::float4v yw = compute_spline_weigths(y);
  //ysum *= compute_spline_weigths(y);
  return sse::sum_row(ysum)[0];
}

template <typename Coord, typename Value>
Value
bicubic_spline(const Coord& x, const Coord& y, const Value knots[16]) {
  // 1st compute spline weigths along x dimension
  Coord xw[4];
  xw[0] = spline_0(x);
  xw[1] = spline_1(x);
  xw[2] = spline_2(x);
  xw[3] = spline_3(x);

  Value yw[4];
  //yw[0] = std::inner_product(xw, xw+4, &knots[ 0], Value());
  //yw[1] = std::inner_product(xw, xw+4, &knots[ 4], Value());
  //yw[2] = std::inner_product(xw, xw+4, &knots[ 8], Value());
  //yw[3] = std::inner_product(xw, xw+4, &knots[12], Value());
  using namespace oyrke::algorithm;

  yw[0] = tiny_algo<4>::inner_product(xw, &knots[ 0], Value());
  yw[1] = tiny_algo<4>::inner_product(xw, &knots[ 4], Value());
  yw[2] = tiny_algo<4>::inner_product(xw, &knots[ 8], Value());
  yw[3] = tiny_algo<4>::inner_product(xw, &knots[12], Value());
  Value result = cubic_spline(y, yw);
  return result;
}

template <typename T>
coord2d
bicubic_spline(const T& x, const T& y, const T xknots[16], const T yknots[16]) {
  // 1st compute spline weigths along x dimension
  T xw[4];
  xw[0] = spline_0(x);
  xw[1] = spline_1(x);
  xw[2] = spline_2(x);
  xw[3] = spline_3(x);

  T yw[4];
  T yy[4];
  //yw[0] = std::inner_product(xw, xw+4, &knots[ 0], Value());
  //yw[1] = std::inner_product(xw, xw+4, &knots[ 4], Value());
  //yw[2] = std::inner_product(xw, xw+4, &knots[ 8], Value());
  //yw[3] = std::inner_product(xw, xw+4, &knots[12], Value());
  using namespace oyrke::algorithm;

  yw[0] = tiny_algo<4>::inner_product(xw, &xknots[ 0], T());
  yy[0] = tiny_algo<4>::inner_product(xw, &yknots[ 0], T());
  yw[1] = tiny_algo<4>::inner_product(xw, &xknots[ 4], T());
  yy[1] = tiny_algo<4>::inner_product(xw, &yknots[ 4], T());
  yw[2] = tiny_algo<4>::inner_product(xw, &xknots[ 8], T());
  yy[2] = tiny_algo<4>::inner_product(xw, &yknots[ 8], T());
  yw[3] = tiny_algo<4>::inner_product(xw, &xknots[12], T());
  yy[3] = tiny_algo<4>::inner_product(xw, &yknots[12], T());

  T xres = cubic_spline(y, yw);
  T yres = cubic_spline(y, yy);
  coord2d result(xres, yres);
  return result;
}

template <typename Coord, typename Value>
Value
bicubic_spline_unrolled(const Coord& x, const Coord& y, const Value knots[16]) {
  // 1st compute spline weigths along x dimension
  static Coord rescale=Coord(1)/Coord(36);
  Coord xw_0 = spline_opt_0(x);
  Coord xw_1 = spline_opt_1(x);
  Coord xw_2 = spline_opt_2(x);
  Coord xw_3 = spline_opt_3(x);

  Value yw_0 = xw_0*knots[0+ 0] + xw_1*knots[1+ 0] + xw_2*knots[2+ 0] + xw_3*knots[3+ 0];
  Value yw_1 = xw_0*knots[0+ 4] + xw_1*knots[1+ 4] + xw_2*knots[2+ 4] + xw_3*knots[3+ 4];
  Value yw_2 = xw_0*knots[0+ 8] + xw_1*knots[1+ 8] + xw_2*knots[2+ 8] + xw_3*knots[3+ 8];
  Value yw_3 = xw_0*knots[0+12] + xw_1*knots[1+12] + xw_2*knots[2+12] + xw_3*knots[3+12];
  Value result = yw_0*spline_opt_0(y) + yw_1*spline_opt_1(y) + yw_2*spline_opt_2(y) + yw_3*spline_opt_3(y);
  return result*rescale;
}


template <typename T>
coord2d
bicubic_spline_unrolled(const T& x, const T& y, const T xknots[16], const T yknots[16]) {
  // 1st compute spline weigths along x dimension
  static T rescale=T(1)/T(36);

  T xw_0 = spline_opt_0(x);
  T xw_1 = spline_opt_1(x);
  T xw_2 = spline_opt_2(x);
  T xw_3 = spline_opt_3(x);

  T yx_0 = xw_0 * xknots[0+ 0] + xw_1 *xknots[1+ 0] + xw_2 * xknots[2+ 0] + xw_3 * xknots[3+ 0];
  T yy_0 = xw_0 * yknots[0+ 0] + xw_1 *yknots[1+ 0] + xw_2 * yknots[2+ 0] + xw_3 * yknots[3+ 0];
  T yx_1 = xw_0 * xknots[0+ 4] + xw_1 *xknots[1+ 4] + xw_2 * xknots[2+ 4] + xw_3 * xknots[3+ 4];
  T yy_1 = xw_0 * yknots[0+ 4] + xw_1 *yknots[1+ 4] + xw_2 * yknots[2+ 4] + xw_3 * yknots[3+ 4];
  T yx_2 = xw_0 * xknots[0+ 8] + xw_1 *xknots[1+ 8] + xw_2 * xknots[2+ 8] + xw_3 * xknots[3+ 8];
  T yy_2 = xw_0 * yknots[0+ 8] + xw_1 *yknots[1+ 8] + xw_2 * yknots[2+ 8] + xw_3 * yknots[3+ 8];
  T yx_3 = xw_0 * xknots[0+12] + xw_1 *xknots[1+12] + xw_2 * xknots[2+12] + xw_3 * xknots[3+12];
  T yy_3 = xw_0 * yknots[0+12] + xw_1 *yknots[1+12] + xw_2 * yknots[2+12] + xw_3 * yknots[3+12];

  T yw_0 = spline_opt_0(y);
  T yw_1 = spline_opt_1(y);
  T yw_2 = spline_opt_2(y);
  T yw_3 = spline_opt_3(y);

  coord2d result(yx_0*yw_0 + yx_1*yw_1 + yx_2*yw_2 + yx_3*yw_3, 
                 yy_0*yw_0 + yy_1*yw_1 + yy_2*yw_2 + yy_3*yw_3);
  return result*rescale;
}

coord2d
bicubic_spline_unrolled_sse(double x, double y, const double xknots[16], const double yknots[16]) {

  namespace sse = oyrke::algorithm::sse;

  sse::double2v xy(x, y); 

  sse::double2v w_0 = spline_opt_0(xy);
  sse::double2v w_1 = spline_opt_1(xy);
  sse::double2v w_2 = spline_opt_2(xy);
  sse::double2v w_3 = spline_opt_3(xy);

  sse::double2v xw_01 = sse::double2v::compose_00(w_0, w_1);
  sse::double2v xw_23 = sse::double2v::compose_00(w_2, w_3);
  sse::double2v yw_01 = sse::double2v::compose_11(w_0, w_1);
  sse::double2v yw_23 = sse::double2v::compose_11(w_2, w_3);

  sse::double2v yx_0 = xw_01 * sse::double2v(xknots+ 0, sse::unaligned) + xw_23 * sse::double2v(xknots+ 2, sse::unaligned);
  sse::double2v yx_1 = xw_01 * sse::double2v(xknots+ 4, sse::unaligned) + xw_23 * sse::double2v(xknots+ 6, sse::unaligned);
  sse::double2v yx_01 = sse::sum_row(yx_0, yx_1);

  sse::double2v yx_2 = xw_01 * sse::double2v(xknots+ 8, sse::unaligned) + xw_23 * sse::double2v(xknots+10, sse::unaligned);
  sse::double2v yx_3 = xw_01 * sse::double2v(xknots+12, sse::unaligned) + xw_23 * sse::double2v(xknots+14, sse::unaligned);
  sse::double2v yx_23 = sse::sum_row(yx_2, yx_3);

  sse::double2v yy_0 = xw_01 * sse::double2v(yknots+ 0, sse::unaligned) + xw_23 * sse::double2v(yknots+ 2, sse::unaligned);
  sse::double2v yy_1 = xw_01 * sse::double2v(yknots+ 4, sse::unaligned) + xw_23 * sse::double2v(yknots+ 6, sse::unaligned);
  sse::double2v yy_01 = sse::sum_row(yy_0, yy_1);

  sse::double2v yy_2 = xw_01 * sse::double2v(yknots+ 8, sse::unaligned) + xw_23 * sse::double2v(yknots+10, sse::unaligned);
  sse::double2v yy_3 = xw_01 * sse::double2v(yknots+12, sse::unaligned) + xw_23 * sse::double2v(yknots+14, sse::unaligned);
  sse::double2v yy_23 = sse::sum_row(yy_2, yy_3);

  sse::double2v xres = yx_01 * yw_01 + yx_23 * yw_23;
  sse::double2v yres = yy_01 * yw_01 + yy_23 * yw_23;

  sse::double2v xyres = sse::double2v(one_sixth*one_sixth) * sse::sum_row(xres, yres);
  return coord2d(xyres.value());
}

coord2d
bicubic_spline_unrolled_sse_v2(double xin, double yin, const double xknots[16], const double yknots[16]) {

  namespace sse = oyrke::algorithm::sse;

    //((a*x + b)*x + c)*x + d
  static const __m128d a01 = {-1.0f/6.0f,  3.0f/6.0f};
  static const __m128d a23 = {-3.0f/6.0f, -1.0f/6.0f};
  static const __m128d b01 = { 3.0f/6.0f, -1.0f     };
  static const __m128d b23 = {      0.0f,  4.0f/6.0f};
  static const __m128d c01 = {-3.0f/6.0f,  3.0f/6.0f};
  static const __m128d c23 = { 3.0f/6.0f,  1.0f/6.0f};
  static const __m128d d01 = { 1.0f/6.0f, 0.0f };
  static const __m128d d23 = {       0.0f, 0.0f};


  sse::double2v x(xin);
  sse::double2v y(yin); 

  sse::double2v xw_01 = ((a01*x + b01)*x + c01)*x + d01;
  sse::double2v yw_01 = ((a01*y + b01)*y + c01)*y + d01;
  sse::double2v xw_23 = ((a23*x + b23)*x + c23)*x; // d23 is 0,0
  sse::double2v yw_23 = ((a23*y + b23)*y + c23)*y; // d23 is 0,0

  sse::double2v yx_0 = xw_01 * sse::double2v(xknots+ 0, sse::unaligned) + xw_23 * sse::double2v(xknots+ 2, sse::unaligned);
  sse::double2v yx_1 = xw_01 * sse::double2v(xknots+ 4, sse::unaligned) + xw_23 * sse::double2v(xknots+ 6, sse::unaligned);
  sse::double2v yx_01 = sse::sum_row(yx_0, yx_1);

  sse::double2v yx_2 = xw_01 * sse::double2v(xknots+ 8, sse::unaligned) + xw_23 * sse::double2v(xknots+10, sse::unaligned);
  sse::double2v yx_3 = xw_01 * sse::double2v(xknots+12, sse::unaligned) + xw_23 * sse::double2v(xknots+14, sse::unaligned);
  sse::double2v yx_23 = sse::sum_row(yx_2, yx_3);

  sse::double2v yy_0 = xw_01 * sse::double2v(yknots+ 0, sse::unaligned) + xw_23 * sse::double2v(yknots+ 2, sse::unaligned);
  sse::double2v yy_1 = xw_01 * sse::double2v(yknots+ 4, sse::unaligned) + xw_23 * sse::double2v(yknots+ 6, sse::unaligned);
  sse::double2v yy_01 = sse::sum_row(yy_0, yy_1);

  sse::double2v yy_2 = xw_01 * sse::double2v(yknots+ 8, sse::unaligned) + xw_23 * sse::double2v(yknots+10, sse::unaligned);
  sse::double2v yy_3 = xw_01 * sse::double2v(yknots+12, sse::unaligned) + xw_23 * sse::double2v(yknots+14, sse::unaligned);
  sse::double2v yy_23 = sse::sum_row(yy_2, yy_3);

  sse::double2v xres = yx_01 * yw_01 + yx_23 * yw_23;
  sse::double2v yres = yy_01 * yw_01 + yy_23 * yw_23;

  sse::double2v xyres = sse::sum_row(xres, yres);
  return coord2d(xyres.value());
}

void
bicubic_spline_unrolled_sse(float& x, float& y, const float xknots[16], const float yknots[16]) {

  namespace sse = oyrke::algorithm::sse;

  sse::float4v xin(x); 
  sse::float4v yin(y); 
  sse::float4v xy = sse::shuffle<0,0,0,0>(xin, yin);  // x,x,y,y

  sse::float4v w_0 = spline_opt_0(xy); 
  sse::float4v w_1 = spline_opt_1(xy);
  sse::float4v w_2 = spline_opt_2(xy);
  sse::float4v w_3 = spline_opt_3(xy);

  //sse::float4v xw_01 = sse::shuffle<0,0,0,0>(w_0, w_1);
  //sse::float4v xw_23 = sse::shuffle<0,0,0,0>(w_2, w_3);
  sse::float4v w_01 = sse::shuffle<0,2,0,2>(w_0, w_1);
  sse::float4v w_23 = sse::shuffle<0,2,0,2>(w_2, w_3);
  sse::float4v xw = sse::shuffle<0,2,0,2>(w_01, w_23);

  //sse::float4v yw_01 = sse::shuffle<2,2,2,2>(w_0, w_1);
  //sse::float4v yw_23 = sse::shuffle<2,2,2,2>(w_2, w_3);
  sse::float4v yw = sse::shuffle<1,3,1,3>(w_01, w_23);

  sse::float4v yx_0 = xw * sse::float4v(xknots+ 0, sse::unaligned);
  sse::float4v yx_1 = xw * sse::float4v(xknots+ 4, sse::unaligned);
  sse::float4v yx_2 = xw * sse::float4v(xknots+ 8, sse::unaligned);
  sse::float4v yx_3 = xw * sse::float4v(xknots+12, sse::unaligned);
  sse::float4v yx = sse::sum_row(yx_0, yx_1, yx_2, yx_3);

  sse::float4v yy_0 = xw * sse::float4v(yknots+ 0, sse::unaligned);
  sse::float4v yy_1 = xw * sse::float4v(yknots+ 4, sse::unaligned);
  sse::float4v yy_2 = xw * sse::float4v(yknots+ 8, sse::unaligned);
  sse::float4v yy_3 = xw * sse::float4v(yknots+12, sse::unaligned);
  sse::float4v yy = sse::sum_row(yy_0, yy_1, yy_2, yy_3);

  sse::float4v xres = yw * yx;
  sse::float4v yres = yw * yy;

  sse::float4v xyres = sse::sum_row(xres, yres) * float(one_sixth*one_sixth);
  //xyres.store0(&x);
  //sse::shuffle<1,1,1,1>(xyres, xyres).store0(&y);
	//return coord2d(xyres.value());
	xyres.store<0>(x);
	xyres.store<1>(y);
	//sse::shuffle<1,1,1,1>(xyres,xyres).store0(&y);
	//sse::float4v::create_duplicated<1>(xyres).store0(&y);
  //x = xyres.at<0>();
  //y = xyres.at<1>();
}



__forceinline void
bicubic_spline_unrolled_sse_v2(float& xin, float& yin, const float xknots[16], const float yknots[16]) {
  /*
      -1  3  -3  1
       3 -6   3  0
      -3  0   3  0
      -1  4   1  0
  */

  //((a*x + b)*x + c)*x + d
  static const __m128 a = {-1.0f/6.0f,  3.0f/6.0f,-3.0f/6.0f,  1.0f/6.0f};
  static const __m128 b = { 3.0f/6.0f, -1.0f     ,      0.0f,  4.0f/6.0f};
  static const __m128 c = {-3.0f/6.0f,  3.0f/6.0f, 3.0f/6.0f,  1.0f/6.0f};
  static const __m128 d = { 1.0f/6.0f, 0.0f,       0.0f, 0.0f};

  namespace sse = oyrke::algorithm::sse;
  sse::float4v x(xin);

  sse::float4v xw = ((x*a + b)*x + c)*x + d;

  sse::float4v yx_0 = xw * sse::float4v(xknots+ 0, sse::unaligned);
  sse::float4v yx_1 = xw * sse::float4v(xknots+ 4, sse::unaligned);
  sse::float4v yx_2 = xw * sse::float4v(xknots+ 8, sse::unaligned);
  sse::float4v yx_3 = xw * sse::float4v(xknots+12, sse::unaligned);
  sse::float4v yx = sse::sum_row(yx_0, yx_1, yx_2, yx_3);

  sse::float4v yy_0 = xw * sse::float4v(yknots+ 0, sse::unaligned);
  sse::float4v yy_1 = xw * sse::float4v(yknots+ 4, sse::unaligned);
  sse::float4v yy_2 = xw * sse::float4v(yknots+ 8, sse::unaligned);
  sse::float4v yy_3 = xw * sse::float4v(yknots+12, sse::unaligned);
  sse::float4v yy = sse::sum_row(yy_0, yy_1, yy_2, yy_3);

  sse::float4v y(yin);
  sse::float4v yw = ((y*a + b)*y + c)*y + d;

  sse::float4v xres = yw * yx;
  sse::float4v yres = yw * yy;

  sse::float4v xyres = sse::sum_row(xres, yres);
  //xin = xyres.at<0>(); 
	//yin = xyres.at<1>();
	xyres.store<0>(xin);
	xyres.store<1>(yin);
	//yin = xyres.at<1>();
	//sse::float4v(xyres.at<1>()).store0(&yin);
  //sse::shuffle<1,1,1,1>(xyres, xyres).store0(&yin);
	//sse::float4v::create_duplicated<1>(xyres).store0(&yin);
  //xin = xres.at<0>() + xres.at<1>() + xres.at<2>() + xres.at<3>();
  //yin = yres.at<0>() + yres.at<1>() + yres.at<2>() + yres.at<3>();
}

__forceinline void
bicubic_spline_unrolled_sse_v3(float& xin, float& yin, const float xknots[16], const float yknots[16]) {
  /*
      -1  3  -3  1
       3 -6   3  0
      -3  0   3  0
      -1  4   1  0
  */

  //((a*x + b)*x + c)*x + d
  static const __m128 a = {-1.0f/6.0f,  3.0f/6.0f,-3.0f/6.0f,  1.0f/6.0f};
  static const __m128 b = { 3.0f/6.0f, -1.0f     ,      0.0f,  4.0f/6.0f};
  static const __m128 c = {-3.0f/6.0f,  3.0f/6.0f, 3.0f/6.0f,  1.0f/6.0f};
  static const __m128 d = { 1.0f/6.0f, 0.0f,       0.0f, 0.0f};

  namespace sse = oyrke::algorithm::sse;
  sse::float4v x(xin);

	sse::float4v kx0 = sse::float4v(xknots+ 0, sse::unaligned);
	sse::float4v kx1 = sse::float4v(xknots+ 4, sse::unaligned);
	sse::float4v kx2 = sse::float4v(xknots+ 8, sse::unaligned);
	sse::float4v kx3 = sse::float4v(xknots+12, sse::unaligned);

  sse::float4v xw = ((x*a + b)*x + c)*x + d;

  sse::float4v yx_0 = xw * kx0;
  sse::float4v yx_1 = xw * kx1;
  sse::float4v yx_2 = xw * kx2;
  sse::float4v yx_3 = xw * kx3;
  sse::float4v yx = sse::sum_row(yx_0, yx_1, yx_2, yx_3);

	sse::float4v ky0 = sse::float4v(yknots+ 0, sse::unaligned);
	sse::float4v ky1 = sse::float4v(yknots+ 4, sse::unaligned);
	sse::float4v ky2 = sse::float4v(yknots+ 8, sse::unaligned);
	sse::float4v ky3 = sse::float4v(yknots+12, sse::unaligned);

  sse::float4v y(yin);
  sse::float4v yw = ((y*a + b)*y + c)*y + d;

	sse::float4v yy_0 = xw * ky0;
  sse::float4v yy_1 = xw * ky1;
  sse::float4v yy_2 = xw * ky2;
  sse::float4v yy_3 = xw * ky3;

	sse::float4v yy = sse::sum_row(yy_0, yy_1, yy_2, yy_3);

  sse::float4v xres = yw * yx;
  sse::float4v yres = yw * yy;

  sse::float4v xyres = sse::sum_row(xres, yres);
  //xin = xyres.at<0>(); 
	//yin = xyres.at<1>();
	xyres.store<0>(xin);
	xyres.store<1>(yin);
	//yin = xyres.at<1>();
	//sse::float4v(xyres.at<1>()).store0(&yin);
  //sse::shuffle<1,1,1,1>(xyres, xyres).store0(&yin);
	//sse::float4v::create_duplicated<1>(xyres).store0(&yin);
  //xin = xres.at<0>() + xres.at<1>() + xres.at<2>() + xres.at<3>();
  //yin = yres.at<0>() + yres.at<1>() + yres.at<2>() + yres.at<3>();
}


__forceinline void
bicubic_spline_unrolled_sse4(float& xin, float& yin, const float xknots[16], const float yknots[16]) {
  /*
      -1  3  -3  1
       3 -6   3  0
      -3  0   3  0
      -1  4   1  0
  */

  //((a*x + b)*x + c)*x + d
  static const __m128 a = {-1.0f/6.0f,  3.0f/6.0f,-3.0f/6.0f,  1.0f/6.0f};
  static const __m128 b = { 3.0f/6.0f, -1.0f     ,      0.0f,  4.0f/6.0f};
  static const __m128 c = {-3.0f/6.0f,  3.0f/6.0f, 3.0f/6.0f,  1.0f/6.0f};
  static const __m128 d = { 1.0f/6.0f, 0.0f,       0.0f, 0.0f};

  namespace sse = oyrke::algorithm::sse;
  sse::float4v x(xin);
  sse::float4v y(yin);

  sse::float4v xw = ((x*a + b)*x + c)*x + d;
	sse::float4v yx = sse::dot4(xw, sse::float4v(xknots+ 0, sse::unaligned),
																	sse::float4v(xknots+ 4, sse::unaligned),
																	sse::float4v(xknots+ 8, sse::unaligned),
																	sse::float4v(xknots+12, sse::unaligned));
	
	sse::float4v yy = sse::dot4(xw, sse::float4v(yknots+ 0, sse::unaligned),
																	sse::float4v(yknots+ 4, sse::unaligned),
																	sse::float4v(yknots+ 8, sse::unaligned),
																	sse::float4v(yknots+12, sse::unaligned));


  sse::float4v yw = ((y*a + b)*y + c)*y + d;

	sse::float4v xyres = sse::dot4(yw, yx, yy);
	xyres.store<0>(xin);
	xyres.store<1>(yin);
}


__forceinline void
bicubic_spline_unrolled_sse4_v2(float& xin, float& yin, const float xknots[16], const float yknots[16]) {
  /*
      -1  3  -3  1
       3 -6   3  0
      -3  0   3  0
      -1  4   1  0
  */

  //((a*x + b)*x + c)*x + d
  static const __m128 a = {-1.0f/6.0f,  3.0f/6.0f,-3.0f/6.0f,  1.0f/6.0f};
  static const __m128 b = { 3.0f/6.0f, -1.0f     ,      0.0f,  4.0f/6.0f};
  static const __m128 c = {-3.0f/6.0f,  3.0f/6.0f, 3.0f/6.0f,  1.0f/6.0f};
  static const __m128 d = { 1.0f/6.0f, 0.0f,       0.0f, 0.0f};

  namespace sse = oyrke::algorithm::sse;
  sse::float4v x(xin);

  sse::float4v xw = ((x*a + b)*x + c)*x + d;
	sse::float4v kx0 = sse::float4v(xknots+ 0, sse::unaligned);
	sse::float4v kx1 = sse::float4v(xknots+ 4, sse::unaligned);
	sse::float4v kx2 = sse::float4v(xknots+ 8, sse::unaligned);
	sse::float4v kx3 = sse::float4v(xknots+12, sse::unaligned);

	sse::float4v yx = sse::dot4(xw, kx0, kx1, kx2, kx3);

	sse::float4v y(yin);
	sse::float4v yw = ((y*a + b)*y + c)*y + d;

	sse::float4v ky0 = sse::float4v(yknots+ 0, sse::unaligned);
	sse::float4v ky1 = sse::float4v(yknots+ 4, sse::unaligned);
	sse::float4v ky2 = sse::float4v(yknots+ 8, sse::unaligned);
	sse::float4v ky3 = sse::float4v(yknots+12, sse::unaligned);
	sse::float4v yy = sse::dot4(xw, ky0, ky1, ky2, ky3);

	sse::float4v xyres = sse::dot4(yw, yx, yy);
	xyres.store<0>(xin);
	xyres.store<1>(yin);
}



template <typename T, typename V>
V
bicubic_spline2(const T& x, const T& y, const V knots[16]) {
  V yw[4];
  yw[0] = cubic_spline(x, &knots[ 0]);
  yw[1] = cubic_spline(x, &knots[ 4]);
  yw[2] = cubic_spline(x, &knots[ 8]);
  yw[3] = cubic_spline(x, &knots[12]);
  V result = cubic_spline(y, yw);
  return result;
}

template <typename T>
coord2d
bicubic_spline2(const T& x, const T& y, const T xknots[16], const T yknots[16]) {
  T yw[4];
  yw[0] = cubic_spline(x, &xknots[ 0]);
  yw[1] = cubic_spline(x, &xknots[ 4]);
  yw[2] = cubic_spline(x, &xknots[ 8]);
  yw[3] = cubic_spline(x, &xknots[12]);
  T xtmp = cubic_spline(y, yw);

  yw[0] = cubic_spline(x, &yknots[ 0]);
  yw[1] = cubic_spline(x, &yknots[ 4]);
  yw[2] = cubic_spline(x, &yknots[ 8]);
  yw[3] = cubic_spline(x, &yknots[12]);

  T ytmp = cubic_spline(y, yw);
  coord2d result(xtmp, ytmp);

  return result;
}


template <typename T, int N>
class interpolator {
  
  static inline void weigthed_sum_recursive(const T *data, T x, T *wsum) {
    size_t M = 1<<N;
    for (size_t i = 0; i < M; ++i) {
      wsum[i] = (1-x)*data[2*i] + x*data[2*i+1];
    }
    //*wsum = (1-x)*data[0] + x*data[1];
    //weigthed_sum_recursive<M-2>(data+2, x, wsum+1);
  }

  //template <>
  //static __forceinline void weigthed_sum_recursive<0>(const T *, T , T *) { }


public:
  static __forceinline void linear_impl(T *data, const T *x) {
    weigthed_sum_recursive(data, x[0], data);
    interpolator<T, N-1>::linear_impl(data, x+1);
  }

  static __forceinline T linear(const T *data, const T *x) {
    T tmp[1<<N];
    weigthed_sum_recursive(data, x[0], tmp);
    interpolator<T,N-1>::linear_impl(tmp, x+1);
    return tmp[0];
  }
};

template <typename T>
class interpolator<T,0> {
public:
  static inline void linear_impl(T *, const T *) { }
  static inline T linear(const T *data, const T *x) { }
};


static const int INTERPOLATE_SAMPLE_SIZE = 10000;
__declspec(align(16)) static float interp_testdata[INTERPOLATE_SAMPLE_SIZE];
__declspec(align(16)) static double dinterp_testdata[INTERPOLATE_SAMPLE_SIZE];

static int 
linear_interp_test() {
  float x = linear_interpolate(interp_testdata[0], interp_testdata[1], interp_testdata[100]);
  return x != 0.0f;
}

static int 
linear_interp2_test() {
  float x = linear_interpolate2(interp_testdata[0], interp_testdata[1], interp_testdata[100]);
  return x != 0.0f;
}


static int 
linear_interp_generic_test() {
  float x = linear_interpolate_generic(interp_testdata, interp_testdata+100);
  return x != 0.0f;
}

static int 
bilinear_interp_plain_test() {
  float x = bilinear_interpolate_plain(interp_testdata[0], interp_testdata[1], interp_testdata[2], interp_testdata[3], 
                                       interp_testdata[100], interp_testdata[101]);
  return x != 0.0f;
} 

static int 
bilinear_interp_opt_test() {
  float x = bilinear_interpolate_opt(interp_testdata[0], interp_testdata[1], interp_testdata[2], interp_testdata[3], 
                                       interp_testdata[100], interp_testdata[101]);
  return x != 0.0f;
}



template <int N>
int
linear_interp_generic_test() {
  float x = linear_interpolate<float, N>(interp_testdata, interp_testdata+100);
  return x != 0.0f;
}

template <int N>
int
linear_interp_genericR_test() {
  float x = interpolator<float, N>::linear(interp_testdata, interp_testdata+100);
  return x != 0.0f;
}


static int
cubic_spline_test() {
  float x = cubic_spline(interp_testdata[0], interp_testdata+100);
  return x != 0.0f;
}

static int
cubic_spline_sse_test() {
  float x = cubic_spline_sse(interp_testdata[0], interp_testdata+128);
  return x != 0.0f;
}

static int
cubic_spline_sse_unaligned_test() {
  float x = cubic_spline_sse_unaligned(interp_testdata[0], interp_testdata+128);
  return x != 0.0f;
}

static int
cubic_spline_sse_overload_test() {
  float x = cubic_spline_sse_v4(interp_testdata[0], interp_testdata+128);
  return x != 0.0f;
}

static int
bicubic_spline_test() {
  float x = bicubic_spline(interp_testdata[0], interp_testdata[1], interp_testdata+100);
  return x != 0.0f;
}


static int
bicubic_spline_sse_test() {
  float x = bicubic_spline_sse(interp_testdata[0], interp_testdata[1], interp_testdata+100);
  return x != 0.0f;
}

static int
bicubic_spline2_test() {
  float x = bicubic_spline2(interp_testdata[0], interp_testdata[1], interp_testdata+100);
  return x != 0.0f;
}


static int
cubic_dspline_test() {
  double x = cubic_spline(dinterp_testdata[0], dinterp_testdata+100);
  return x != 0.0;
}


#if 0
static int
cubic_dspline_sse_test() {
  double x = cubic_spline_sse(dinterp_testdata[0], dinterp_testdata+128);
  return x != 0.0;
}

static int
cubic_dspline_sse_unaligned_test() {
  double x = cubic_spline_sse_unaligned(dinterp_testdata[0], dinterp_testdata+128);
  return x != 0.0;
}

static int
cubic_dspline_sse_overload_test() {
  double x = cubic_spline_sse_v4(dinterp_testdata[0], dinterp_testdata+128);
  return x != 0.0;
}


static int
bicubic_dspline_sse_test() {
  double x = bicubic_spline_sse(dinterp_testdata[0], dinterp_testdata[1], dinterp_testdata+100);
  return x != 0.0;
}

#endif


static int
bicubic_dspline_test() {
  double x = bicubic_spline(dinterp_testdata[0], dinterp_testdata[1], dinterp_testdata+100);
  return x != 0.0;
}


static int
bicubic_dspline2_test() {
  double x = bicubic_spline2(dinterp_testdata[0], dinterp_testdata[1], dinterp_testdata+100);
  return x != 0.0;
}


static int
bicubic_coordspline_test() {
  coord2d *knots = reinterpret_cast<coord2d*>(dinterp_testdata)+100;
  coord2d x = bicubic_spline(dinterp_testdata[0], dinterp_testdata[1], knots);
  return x.x() != 0.0;
}


static int
bicubic_coordspline_opt_test() {
  coord2d x = bicubic_spline(dinterp_testdata[0], dinterp_testdata[1], 
    dinterp_testdata+100, dinterp_testdata+200);
  return x.x() != 0.0;
}


static int
bicubic_coordspline_unrolled_test() {
  coord2d x = bicubic_spline_unrolled(dinterp_testdata[0], dinterp_testdata[1], 
    dinterp_testdata+100, dinterp_testdata+200);
  return x.x() != 0.0;
}

static int
bicubic_coordspline_unrolled_sse_double_test() {
  coord2d x = bicubic_spline_unrolled_sse(dinterp_testdata[0], dinterp_testdata[1], 
    dinterp_testdata+100, dinterp_testdata+200);
  return x.x() != 0.0;
}


static int
bicubic_coordspline_unrolled_sse_v2_double_test() {
  coord2d x = bicubic_spline_unrolled_sse_v2(dinterp_testdata[0], dinterp_testdata[1], 
    dinterp_testdata+100, dinterp_testdata+200);
  return x.x() != 0.0;
}


static int
bicubic_coordspline_unrolled_sse_float_test() {
  float x = interp_testdata[0];
  float y = interp_testdata[1];

  bicubic_spline_unrolled_sse(x, y, interp_testdata+100, interp_testdata+200);
  return x != 0.0f;
}


static int
bicubic_coordspline_unrolled_sse_v2_float_test() {
  float x = interp_testdata[0];
  float y = interp_testdata[1];

  bicubic_spline_unrolled_sse_v2(x, y, interp_testdata+100, interp_testdata+200);
  return x != 0.0f;
}


static int
bicubic_coordspline_unrolled_sse_v3_float_test() {
  float x = interp_testdata[0];
  float y = interp_testdata[1];

  bicubic_spline_unrolled_sse_v3(x, y, interp_testdata+100, interp_testdata+200);
  return x != 0.0f;
}


static int
bicubic_coordspline_unrolled_sse4_float_test() {
  float x = interp_testdata[0];
  float y = interp_testdata[1];

  bicubic_spline_unrolled_sse4(x, y, interp_testdata+100, interp_testdata+200);
  return x != 0.0f;
}


static int
bicubic_coordspline_unrolled_sse4_v2_float_test() {
  float x = interp_testdata[0];
  float y = interp_testdata[1];

  bicubic_spline_unrolled_sse4_v2(x, y, interp_testdata+100, interp_testdata+200);
  return x != 0.0f;
}


static int
bicubic_coordspline_unrolled_xy_test() {
  coord2d *knots = reinterpret_cast<coord2d*>(dinterp_testdata)+100;
  coord2d x = bicubic_spline_unrolled(dinterp_testdata[0], dinterp_testdata[1], knots);
  return x.x() != 0.0;
}



static int
bicubic_coordspline2_test() {
  coord2d *knots = reinterpret_cast<coord2d*>(dinterp_testdata)+100;
  coord2d x = bicubic_spline2(dinterp_testdata[0], dinterp_testdata[1], knots);
  return x.x() != 0.0;
}


static int
bicubic_coordspline2_opt_test() {
  coord2d x = bicubic_spline2(dinterp_testdata[0], dinterp_testdata[1], dinterp_testdata+100, dinterp_testdata+200);
  return x.x() != 0.0;
}


static int 
linear_interp_looptest() {
  float x= 0.0f;
  for (int i = 0; i < 100; ++i) {
    x += linear_interpolate(interp_testdata[i], interp_testdata[i+1], interp_testdata[i+100]);
  }
  return x != 0.0f;
}

static int 
linear_interp2_looptest() {
  float x = 0.0f;
  for (int i = 0; i < 100; ++i) {
    x += linear_interpolate2(interp_testdata[i+0], interp_testdata[i+1], interp_testdata[i+100]);
  }
  return x != 0.0f;
}


static int 
linear_interp_generic_looptest() {
  float x = 0.0f;
  for (int i = 0; i < 100; ++i) {
    x += linear_interpolate_generic(interp_testdata+i, interp_testdata+i+100);
  }
  return x != 0.0f;
}

static int 
bilinear_interp_plain_looptest() {
  float x = 0.0f;
  for (int i = 0; i < 100; ++i) {
    x += bilinear_interpolate_plain(interp_testdata[i+0], interp_testdata[i+1], interp_testdata[i+2], interp_testdata[i+3], 
                                       interp_testdata[i+100], interp_testdata[i+101]);
  }
  return x != 0.0f;
}


static int 
bilinear_interp_opt_looptest() {
  float x = 0.0f;
  for (int i = 0; i < 100; ++i) {
    x += bilinear_interpolate_opt(interp_testdata[i], interp_testdata[i+1], interp_testdata[i+2], interp_testdata[i+3], 
                                       interp_testdata[i+100], interp_testdata[i+101]);
  }
  return x != 0.0f;
}



template <int N>
int
linear_interp_generic_looptest() {
  float x = 0.0f;
  for (int i = 0; i < 100; ++i) {
    x += linear_interpolate<float, N>(interp_testdata+i, interp_testdata+i+100);
  }
  return x != 0.0f;
}



void
test_interpolate_performance() {
  srand((unsigned)time(0));

  double max_test_time = 0.5;
  srand((unsigned)time(0));

  for (int i = 0; i < INTERPOLATE_SAMPLE_SIZE; ++i) {
    interp_testdata[i] = float(rand()) / RAND_MAX;
    dinterp_testdata[i] = interp_testdata[i];
  }

  using namespace oyrke::test;
  namespace sse = oyrke::algorithm::sse;

  run_adaptive_performance_test("Linear interpolate (opt)            : ", max_test_time, linear_interp_test);
  run_adaptive_performance_test("Linear interpolate (opt2)           : ", max_test_time, linear_interp2_test);
  run_adaptive_performance_test("Linear interpolate (generic)        : ", max_test_time, linear_interp_generic_test);
  run_adaptive_performance_test("Bilinear interpolate (plain)        : ", max_test_time, bilinear_interp_plain_test);
  run_adaptive_performance_test("Bilinear interpolate (opt)          : ", max_test_time, bilinear_interp_opt_test);
  //run_adaptive_performance_test("1-inear interpolate (generic)     : ", max_test_time, linear_interp_generic_test<1>);
  //run_adaptive_performance_test("2-linear interpolate (generic)    : ", max_test_time, linear_interp_generic_test<2>);
  //run_adaptive_performance_test("3-linear interpolate (generic)    : ", max_test_time, linear_interp_generic_test<3>);
  //run_adaptive_performance_test("4-linear interpolate (generic)    : ", max_test_time, linear_interp_generic_test<4>);
  //run_adaptive_performance_test("5-linear interpolate (generic)    : ", max_test_time, linear_interp_generic_test<5>);
  //run_adaptive_performance_test("6-linear interpolate (generic)    : ", max_test_time, linear_interp_generic_test<6>);
  //run_adaptive_performance_test("7-linear interpolate (generic)    : ", max_test_time, linear_interp_generic_test<7>);
  //run_adaptive_performance_test("8-linear interpolate (generic)    : ", max_test_time, linear_interp_generic_test<8>);
  //run_adaptive_performance_test("9-linear interpolate (generic)    : ", max_test_time, linear_interp_generic_test<9>);
  //run_adaptive_performance_test("10-linear interpolate (generic)   : ", max_test_time, linear_interp_generic_test<10>);

  //run_adaptive_performance_test("1-inear interpolate (generic R)   : ", max_test_time, linear_interp_genericR_test<1>);
  //run_adaptive_performance_test("2-linear interpolate (generic R)  : ", max_test_time, linear_interp_genericR_test<2>);
  //run_adaptive_performance_test("3-linear interpolate (generic R)  : ", max_test_time, linear_interp_genericR_test<3>);
  //run_adaptive_performance_test("4-linear interpolate (generic R)  : ", max_test_time, linear_interp_genericR_test<4>);
  //run_adaptive_performance_test("5-linear interpolate (generic R)  : ", max_test_time, linear_interp_genericR_test<5>);
  //run_adaptive_performance_test("6-linear interpolate (generic R)  : ", max_test_time, linear_interp_genericR_test<6>);
  //run_adaptive_performance_test("7-linear interpolate (generic R)  : ", max_test_time, linear_interp_genericR_test<7>);
  //run_adaptive_performance_test("8-linear interpolate (generic R)  : ", max_test_time, linear_interp_genericR_test<8>);
  //run_adaptive_performance_test("9-linear interpolate (generic R)  : ", max_test_time, linear_interp_genericR_test<9>);
  //run_adaptive_performance_test("10-linear interpolate (generic R) : ", max_test_time, linear_interp_genericR_test<10>);

  //run_adaptive_performance_test("Linear interpolate x100 (opt)        : ", max_test_time, linear_interp_looptest);
  //run_adaptive_performance_test("Linear interpolate x100 (opt2)       : ", max_test_time, linear_interp2_looptest);
  //run_adaptive_performance_test("Linear interpolate x100 (generic)    : ", max_test_time, linear_interp_generic_looptest);
  //run_adaptive_performance_test("Bilinear interpolate x100 (plain)    : ", max_test_time, bilinear_interp_plain_looptest);
  //run_adaptive_performance_test("Bilinear interpolate x100 (opt)      : ", max_test_time, bilinear_interp_opt_looptest);
  //run_adaptive_performance_test("1-inear interpolate x100 (generic)   : ", max_test_time, linear_interp_generic_looptest<1>);
  //run_adaptive_performance_test("2-linear interpolate x100 (generic)  : ", max_test_time, linear_interp_generic_looptest<2>);
  //run_adaptive_performance_test("3-linear interpolate x100 (generic)  : ", max_test_time, linear_interp_generic_looptest<3>);
  //run_adaptive_performance_test("4-linear interpolate x100 (generic)  : ", max_test_time, linear_interp_generic_looptest<4>);
  //run_adaptive_performance_test("5-linear interpolate x100 (generic)  : ", max_test_time, linear_interp_generic_looptest<5>);
  //run_adaptive_performance_test("6-linear interpolate x100 (generic)  : ", max_test_time, linear_interp_generic_looptest<6>);
  //run_adaptive_performance_test("7-linear interpolate x100 (generic)  : ", max_test_time, linear_interp_generic_looptest<7>);
  //run_adaptive_performance_test("8-linear interpolate x100 (generic)  : ", max_test_time, linear_interp_generic_looptest<8>);
  //run_adaptive_performance_test("9-linear interpolate x100 (generic)  : ", max_test_time, linear_interp_generic_looptest<9>);
  //run_adaptive_performance_test("10-linear interpolate x100 (generic) : ", max_test_time, linear_interp_generic_looptest<10>);
  
  run_adaptive_performance_test("bicubic spline                       : ", max_test_time, bicubic_spline_test);
  run_adaptive_performance_test("bicubic spline (recompute splines)   : ", max_test_time, bicubic_spline2_test);
  run_adaptive_performance_test("bicubic spline SSE                   : ", max_test_time, bicubic_spline_sse_test);

  run_adaptive_performance_test("cubic spline (double)                : ", max_test_time, cubic_dspline_test);
  //run_adaptive_performance_test("cubic spline SSE                     : ", max_test_time, cubic_spline_sse_test);
  //run_adaptive_performance_test("cubic spline SSE unaligned           : ", max_test_time, cubic_spline_sse_unaligned_test);
  //run_adaptive_performance_test("cubic spline SSE math overload       : ", max_test_time, cubic_spline_sse_overload_test);

	run_adaptive_performance_test("cubic spline                         : ", max_test_time, cubic_spline_test);
  run_adaptive_performance_test("cubic spline SSE                     : ", max_test_time, cubic_spline_sse_test);
  run_adaptive_performance_test("cubic spline SSE unaligned           : ", max_test_time, cubic_spline_sse_unaligned_test);
  run_adaptive_performance_test("cubic spline SSE math overload       : ", max_test_time, cubic_spline_sse_overload_test);
  
  run_adaptive_performance_test("bicubic_spline<d,d>(double x, y, knots)                    : ", max_test_time, bicubic_dspline_test);
  run_adaptive_performance_test("bicubic_spline2<d,d>(double x, y, knots)                   : ", max_test_time, bicubic_dspline2_test);
  //run_adaptive_performance_test("bicubic spline SSE                   : ", max_test_time, bicubic_spline_sse_test);
  run_adaptive_performance_test("bicubic_spline<d, 2d>(double x, y, xyknots)                : ", max_test_time, bicubic_coordspline_test);
  run_adaptive_performance_test("bicubic_spline2<d, 2d>(double x, y, xyknots)               : ", max_test_time, bicubic_coordspline2_test);
  run_adaptive_performance_test("bicubic spline<d>(double x, y, xknots, yknots)             : ", max_test_time, bicubic_coordspline_opt_test);
  run_adaptive_performance_test("bicubic_spline2<d>(double x, y, xknots, yknots)            : ", max_test_time, bicubic_coordspline2_opt_test);
  run_adaptive_performance_test("bicubic_spline_unrolled<d>(double x, y, xknots, yknots)    : ", max_test_time, bicubic_coordspline_unrolled_test);
  run_adaptive_performance_test("bicubic_spline_unrolled<d,2d>(double x, y, xyknots)        : ", max_test_time, bicubic_coordspline_unrolled_xy_test);
  run_adaptive_performance_test("bicubic_spline_unrolled_sse(double x, y, xknots, yknots)   : ", max_test_time, bicubic_coordspline_unrolled_sse_double_test);
  run_adaptive_performance_test("bicubic_spline_unrolled_sse_v2(double x, y, xknots, yknots): ", max_test_time, bicubic_coordspline_unrolled_sse_v2_double_test);
  run_adaptive_performance_test("bicubic_spline_unrolled_sse(float x, y, xknots, yknots)    : ", max_test_time, bicubic_coordspline_unrolled_sse_float_test);
  run_adaptive_performance_test("bicubic_spline_unrolled_sse_v2(float x, y, xknots, yknots) : ", max_test_time, bicubic_coordspline_unrolled_sse_v2_float_test);
  run_adaptive_performance_test("bicubic_spline_unrolled_sse_v3(float x, y, xknots, yknots) : ", max_test_time, bicubic_coordspline_unrolled_sse_v3_float_test);
  run_adaptive_performance_test("bicubic_spline_unrolled_sse4(float x, y, xknots, yknots)   : ", max_test_time, bicubic_coordspline_unrolled_sse4_float_test);
  run_adaptive_performance_test("bicubic_spline_unrolled_sse4_v2(float x, y, xknots, yknots): ", max_test_time, bicubic_coordspline_unrolled_sse4_v2_float_test);
}

// local variable is initialized but not referenced
#pragma warning(push)
#pragma warning(disable: 4189)

template <typename U, typename V> 
struct cast {
  V operator ()(const U& x) const { return V(x); }
};

void
test_interpolate() {
  float y = linear_interpolate(1.0f, 5.0f, 0.25f); // expect 2

  float y2_plain = bilinear_interpolate_plain(1.0f, 2.0f, 4.0f, 8.0f, 0.25f, 0.75f); // 4.0625
  float y2_opt   = bilinear_interpolate_opt(1.0f, 2.0f, 4.0f, 8.0f, 0.25f, 0.75f); // 4.0625

  float data[2] = { 1.0f, 5.0f };
  float x = 0.25f;

  __declspec(align(16)) float data2[4] = { 1.0f, 2.0f, 4.0f, 8.0f };
  float x2[2]    = { 0.25f, 0.75f };

  float yy = linear_interpolate<float, 1>(data, &x);
  float yy2 = linear_interpolate<float, 2>(data2, x2);

  yy = interpolator<float, 1>::linear(data, &x);
  yy2= interpolator<float, 2>::linear(data2, x2);

  float cx = cubic_spline(0.5f, data2);
  float cx_sse = cubic_spline_sse(0.5f, data2);

  double data16[] = {1.0,  2.0,  3.0,  4.0,
                     5.0,  6.0,  7.0,  8.0,
                     9.0, 10.0, 11.0, 12.0,
                    13.0, 14.0, 15.0, 16.0};
  double bcx = bicubic_spline(0.2, 0.6, data16);

  coord2d coord16[] = {coord2d( 1.0, 51.0), coord2d( 2.0, 52.0), coord2d( 3.0, 53.0), coord2d( 4.0, 54.0),
                       coord2d( 5.0, 55.0), coord2d( 6.0, 56.0), coord2d( 7.0, 57.0), coord2d( 8.0, 58.0),
                       coord2d( 9.0, 59.0), coord2d(10.0, 60.0), coord2d(11.0, 61.0), coord2d(12.0, 62.0),
                       coord2d(13.0, 63.0), coord2d(14.0, 64.0), coord2d(15.0, 65.0), coord2d(16.0, 66.0)};

  double data16y[16];
  std::transform(data16, data16+16, data16y, std::bind1st(std::plus<double>(), 50.0));
  coord2d ccx = bicubic_spline(0.2, 0.6, coord16);
  coord2d ccx2= bicubic_spline(0.2, 0.6, data16, data16y);
  coord2d ccx3= bicubic_spline_unrolled(0.2, 0.6, coord16);
  coord2d ccx4= bicubic_spline_unrolled(0.2, 0.6, data16, data16y);
  coord2d ccx5= bicubic_spline_unrolled_sse(0.2, 0.6, data16, data16y);

  float fdata16x[16];
  float fdata16y[16];


  std::transform(data16, data16+16, fdata16x, cast<double, float>());
  std::transform(data16y, data16y+16, fdata16y, cast<double, float>());

  float xc=0.2f;
  float yc=0.6f;
  bicubic_spline_unrolled_sse(xc, yc, fdata16x, fdata16y);

  __declspec(align(16)) float data16f[] = {1.0f,  2.0,  3.0,  4.0,
                     5.0,  6.0,  7.0,  8.0,
                     9.0, 10.0, 11.0, 12.0,
                    13.0, 14.0, 15.0, 16.0};
  float bcx_sse = bicubic_spline_sse(0.2, 0.6, data16f);
}
