#include <Utility/sse.h>
#include <Utility/math.h>
#include <cmath>


#ifdef NDEBUG
volatile int global_counter = 0;

void assert(bool ok) {
	if (!ok) {
		global_counter++;
	}
}
#else
#		include <assert.h>
#endif

void
test_sse_float() {
  using namespace oyrke::algorithm::sse;

  float x = 4.0f;
  float4v same(x);
  assert(same[0] == x);
  assert(same[3] == x);

  float data1[4] = {1.0f, 8.0f, 64.0f, 512.0f};
  float4v pow8(data1, unaligned);
  assert(pow8[0] == data1[0]);
  assert(pow8[1] == data1[1]);
  assert(pow8[2] == data1[2]);
  assert(pow8[3] == data1[3]);

	float4v args4(1.0f, 8.0f, 64.0f, 512.0f);
	assert(args4.equals(pow8));
	float4v sum=sum_row(pow8);
  assert(sum[0] == 585.0f); // 1+8+64+512

  float data2[] = {10.0f, 20.0f, 30.0f, 40.0f};
  float4v sum2 = sum_row(pow8, float4v(data2, unaligned));
  assert(sum2[0] == 585.0f);
  assert(sum2[1] == 100.0f);

  float data3[] = {100.0f, 300.0f, 600.0f, 1100.0f};
  float data4[] = {21.0f, 22.0f, 23.0f, 24.0f};

  float4v sum4 = sum_row(pow8, float4v(data2, unaligned), float4v(data3, unaligned), float4v(data4, unaligned));
  assert(sum4[0] == 585.0f);
  assert(sum4[1] == 100.0f);
  assert(sum4[2] == 2100.0f);
  assert(sum4[3] == 90.0f);

	float4v abstest(1.0, -3.14, -4.987654, 6.28);
	float4v abs2 = math::abs(abstest);
  assert(abs2[0] == 1.0f);
  assert(abs2[1] == 3.14f);
  assert(abs2[2] == 4.987654f);
  assert(abs2[3] == 6.28f);

	abs2 = -abstest;
  assert(abs2[0] == -1.0f);

	float4v rcp = math::reciprocal(abstest);
	float4v rcp_nr = math::reciprocal_nr(abstest);
	float4v rcp_cc = float4v(1.0f)/abstest;

	float4v rsqrt    = math::sqrt_reciprocal(abstest);
	float4v rsqrt_nr = math::sqrt_reciprocal_nr(abstest);
	float4v rsqrt_cc = float4v(1.0f)/math::sqrt(abstest);


	float4v nans=float4v::zeros()/float4v::zeros();
	bool ok1 = nans.at<2>() == nans.at<2>();
	bool ok2 = _mm_comieq_ss(nans.value(), nans.value()) != 0;
	//assert(!ok1);
	assert(ok2);
}


void
test_sse_double() {
  using namespace oyrke::algorithm::sse;

  double x = 4.0;
  double2v same(x);
  assert(same[0] == x);
  assert(same[1] == x);

  double data1[2] = {1.0, 8.0};
  double2v pow8(data1, unaligned);
  assert(pow8[0] == data1[0]);
  assert(pow8[1] == data1[1]);

  double2v sum=sum_row(pow8);
  assert(sum[0] == 9.0); // 1+8

  double data2[] = {10.0, 20.0};
  double2v sum2 = sum_row(pow8, double2v(data2, unaligned));
  assert(sum2[0] == 9.0f);
  assert(sum2[1] == 30.0f);
}

using namespace oyrke::algorithm::sse;

void
test_sse_mask() {
	mask4v mask = mask4v::mask<2>();
	mask = mask4v::mask<0>();
	mask = mask4v::zeros();
	mask = mask4v::ones();
}

bool
almost_equal(float x, float y, float eps = 1e-6) {
	float diff = abs(x-y);
	float sum  = abs(x) + abs(y);
	bool  eq   = diff < eps || diff/sum < eps;
	return eq;
}

mask4v
almost_equal(const float4v& x, const float4v& y, float eps = 1e-6) {
	float4v diff = math::abs(x-y);
	float4v sum  = math::abs(x) + math::abs(y);
	mask4v  eq   = (diff < eps) || (diff/sum < eps);
	return eq;
}

void
test_sse_func(const float4v& x, float4v (*func4v)(const float4v&), float (*funcf)(float)) {
	float4v ys = func4v(x);
	float4v ys_correct(funcf(x[0]), funcf(x[1]), funcf(x[2]), funcf(x[3]));
	assert(almost_equal(ys, ys_correct).all_true());
}

void
test_sse_func_loop(float4v (*func4v)(const float4v&), float (*funcf)(float)) {
	for (int i = -10; i<10; ++i) {
		float4v x = float4v(0.1, 0.3, 0.7, 1.1) + float(i);
		test_sse_func(x, func4v, funcf);
	}
}

void
test_sse_math() {

	test_sse_func_loop(math::sin, sin);
	test_sse_func_loop(math::cos, cos);

	{
		float4v x = float4v(0.1, 0.3, -3.2, -1.1);
		float4v sin_correct = math::sin(x);
		float4v cos_correct = math::cos(x);
		float4v ysin, ycos;
		math::sincos(x, ysin, ycos);
		assert(almost_equal(ysin, sin_correct).all_true());
		assert(almost_equal(ycos, cos_correct).all_true());
	}
	//float yc = math::cos(x);
	//assert(almost_equal(cos(x.at<0>()), yc.at<0>()));
	//assert(almost_equal(cos(x.at<1>()), yc.at<1>()));
	//assert(almost_equal(cos(x.at<2>()), yc.at<2>()));
	//assert(almost_equal(cos(x.at<3>()), yc.at<3>()));
	{
		float4v x(0.001, 0.125, 8, 1000);
		float4v y1 = math::cbrt(x);
		float4v y2 = math::cbrt_halleys(x);

		x = float4v(1e-12, 1e-24, 1e-36, 1e-45);
		y1 = math::cbrt(x);
		y2 = math::cbrt_halleys(x);
	}

	{
		float4v x(2.0, 3.0, 4.0, 5.0);
		float4v z = math::pow(x, x);
		z = math::pow(x, 7.1-x);
		float4v z2 = math::pow_v2(x, 7.1-x);
		float z0c = pow(x.at<0>(), 7.1f-x.at<0>());
		float z1c = pow(x.at<1>(), 7.1f-x.at<1>());
		float z2c = pow(x.at<2>(), 7.1f-x.at<2>());
		float z3c = pow(x.at<3>(), 7.1f-x.at<3>());
		float4v diff = z2-z;
	}
}

#include <limits>

void
test_sse_ieee_specials() {
	typedef std::numeric_limits<float> limits_t;
	float4v specials(limits_t::quiet_NaN(), limits_t::signaling_NaN(), limits_t::infinity(), -limits_t::infinity());

	int4v expnt;
	float4v y;
	y = math::frexp(specials, expnt);
	assert((y.reinterpret_bits() == specials.reinterpret_bits()).all_true());

	float4v values(0.0f, limits_t::denorm_min(), limits_t::min(), limits_t::max());
	y = math::frexp(values, expnt);

	int e0=0, e1=0, e2=0, e3=0;
	float ys0=frexp(values[0], &e0);
	float ys1=frexp(values[1], &e1);
	float ys2=frexp(values[2], &e2);
	float ys3=frexp(values[3], &e3);
	assert(y[0] == ys0);
	assert(y[1] == ys1);
	assert(y[2] == ys2);
	assert(y[3] == ys3);
	assert(expnt[0] == e0);
	assert(expnt[1] == e1);
	assert(expnt[2] == e2);
	assert(expnt[3] == e3);
	float4v negvalues = -values;
	int4v expnt2;
	float4v y2 = math::frexp(negvalues, expnt2);
	ys0=frexp(negvalues[0], &e0);
	ys1=frexp(negvalues[1], &e1);
	ys2=frexp(negvalues[2], &e2);
	ys3=frexp(negvalues[3], &e3);

	assert(y2[0] == ys0);
	assert(y2[1] == ys1);
	assert(y2[2] == ys2);
	assert(y2[3] == ys3);
	assert(expnt[0] == e0);
	assert(expnt[1] == e1);
	assert(expnt[2] == e2);
	assert(expnt[3] == e3);
	assert((y2 == -y).all_true());
	assert((expnt2 == expnt).all_true());

	float4v normal(1.0, 2.0, 3.0, 10.0);
	int4v exp_norm;
	float4v mant_norm = math::frexp(normal, exp_norm);
	float4v rr_normal = math::ldexp(mant_norm, exp_norm);
	assert((rr_normal == normal).all_true());

	float4v rr_values    = math::ldexp(y, expnt);
	float4v rr_negvalues = math::ldexp(y2, expnt2);

	assert((rr_values == values).all_true());
	assert((rr_negvalues == negvalues).all_true());
}

#include <immintrin.h>

void
test_avx() {
  float data[] = { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8};
  __m256 x = _mm256_loadu_ps(data);
  __m256 twox = _mm256_add_ps(x, x);
}


void
test_sse_sort() {
    float4v x(4.0f, 3.0f, 2.0f, 1.0f);
    float4v y = math::sort_4_asc(x);

    x = float4v(4.0f, 2.0f, 1.0f, 3.0f);
    y = math::sort_4_asc(x);
}


void
test_sse() {
    test_sse_sort();
  test_sse_float();
  test_sse_double();
    test_sse_mask();
	test_sse_math(); 
	test_sse_ieee_specials();
  test_avx();
}
