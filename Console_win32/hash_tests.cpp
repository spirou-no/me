#include <utility/performance_testing.h>
#include <utility/sprintf.h>
#include <utility/sse.h>
#include <utility/math.h>
#include <utility/tiny_algo.h>

#include <math.h>
#include <cmath>
#include <time.h>
#include <functional>

/*
Empty           :       974.095 ms   267810742        3.64 ns
Hash int[2]     :      1067.609 ms   278854603        3.83 ns
Hash int[3]     :       965.901 ms   258532475        3.74 ns
Hash int[6]     :      1048.198 ms   273790065        3.83 ns
Hash int[7]     :      1030.614 ms    50872257       20.26 ns
Hash int[8]     :       964.708 ms    42400467       22.75 ns
Hash int[9]     :       998.510 ms    43488791       22.96 ns
Hash int[10]    :      1016.837 ms    46341351       21.94 ns
*/

/*
Empty                :      1004.435 ms   262438643        3.83 ns
Hash int[2]          :      1013.877 ms   259202308        3.91 ns
Hash int[3]          :       900.895 ms   248350651        3.63 ns
Hash int[6]          :       978.496 ms   269504274        3.63 ns
Hash int[7]          :      1012.726 ms    51536723       19.65 ns
Hash int[8]          :       870.739 ms    38432619       22.66 ns
Hash int[9]          :      1011.629 ms    45968157       22.01 ns
Hash int[10]         :      1003.838 ms    45758024       21.94 ns
Hash int[2]  plain   :       918.889 ms    93205424        9.86 ns
Hash int[3]  plain   :      1058.051 ms    78095642       13.55 ns
Hash int[6]  plain   :      1013.640 ms    40192691       25.22 ns
Hash int[7]  plain   :      1052.210 ms    34466457       30.53 ns
Hash int[8]  plain   :       977.595 ms    29040770       33.66 ns
Hash int[9]  plain   :      1064.640 ms    26585017       40.05 ns
Hash int[10] plain   :      1026.654 ms    22917418       44.80 ns


=== Force inline ===
Empty                :       992.301 ms   271638407        3.65 ns
Hash int[2]          :       764.815 ms   205739895        3.72 ns
Hash int[3]          :      1019.098 ms   257241958        3.96 ns
Hash int[6]          :       993.982 ms   273940926        3.63 ns
Hash int[7]          :      1002.705 ms   277657375        3.61 ns
Hash int[8]          :      1010.092 ms   268705963        3.76 ns
Hash int[9]          :       884.328 ms   244711709        3.61 ns
Hash int[10]         :      1009.221 ms   279080057        3.62 ns
Hash int[2]  plain   :       989.166 ms   100000000        9.89 ns
Hash int[3]  plain   :       881.133 ms    67392990       13.07 ns
Hash int[6]  plain   :      1056.846 ms    40197340       26.29 ns
Hash int[7]  plain   :      1010.131 ms    34403782       29.36 ns
Hash int[8]  plain   :       873.248 ms    26203961       33.33 ns
Hash int[9]  plain   :       990.300 ms    25943206       38.17 ns
Hash int[10] plain   :       938.962 ms    21068118       44.57 ns
Hash int[2]  dyn     :       992.262 ms   235609691        4.21 ns
Hash int[3]  dyn     :      1006.875 ms   209291807        4.81 ns
Hash int[6]  dyn     :       932.191 ms   100000000        9.32 ns
Hash int[7]  dyn     :       962.592 ms   100000000        9.63 ns
Hash int[8]  dyn     :       943.468 ms    87742763       10.75 ns
Hash int[9]  dyn     :       965.496 ms    80309995       12.02 ns
Hash int[10] dyn     :      1044.525 ms    75960610       13.75 ns
Hash int[3] /100     :       905.772 ms     2043185      443.31 ns
Hash int[6] /100     :       929.870 ms     1000000      929.87 ns
Hash int[7] /100     :      1040.003 ms      930167     1118.08 ns
Hash int[10]/100     :       969.603 ms      616927     1571.66 ns


Empty                :      7884.086 ms  2147483647        3.67 ns
Hash int[2]          :      7955.177 ms  2147483647        3.70 ns
Hash int[3]          :      7917.969 ms  2147483647        3.69 ns
Hash int[6]          :      7860.363 ms  2147483647        3.66 ns
Hash int[7]          :      7970.183 ms  2147483647        3.71 ns
Hash int[8]          :      7865.193 ms  2147483647        3.66 ns
Hash int[9]          :      7899.646 ms  2147483647        3.68 ns
Hash int[10]         :      7914.187 ms  2147483647        3.69 ns
Hash int[2]  plain   :      9652.351 ms  1000000000        9.65 ns
Hash int[3]  plain   :     10076.888 ms   757245037       13.31 ns
Hash int[6]  plain   :      9805.752 ms   384101766       25.53 ns
Hash int[7]  plain   :     10189.697 ms   337718062       30.17 ns
Hash int[8]  plain   :      9906.050 ms   290439700       34.11 ns
Hash int[9]  plain   :     10070.292 ms   256129005       39.32 ns
Hash int[10] plain   :      9906.798 ms   228365008       43.38 ns
Hash int[2]  dyn     :      9236.211 ms  2147483647        4.30 ns
Hash int[3]  dyn     :     10032.301 ms  2043898703        4.91 ns
Hash int[6]  dyn     :      8534.163 ms  1000000000        8.53 ns
Hash int[7]  dyn     :      9856.913 ms   989385419        9.96 ns
Hash int[8]  dyn     :     10202.101 ms   923669881       11.05 ns
Hash int[9]  dyn     :     10270.782 ms   833987638       12.32 ns
Hash int[10] dyn     :     10240.531 ms   759350371       13.49 ns
Hash int[3] /100     :     10098.830 ms    22485700      449.12 ns
Hash int[6] /100     :      9350.508 ms     9805997      953.55 ns
Hash int[7] /100     :      9755.947 ms     8873675     1099.43 ns
Hash int[10]/100     :     10242.187 ms     6369148     1608.09 ns


=== Plain inline ===
Empty                :      2147.491 ms   542062507        3.96 ns
Hash int[2]          :      2095.962 ms   556751981        3.76 ns
Hash int[3]          :      2196.997 ms   527468697        4.17 ns
Hash int[6]          :      1620.843 ms   421060772        3.85 ns
Hash int[7]          :      2012.620 ms   100000000       20.13 ns
Hash int[8]          :      1565.669 ms    67539604       23.18 ns
Hash int[9]          :      2188.543 ms    86427263       25.32 ns
Hash int[10]         :      1820.355 ms    80602319       22.58 ns
Hash int[2]  plain   :      1976.960 ms   200315340        9.87 ns
Hash int[3]  plain   :      1349.799 ms   100000000       13.50 ns
Hash int[6]  plain   :      2246.202 ms    79288900       28.33 ns
Hash int[7]  plain   :      2055.260 ms    68804853       29.87 ns
Hash int[8]  plain   :      2073.593 ms    60111270       34.50 ns
Hash int[9]  plain   :      2053.893 ms    52696861       38.98 ns
Hash int[10] plain   :      2070.966 ms    47266169       43.81 ns
Hash int[2]  dyn     :      1831.828 ms   435403669        4.21 ns
Hash int[3]  dyn     :      2107.596 ms   417977781        5.04 ns
Hash int[6]  dyn     :      1993.318 ms   225563926        8.84 ns
Hash int[7]  dyn     :      1705.730 ms   100000000       17.06 ns
Hash int[8]  dyn     :      1795.907 ms   100000000       17.96 ns
Hash int[9]  dyn     :      1986.656 ms    92246566       21.54 ns
Hash int[10] dyn     :      2063.584 ms    89283496       23.11 ns
Hash int[3] /100     :      1847.885 ms     4174077      442.71 ns
Hash int[6] /100     :      1943.786 ms     2042701      951.58 ns
Hash int[7] /100     :      2067.442 ms     1801644     1147.53 ns
Hash int[10]/100     :      1566.677 ms     1000000     1566.68 ns


=== Hasher: seed = a*(*src^invert + seed)
Empty                :      2049.877 ms   553039865        3.71 ns
Hash int[2]          :      2005.473 ms   559910183        3.58 ns
Hash int[3]          :      2068.587 ms   560422562        3.69 ns
Hash int[6]          :      1669.615 ms   100000000       16.70 ns
Hash int[7]          :      1726.372 ms   100000000       17.26 ns
Hash int[8]          :      1849.446 ms   100000000       18.49 ns
Hash int[9]          :      1989.488 ms    78822380       25.24 ns
Hash int[10]         :      2057.337 ms    77415988       26.58 ns
Hash int[2]  plain   :      2050.578 ms   209497102        9.79 ns
Hash int[3]  plain   :      2046.720 ms   153074581       13.37 ns
Hash int[6]  plain   :      2002.152 ms    85107604       23.52 ns
Hash int[7]  plain   :      2063.438 ms    74683935       27.63 ns
Hash int[8]  plain   :      2003.426 ms    64630458       31.00 ns
Hash int[9]  plain   :      2061.569 ms    58205734       35.42 ns
Hash int[10] plain   :      2001.526 ms    46750401       42.81 ns
Hash int[2]  dyn     :      2056.353 ms   418935604        4.91 ns
Hash int[3]  dyn     :      1999.773 ms   373407600        5.36 ns
Hash int[6]  dyn     :      1384.165 ms   100000000       13.84 ns
Hash int[7]  dyn     :      1466.471 ms   100000000       14.66 ns
Hash int[8]  dyn     :      1518.148 ms   100000000       15.18 ns
Hash int[9]  dyn     :      2014.062 ms   100000000       20.14 ns
Hash int[10] dyn     :      1999.031 ms    98373470       20.32 ns
Hash int[3] /100     :      2059.016 ms     5026151      409.66 ns
Hash int[6] /100     :      1998.817 ms     2572207      777.08 ns
Hash int[7] /100     :      2054.098 ms     2197305      934.83 ns
Hash int[10]/100     :      2000.474 ms     1615611     1238.22 ns

=== Hasher: seed = a*(*src^invert + seed) Force inline
Empty                :      2047.223 ms   552581834        3.70 ns
Hash int[2]          :      2007.311 ms   559878656        3.59 ns
Hash int[3]          :      2062.697 ms   558736875        3.69 ns
Hash int[6]          :      1997.105 ms   559410981        3.57 ns
Hash int[7]          :      2060.857 ms   560268618        3.68 ns
Hash int[8]          :      2000.811 ms   559721073        3.57 ns
Hash int[9]          :      2059.414 ms   559582823        3.68 ns
Hash int[10]         :      1991.084 ms   557219556        3.57 ns
Hash int[2]  plain   :      2060.501 ms   209448559        9.84 ns
Hash int[3]  plain   :      1999.176 ms   154041146       12.98 ns
Hash int[6]  plain   :      2060.013 ms    85164605       24.19 ns
Hash int[7]  plain   :      2003.236 ms    74617147       26.85 ns
Hash int[8]  plain   :      2079.746 ms    64693358       32.15 ns
Hash int[9]  plain   :      2000.958 ms    58186906       34.39 ns
Hash int[10] plain   :      2048.624 ms    46437571       44.12 ns
Hash int[2]  dyn     :      2000.667 ms   420433924        4.76 ns
Hash int[3]  dyn     :      2060.874 ms   373688856        5.51 ns
Hash int[6]  dyn     :      1999.765 ms   224116638        8.92 ns
Hash int[7]  dyn     :      1888.287 ms   186649914       10.12 ns
Hash int[8]  dyn     :      2057.440 ms   176540505       11.65 ns
Hash int[9]  dyn     :      1998.646 ms   151275396       13.21 ns
Hash int[10] dyn     :      1464.839 ms   100000000       14.65 ns
Hash int[3] /100     :      2052.358 ms     5032826      407.79 ns
Hash int[6] /100     :      2042.069 ms     2548666      801.23 ns
Hash int[7] /100     :      1999.990 ms     2201624      908.42 ns
Hash int[10]/100     :      1905.829 ms     1540044     1237.52 ns
*/

static const int fixed_pos[] = {102, 14, 23, 20, 95, 55, 61, 70, 3, 83};

template <size_t N>
struct hash_const_index {
  static int run_it() {
    return (int) oyrke::algorithm::tiny_algo<N>::hash(fixed_pos);
  }
};


static const int HASH_SAMPLE_SIZE = 10000;
static int hash_sample_data[HASH_SAMPLE_SIZE];
__declspec(align(16)) static float real_sample_data[HASH_SAMPLE_SIZE];

template <size_t N>
struct hash_dynamic_index {
  static int run_it() {
    return (int) oyrke::algorithm::tiny_algo<N>::hash(hash_sample_data);
  }
};

const int HASHER_LOOPS = 100;

template <size_t N>
struct hash_dynloop_index {
  static int run_it() {
    size_t hashcode = 0;
    for (int i = 0; i < HASHER_LOOPS; ++i) {
      int *pos = hash_sample_data + i;
      hashcode += oyrke::algorithm::tiny_algo<N>::hash(pos);
    }
    return (int) hashcode;
  }
};


struct magic_numbers {
  static const size_t seed        =        123u;
  static const size_t multiplier  = 1592763461u;
  static const size_t additive    =         13u;
  static const size_t mult_2      =  910047527u;
  static const size_t invert      = 0x945da3b2u;
};

template <typename InIter>
static __forceinline /*inline*/ size_t hash(InIter first, InIter last) {
//  BOOST_STATIC_ASSERT(std::numeric_limits<T>::is_specialized);
//  BOOST_STATIC_ASSERT(std::numeric_limits<T>::is_integer);

  size_t seed = magic_numbers::seed;
  for (; first != last; ++first) {
    //seed = magic_numbers::multiplier * (size_t(*first) + seed) + magic_numbers::additive;
    //seed = magic_numbers::multiplier * size_t(*first) + magic_numbers::mult_2 * (*first ^ seed);
    seed = magic_numbers::multiplier * ((size_t(*first) ^ magic_numbers::invert) + seed);
  }
  return seed;
}


template <size_t N>
struct hash_fixplain_index {
  static int run_it() {
    return (int) hash(fixed_pos, fixed_pos+N);
  }
};


static int
empty_test() {
  return 1;
}

template <typename BinaryOp, int N>
int binary_function_test() {
  int sum = 0;
  for (int i = 0; i < N; ++i) {
    sum += static_cast<int>(BinaryOp()(hash_sample_data[i],hash_sample_data[i+1]));
  }
  return sum;
}

template <typename UnaryOp, int N>
int unary_function_test() {
  int sum = 0;
  for (int i = 0; i < N; ++i) {
    sum += UnaryOp()(hash_sample_data[i]);
  }
  return sum;
}


template <typename T, int N>
struct divide_const { int operator()(int x) const { return x/N; } };
template <typename T, int N>
struct right_shift_const { int operator()(int x) const { return x >> N; } };

template <typename T>
struct shift_left : std::binary_function<T,T,T> {
  T operator()(const T& lhs, const T& rhs) const { return lhs << rhs; }
};

template <typename T>
struct shift_right : std::binary_function<T,T,T> {
  T operator()(const T& lhs, const T& rhs) const { return lhs >> rhs; }
};

struct index_manip_base {
  virtual bool compute(double& x, double& y) const = 0;
};

struct index_manip :  public index_manip_base {
  double xmin_, ymin_;
  double inv_xinc_, inv_yinc_;
  double dnx_, dny_;
  int nx_, ny_;

  int buffer_index(double x, double y) const {
    double norm_x = (x-xmin_)*inv_xinc_;
    double norm_y = (y-ymin_)*inv_yinc_;

    bool ok = (norm_x >= 0.0) && (norm_x <= dnx_)
           && (norm_y >= 0.0) && (norm_y <= dny_);
    int base_ix = -1;
    if (ok) {
      int norm_i = int(norm_x);
      int norm_j = int(norm_y);
      base_ix = norm_i + norm_j*nx_;
    }
    return base_ix;
  }

  virtual void init() {
    dnx_ = nx_-1;
    dny_ = ny_-1;
  }
  bool compute(double& x, double& y) const;
};

struct bilinear_index_manip :  public index_manip {
  // all inherited members ok
  double *knots;
  bool compute(double& x, double& y) const;
};

struct bilinear_index_manip_sse :  public index_manip {
  // all inherited members ok
  oyrke::algorithm::sse::double2v min_, inv_inc_, lower_limit_, upper_limit_, *knots_sse;
  void init() {
    namespace sse = oyrke::algorithm::sse;
    index_manip::init();
    min_ = sse::double2v(xmin_, ymin_);
    upper_limit_ = sse::double2v(dnx_, dny_);
    inv_inc_ = sse::double2v(inv_xinc_, inv_yinc_);
  }
  
  bool compute(double& x, double& y) const;
};

struct bilinear_index_manip_ssenan :  public bilinear_index_manip_sse {
  // all inherited members ok
  bool compute(double& x, double& y) const;
  void compute_pos(double pos[2]) const;
};

#if 0
__forceinline int double2int(double flt)
{
  int intgr;

  _asm fld flt
  _asm fistp intgr         
         
  return intgr ;
}
#endif

bool
index_manip::compute(double& x, double& y) const {
  double norm_x = (x-xmin_)*inv_xinc_;
  double norm_y = (y-ymin_)*inv_yinc_;

  bool ok = (norm_x >= 0.0) && (norm_x <= dnx_)
         && (norm_y >= 0.0) && (norm_y <= dny_);
  if (ok) {
    int norm_i = int(norm_x);
    int norm_j = int(norm_y);
    //int norm_i = double2int(norm_x);
    //int norm_j = double2int(norm_y);
    x=norm_x - norm_i;
    y=norm_y - norm_j;
  }
  return ok;
}



//template <typename T>
//__forceinline
//T bilinear_interpolate(T x, T y, T f00, T f10, T f01, T f11) {
//  return (1-x)*(f00*(1-y) + f01*y) + x*(f10*(1-y) + f11*y);
//}

template <typename T>
__forceinline
T bilinear_interpolate(T x, T y, T f00, T f10, T f01, T f11) {
  T w0 = f00 + (f01 - f00)*y;
  T w1 = f10 + (f11 - f10)*y;
  return w0 + (w1 - w0)*x;
  // return (1-x) * w0 + x * w1;
  // return (1-x)*(f00*(1-y) + f01*y) + x*(f10*(1-y) + f11*y);
}



__forceinline
oyrke::algorithm::sse::double2v bilinear_interpolate_sse(
  const oyrke::algorithm::sse::double2v& xy, 
  const oyrke::algorithm::sse::double2v& f00, const oyrke::algorithm::sse::double2v& f01, 
  const oyrke::algorithm::sse::double2v& f10, const oyrke::algorithm::sse::double2v& f11
) {
  namespace sse = oyrke::algorithm::sse;

  sse::double2v X=sse::double2v::compose_00(xy, xy);
  sse::double2v Y=sse::double2v::compose_11(xy, xy);
  return bilinear_interpolate(X, Y, f00, f01, f10, f11);
}


bool
bilinear_index_manip::compute(double& x, double& y) const {
  double norm_x = (x-xmin_)*inv_xinc_;
  double norm_y = (y-ymin_)*inv_yinc_;

  bool ok = (norm_x >= 0.0) && (norm_x <= dnx_)
         && (norm_y >= 0.0) && (norm_y <= dny_);
  if (ok) {
    int norm_i = int(norm_x);
    int norm_j = int(norm_y);
    //int norm_i = double2int(norm_x);
    //int norm_j = double2int(norm_y);
    x=norm_x - norm_i;
    y=norm_y - norm_j;
    int base_ix = norm_i + norm_j*nx_;

    double x1 = bilinear_interpolate(x, y, knots[base_ix],     knots[base_ix+1],     knots[base_ix+nx_],       knots[base_ix+nx_ + 1]);
    double y1 = bilinear_interpolate(x, y, knots[base_ix+400], knots[base_ix+1+400], knots[base_ix+nx_ + 400], knots[base_ix+nx_ + 1 + 400]);
    //double x1 = bilinear_interpolate(knots[norm_i],   knots[norm_i+2],   knots[norm_i+norm_j],   knots[norm_i + norm_j+2], x, y);
    //double y1 = bilinear_interpolate(knots[norm_i+1], knots[norm_i+2+1], knots[norm_i+norm_j+1], knots[norm_i + norm_j+2+1], x, y);
    x = x1;
    y = y1;
  }
  return ok;
}

bool
bilinear_index_manip_sse::compute(double& x, double& y) const {
  namespace sse = oyrke::algorithm::sse;
  sse::double2v in(x, y);
  sse::double2v normalized = (in - min_)*inv_inc_;

  bool ok = (normalized >= sse::double2v::zeros()).and20(normalized < upper_limit_).all_true();
  if (ok) {
    sse::int4v index = normalized.trunc_int();
    sse::double2v fraction= normalized - sse::double2v(index);
    int i = index[0];
    int j = index[1];
    int base_ix = i + j*nx_;
    sse::double2v result = bilinear_interpolate_sse(fraction, knots_sse[base_ix], knots_sse[base_ix+1], knots_sse[base_ix+nx_], knots_sse[base_ix+nx_+1]);
    x = result[0];
    y = result[1];
  }
  return ok;
}


bool
bilinear_index_manip_ssenan::compute(double& x, double& y) const {
  namespace sse = oyrke::algorithm::sse;

  sse::double2v in(x, y);
  sse::double2v normalized = (in - min_)*inv_inc_;
  normalized = normalized.max(sse::double2v(1.0)).min(upper_limit_);  // clip index
  
  sse::int4v index = normalized.trunc_int();
  //sse::double2v fraction= normalized - sse::double2v(index);
  sse::double2v fraction= normalized - normalized.floor();
  //int i = index[0];
  //int j = index[1];
  int i = index.at<0>();
  int j = index.at<1>();
  int base_ix = i + j*nx_;
  sse::double2v result = bilinear_interpolate_sse(fraction, knots_sse[base_ix], knots_sse[base_ix+1], knots_sse[base_ix+nx_], knots_sse[base_ix+nx_+1]);

  result.store<0>(x);
	result.store<1>(y);

	bool ok = _mm_comieq_sd(result.value(), result.value()) != 0;
	return ok;
}

void
bilinear_index_manip_ssenan::compute_pos(double pos[2]) const {
  namespace sse = oyrke::algorithm::sse;

  sse::double2v in(pos, sse::unaligned);
  sse::double2v normalized = (in - min_)*inv_inc_;
  normalized = normalized.max(sse::double2v::zeros()).min(upper_limit_);  // clip index
    
  sse::int4v index = normalized.trunc_int();
  sse::double2v fraction= normalized - sse::double2v(index);
  //int i = index[0];
  //int j = index[1];
  int i = index.at<0>(); //index.get0();
  int j = index.at<1>(); //index.get1();
  int base_ix = i + j*nx_;
  sse::double2v result = bilinear_interpolate_sse(fraction, knots_sse[base_ix], knots_sse[base_ix+1], knots_sse[base_ix+nx_], knots_sse[base_ix+nx_+1]);

  result.store(pos, sse::unaligned);
}


#if 0
double2v
bilinear_interpolate(const array2d<double2v>& knots, int i, int j, const double2v& fraction) {
  double2v w10(0.0, 1.0);
  double2v w00(1.0, 1.0);
  double2v w01(1.0, 0.0);
  
  return (w00-fraction)*knots(i, j)
       + (    fraction)*knots(i+1, j+1)
       + (w10-fraction)*knots(i+1, j)
       + (w01-fraction)*knots(i  , j+1);

  double2v fraction_m1(double2v(1.0) - fraction);

  return fraction * knots(i+1, j+1)
    + fraction_m1 * knots(i  , j)
    + combine(fraction, fraction_m1) * knots(i+1, j)
    + combine(fraction_m1, fraction) * knots(i  , j+1);
}

bool
index_manip::compute(double& x, double& y) const {
  double2v xy(x, y);
  // double2v min_, inv_inc_, low_limit_, high_limit_
  double2v norm_ix  = (xy - min_)*inv_inc_;
  bool ok = (norm_ix >= low_limit_ && norm_ix_<high_limit_).all_true();
  if (ok) {
    int4v node(norm_ix);
    double2v fraction = norm_ix - node;
    double2v result = bilinear_interpolate(knots, node, fraction);
    ok = result.isnan().all_false();
  }
  double norm_x = (x-xmin_)*inv_xinc_;
  double norm_y = (y-ymin_)*inv_yinc_;

  bool ok = (norm_x >= 0.0) && (norm_x <= dnx_)
         && (norm_y >= 0.0) && (norm_y <= dny_);
  if (ok) {
    int norm_i = int(norm_x);
    int norm_j = int(norm_y);
    //int norm_i = double2int(norm_x);
    //int norm_j = double2int(norm_y);
    x=norm_x - norm_i;
    y=norm_y - norm_j;
  }
  return ok;
}
#endif

static index_manip_base *coord_normalize_tester;

int
index_normalize_test() {
  double x=real_sample_data[0], y=real_sample_data[1];
  return int(coord_normalize_tester->compute(x, y));
}

int
index_normalize_pos_test() {
  double pos[]={real_sample_data[0], real_sample_data[1]};
  reinterpret_cast<bilinear_index_manip_ssenan*>(coord_normalize_tester)->compute_pos(pos);
  return int(pos[0]);
}

static char *pointer;

int delete_test() {
  delete pointer;
  return 1;
}

template <template <size_t N> class HASH_TEST>
void run_hash_test(const char *format) {
  using namespace oyrke::test;
  using namespace oyrke::utility;
  double max_test_time = 0.5;
  run_adaptive_performance_test(sprintf(format, 2), max_test_time,  HASH_TEST<2>::run_it);
  run_adaptive_performance_test(sprintf(format, 3), max_test_time,  HASH_TEST<3>::run_it);
  run_adaptive_performance_test(sprintf(format, 4), max_test_time,  HASH_TEST<4>::run_it);
  run_adaptive_performance_test(sprintf(format, 5), max_test_time,  HASH_TEST<5>::run_it);
  run_adaptive_performance_test(sprintf(format, 6), max_test_time,  HASH_TEST<6>::run_it);
  run_adaptive_performance_test(sprintf(format, 7), max_test_time,  HASH_TEST<7>::run_it);
  run_adaptive_performance_test(sprintf(format, 8), max_test_time,  HASH_TEST<8>::run_it);
  run_adaptive_performance_test(sprintf(format, 9), max_test_time,  HASH_TEST<9>::run_it);
  run_adaptive_performance_test(sprintf(format, 10), max_test_time, HASH_TEST<10>::run_it);
}


void
init_testdata() {
  srand((unsigned)time(0));

  size_t value = 1231;
  for (int i = 0; i < HASH_SAMPLE_SIZE; ++i) {
    value += i;
    value = oyrke::algorithm::tiny_algo<1>::hash(&value);
    hash_sample_data[i] = int(value);
		float x = 2.0f*(float(rand()) / RAND_MAX - 0.5f);
    real_sample_data[i] = 1.0e3f * x * abs(x); 
  }
}


int
test_empty_compare_0() {
	return real_sample_data[0] != 0.0f;
}

int
read_32_floats_test() {
  using namespace oyrke::algorithm;
  //float sum = tiny_algo<float, 32>::accumulate(real_sample_data, 0.0f);
  //int sum = std::accumulate((int *)real_sample_data, (int *)real_sample_data+32, 0);
  int sum = tiny_algo<32>::accumulate((int*)real_sample_data, 0);
  return sum != 0;
}

template <typename ALIGNMENT, size_t offset>
int
read_32_floats_sse_test() {
  using namespace oyrke::algorithm;
  //float sum = tiny_algo<float, 32>::accumulate(real_sample_data, 0.0f);
  //int sum = std::accumulate((int *)real_sample_data, (int *)real_sample_data+32, 0);
  //int sum = tiny_algo<int, 32>::accumulate((int*)real_sample_data, 0);
  using namespace oyrke::algorithm::sse;
  
	ALIGNMENT align = static_cast<ALIGNMENT>(0);
	float *base = real_sample_data + offset;
  float4v sum(base, align);
  sum = sum + float4v(base+ 4, align);
  sum = sum + float4v(base+ 8, align);
  sum = sum + float4v(base+12, align);
  sum = sum + float4v(base+16, align);
  sum = sum + float4v(base+20, align);
  sum = sum + float4v(base+24, align);
  sum = sum + float4v(base+28, align);
  return sum[0] != 0.0f;
}

int
read_32_floats_sse_args4_test() {
  using namespace oyrke::algorithm;
  //float sum = tiny_algo<float, 32>::accumulate(real_sample_data, 0.0f);
  //int sum = std::accumulate((int *)real_sample_data, (int *)real_sample_data+32, 0);
  //int sum = tiny_algo<int, 32>::accumulate((int*)real_sample_data, 0);
  using namespace oyrke::algorithm::sse;
  int i=0;
  float4v sum(real_sample_data[i], real_sample_data[i+1],	real_sample_data[i+2], real_sample_data[i+3]); i+=4;
  sum += float4v(real_sample_data[i], real_sample_data[i+1],	real_sample_data[i+2], real_sample_data[i+3]); i+=4;
  sum += float4v(real_sample_data[i], real_sample_data[i+1],	real_sample_data[i+2], real_sample_data[i+3]); i+=4;
  sum += float4v(real_sample_data[i], real_sample_data[i+1],	real_sample_data[i+2], real_sample_data[i+3]); i+=4;
  sum += float4v(real_sample_data[i], real_sample_data[i+1],	real_sample_data[i+2], real_sample_data[i+3]); i+=4;
  sum += float4v(real_sample_data[i], real_sample_data[i+1],	real_sample_data[i+2], real_sample_data[i+3]); i+=4;
  sum += float4v(real_sample_data[i], real_sample_data[i+1],	real_sample_data[i+2], real_sample_data[i+3]); i+=4;
  sum += float4v(real_sample_data[i], real_sample_data[i+1],	real_sample_data[i+2], real_sample_data[i+3]); i+=4;
  return sum.at<0>() != 0.0f;
}

namespace sse = oyrke::algorithm::sse;

int
test_plus_sse() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y(&real_sample_data[4], sse::unaligned);
	sse::float4v z = x + y;
  return z.at<0>() != 0.0f;
}

int
test_multiply_sse() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y(&real_sample_data[4], sse::unaligned);
	sse::float4v z = x * y;
  return z.at<0>() != 0.0f;
}


int
test_divide_sse() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y(&real_sample_data[4], sse::unaligned);
	sse::float4v z = x + y;
  return z.at<0>() != 0.0f;
}



int
test_rcp() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = sse::math::reciprocal(x);
  return y.at<0>() != 0.0f;
}

int
test_rcp_nr() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = sse::math::reciprocal_nr(x);
  return y.at<0>() != 0.0f;
}

int
test_rcp_exact() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = 1.0f/x;
  return y.at<0>() != 0.0f;
}

int
test_rsqrt() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = sse::math::sqrt_reciprocal(x);
  return y.at<0>() != 0.0f;
}

int
test_rsqrt_nr() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = sse::math::sqrt_reciprocal_nr(x);
  return y.at<0>() != 0.0f;
}

int
test_rsqrt_exact() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = 1.0f/sse::math::sqrt(x);
  return y.at<0>() != 0.0f;
}

int
test_sqrt_exact() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = sse::math::sqrt(x);
  return y.at<0>() != 0.0f;
}


int
test_sin() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = sse::math::sin(x);
  return y.at<0>() != 0.0f;
}


int
test_cos() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = sse::math::cos(x);
  return y.at<0>() != 0.0f;
}

int
test_sincos() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v ysin;
	sse::float4v ycos;
	sse::math::sincos(x, ysin, ycos);
  return ysin.at<0>() != 0.0f;
}


int
test_cbrt() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = sse::math::cbrt(x);
  return y.at<0>() != 0.0f;
}

int
test_cbrt_halleys() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y = sse::math::cbrt_halleys(x);
  return y.at<0>() != 0.0f;
}


#if 0
int
test_sin_cephes() {
	float x=real_sample_data[0];
	sse::float4v y = sse::math::sin_cephes(x);
  return y.at<0>() != 0.0f;
}
#endif

int
test_sin_std() {
	float x=real_sample_data[0];
	sse::float4v y = std::sin(x);
	return y.at<0>() != 0.0f;
}

int
test_pow_int4v() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	
	sse::float4v z = sse::math::pow(x*0.01, sse::int4v(10));
	return z.at<0>() != 0.0f;
}

int
test_pow_int() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	
	sse::float4v z = sse::math::pow(x*0.01, 10);
	return z.at<0>() != 0.0f;
}


int
test_pow() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y(&real_sample_data[10], sse::unaligned);
	
	sse::float4v z = sse::math::pow(x*0.01, y*0.01);
	return z.at<0>() != 0.0f;
}

int
test_pow_int4v_int() {
	sse::int4v x(3);
	
	sse::int4v z = sse::math::pow(x, 10);
	return z.at<0>() != 0;
}

int
test_pow_int4v_int4v() {
	sse::int4v x(3);
	sse::int4v y(10);
	
	sse::int4v z = sse::math::pow(x, y);
	return z.at<0>() != 0;
}

int
test_pow_int_template() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v z = sse::math::pow<100>(1.0+x*1e-5);
	return z.at<0>() != 0.0f;
}


int
test_pow_int4v_int_template() {
	sse::int4v x(3);
	
	sse::int4v z = sse::math::pow<10>(x);
	return z.at<0>() != 0;
}


int
test_pow_v2() {
	sse::float4v x(&real_sample_data[0], sse::unaligned);
	sse::float4v y(&real_sample_data[10], sse::unaligned);
	
	sse::float4v z = sse::math::pow_v2(x*0.01, y*0.01);
	return z.at<0>() != 0.0f;
}

int
test_pow_std() {
	float x = real_sample_data[0];
	float y = real_sample_data[10]; 
	
	float z = std::pow(x*0.01f, y*0.01f);
	return z != 0.0f;
}

int
test_pow_int_std() {
	float x = real_sample_data[0];
	
	float z = std::pow(x*0.01f, 10);
	return z != 0.0f;
}


int
test_frexp_std() {
	float x=real_sample_data[0];
	int e;
	float y=std::frexp(x, &e);
	return y != 0.0f;
}


int
test_frexp_sse() {
	sse::float4v x(real_sample_data, sse::unaligned);
	sse::int4v e;
	sse::float4v y=sse::math::frexp(x, e);
	return y.at<0>() != 0.0f;
}

int
test_ldexp_std() {
	float x=real_sample_data[0];
	float z = std::ldexp(x, 10);
	return z != 0.0f;
}


int
test_ldexp_sse() {
	sse::float4v x(real_sample_data, sse::unaligned);
	sse::float4v y = sse::math::ldexp(x, sse::int4v(10));
	return y.at<0>() != 0.0f;
}


void
test_hash_performance() {
  init_testdata();

  using namespace oyrke::test;
//  run_adaptive_performance_test("Empty                 : ", max_test_time, empty_test);

  run_hash_test<hash_const_index>(   "Hash int[%2d]          : ");
  run_hash_test<hash_fixplain_index>("Hash int[%2d]  plain   : ");
  run_hash_test<hash_dynamic_index>( "Hash int[%2d]  dyn     : ");
  run_hash_test<hash_dynloop_index>( "Hash int[%2d] x100     : ");

  /*
  run_adaptive_performance_test("Hash int[2]          : ", max_test_time, test_tiny_algo_test_2);
  run_adaptive_performance_test("Hash int[3]          : ", max_test_time, test_tiny_algo_test_3);
  run_adaptive_performance_test("Hash int[6]          : ", max_test_time, test_tiny_algo_test_6);
  run_adaptive_performance_test("Hash int[7]          : ", max_test_time, test_tiny_algo_test_7);
  run_adaptive_performance_test("Hash int[8]          : ", max_test_time, test_tiny_algo_test_8);
  run_adaptive_performance_test("Hash int[9]          : ", max_test_time, test_tiny_algo_test_9);
  run_adaptive_performance_test("Hash int[10]         : ", max_test_time, test_tiny_algo_test_10);

  run_adaptive_performance_test("Hash int[2]  plain   : ", max_test_time, test_tiny_algo_test_2_plain);
  run_adaptive_performance_test("Hash int[3]  plain   : ", max_test_time, test_tiny_algo_test_3_plain);
  run_adaptive_performance_test("Hash int[6]  plain   : ", max_test_time, test_tiny_algo_test_6_plain);
  run_adaptive_performance_test("Hash int[7]  plain   : ", max_test_time, test_tiny_algo_test_7_plain);
  run_adaptive_performance_test("Hash int[8]  plain   : ", max_test_time, test_tiny_algo_test_8_plain);
  run_adaptive_performance_test("Hash int[9]  plain   : ", max_test_time, test_tiny_algo_test_9_plain);
  run_adaptive_performance_test("Hash int[10] plain   : ", max_test_time, test_tiny_algo_test_10_plain);

  run_adaptive_performance_test("Hash int[2]  dyn     : ", max_test_time, test_tiny_algo_test_2_dyn);
  run_adaptive_performance_test("Hash int[3]  dyn     : ", max_test_time, test_tiny_algo_test_3_dyn);
  run_adaptive_performance_test("Hash int[6]  dyn     : ", max_test_time, test_tiny_algo_test_6_dyn);
  run_adaptive_performance_test("Hash int[7]  dyn     : ", max_test_time, test_tiny_algo_test_7_dyn);
  run_adaptive_performance_test("Hash int[8]  dyn     : ", max_test_time, test_tiny_algo_test_8_dyn);
  run_adaptive_performance_test("Hash int[9]  dyn     : ", max_test_time, test_tiny_algo_test_9_dyn);
  run_adaptive_performance_test("Hash int[10] dyn     : ", max_test_time, test_tiny_algo_test_10_dyn);

  run_adaptive_performance_test("Hash int[3] /100     : ", max_test_time, test_tiny_algo_test_3_loop);
  run_adaptive_performance_test("Hash int[6] /100     : ", max_test_time, test_tiny_algo_test_6_loop);
  run_adaptive_performance_test("Hash int[7] /100     : ", max_test_time, test_tiny_algo_test_7_loop);
  run_adaptive_performance_test("Hash int[10]/100     : ", max_test_time, test_tiny_algo_test_10_loop);


  === Hasher w/ template <template <...> >
  Empty                 :      2193.739 ms   553852354        3.96 ns
  Hash int[ 2]          :      1986.278 ms   550778458        3.61 ns
  Hash int[ 3]          :      2013.426 ms   559920255        3.60 ns
  Hash int[ 4]          :      2005.417 ms   558005523        3.59 ns
  Hash int[ 5]          :      2033.657 ms   552255742        3.68 ns
  Hash int[ 6]          :      2013.656 ms   560094164        3.60 ns
  Hash int[ 7]          :      2030.307 ms   550898825        3.69 ns
  Hash int[ 8]          :      2010.678 ms   559035743        3.60 ns
  Hash int[ 9]          :      2063.559 ms   559296478        3.69 ns
  Hash int[10]          :      2010.494 ms   559571888        3.59 ns
  Hash int[ 2]  plain   :      2048.468 ms   259115341        7.91 ns
  Hash int[ 3]  plain   :      1987.979 ms   191928539       10.36 ns
  Hash int[ 4]  plain   :      1928.162 ms   148426725       12.99 ns
  Hash int[ 5]  plain   :      2054.134 ms   124599850       16.49 ns
  Hash int[ 6]  plain   :      1914.490 ms   100000000       19.14 ns
  Hash int[ 7]  plain   :      2015.585 ms    88336387       22.82 ns
  Hash int[ 8]  plain   :      2004.101 ms    78705135       25.46 ns
  Hash int[ 9]  plain   :      2051.806 ms    60388933       33.98 ns
  Hash int[10]  plain   :      1989.950 ms    50883936       39.11 ns
  Hash int[ 2]  dyn     :      2063.432 ms   420194558        4.91 ns
  Hash int[ 3]  dyn     :      1991.696 ms   369683828        5.39 ns
  Hash int[ 4]  dyn     :      2044.839 ms   302966272        6.75 ns
  Hash int[ 5]  dyn     :      2059.068 ms   238518584        8.63 ns
  Hash int[ 6]  dyn     :      2054.904 ms   222490356        9.24 ns
  Hash int[ 7]  dyn     :      2002.078 ms   196438415       10.19 ns
  Hash int[ 8]  dyn     :      1919.095 ms   168433125       11.39 ns
  Hash int[ 9]  dyn     :      2045.496 ms   150041895       13.63 ns
  Hash int[10]  dyn     :      2001.663 ms   141613347       14.13 ns
  Hash int[ 2] /100     :      2006.609 ms     6722339      298.50 ns
  Hash int[ 3] /100     :      2035.270 ms     4960044      410.33 ns
  Hash int[ 4] /100     :      2011.926 ms     3921721      513.02 ns
  Hash int[ 5] /100     :      2046.781 ms     2899245      705.97 ns
  Hash int[ 6] /100     :      1998.719 ms     2557661      781.46 ns
  Hash int[ 7] /100     :      1892.877 ms     2071666      913.70 ns
  Hash int[ 8] /100     :      2051.525 ms     1932059     1061.83 ns
  Hash int[ 9] /100     :      1994.788 ms     1735527     1149.38 ns
  Hash int[10] /100     :      1924.915 ms     1545498     1245.50 ns
  */
  ::scanf("\n");
}



void
test_basic_performance() {
  double max_test_time = 0.5;
  init_testdata();

  using namespace oyrke::test;
  run_adaptive_performance_test("Empty                 : ", max_test_time, empty_test);
  run_adaptive_performance_test("Integer add x1        : ", max_test_time, binary_function_test<std::plus<int>, 1>);
  run_adaptive_performance_test("Integer division x1   : ", max_test_time, binary_function_test<std::divides<int>, 1>);
  run_adaptive_performance_test("Integer modulo   x1   : ", max_test_time, binary_function_test<std::modulus<int>, 1>);
  run_adaptive_performance_test("Integer multiply x1   : ", max_test_time, binary_function_test<std::multiplies<int>, 1>);
  run_adaptive_performance_test("Integer division x100 : ", max_test_time, binary_function_test<std::divides<int>, 100>);
  run_adaptive_performance_test("Integer modulo   x100 : ", max_test_time, binary_function_test<std::modulus<int>, 100>);
  run_adaptive_performance_test("Integer multiply x100 : ", max_test_time, binary_function_test<std::multiplies<int>, 100>);
  run_adaptive_performance_test("Integer rshift x100   : ", max_test_time, binary_function_test<shift_right<int>, 100>);
  run_adaptive_performance_test("Integer divide/32 x100: ", max_test_time, unary_function_test<divide_const<int, 32>, 100>);
  run_adaptive_performance_test("Right shift 5    x100 : ", max_test_time, unary_function_test<right_shift_const<int, 5>, 100>);

  run_adaptive_performance_test("Double add   x1       : ", max_test_time, binary_function_test<std::plus<double>, 1>);
  run_adaptive_performance_test("Double division x1    : ", max_test_time, binary_function_test<std::divides<double>, 1>);
  run_adaptive_performance_test("Double multiply x1    : ", max_test_time, binary_function_test<std::multiplies<double>, 1>);

  run_adaptive_performance_test("Delete nullptr x1    : ", max_test_time, delete_test);

  run_adaptive_performance_test("add 32 floats                        : ", max_test_time, read_32_floats_test);
  run_adaptive_performance_test("add 32 floats SSE aligned            : ", max_test_time, read_32_floats_sse_test<sse::ALIGNED, 0>);
  run_adaptive_performance_test("add 32 floats SSE unaligned(0)       : ", max_test_time, read_32_floats_sse_test<sse::UNALIGNED, 0>);
  run_adaptive_performance_test("add 32 floats SSE unaligned(1)       : ", max_test_time, read_32_floats_sse_test<sse::UNALIGNED, 1>);
  run_adaptive_performance_test("add 32 floats SSE 4 args             : ", max_test_time, read_32_floats_sse_args4_test);

  run_adaptive_performance_test("compare x!= 0.0                      : ", max_test_time, test_empty_compare_0);
  run_adaptive_performance_test("x+y SSE                              : ", max_test_time, test_plus_sse);
  run_adaptive_performance_test("x*y SSE                              : ", max_test_time, test_multiply_sse);
  run_adaptive_performance_test("x/y SSE                              : ", max_test_time, test_divide_sse);
  run_adaptive_performance_test("reciprocal x.rcp                     : ", max_test_time, test_rcp);
  run_adaptive_performance_test("reciprocal NR                        : ", max_test_time, test_rcp_nr);
  run_adaptive_performance_test("reciprocal exact                     : ", max_test_time, test_rcp_exact);
  run_adaptive_performance_test("reciprocal sqrt                      : ", max_test_time, test_rsqrt);
  run_adaptive_performance_test("reciprocal sqrt NR                   : ", max_test_time, test_rsqrt_nr);
  run_adaptive_performance_test("reciprocal sqrt exact                : ", max_test_time, test_rsqrt_exact);
  run_adaptive_performance_test("sqrt exact                           : ", max_test_time, test_sqrt_exact);
  run_adaptive_performance_test("sin                                  : ", max_test_time, test_sin);
  run_adaptive_performance_test("cos                                  : ", max_test_time, test_cos);
  run_adaptive_performance_test("sincos                               : ", max_test_time, test_sincos);
  //run_adaptive_performance_test("sin cephes                           : ", max_test_time, test_sin_cephes);
  run_adaptive_performance_test("std::sin                             : ", max_test_time, test_sin_std);
  run_adaptive_performance_test("cbrt                                 : ", max_test_time, test_cbrt);
  run_adaptive_performance_test("cbrt halleys                         : ", max_test_time, test_cbrt_halleys);
  run_adaptive_performance_test("pow                                  : ", max_test_time, test_pow);
  run_adaptive_performance_test("pow v2                               : ", max_test_time, test_pow_v2);
  run_adaptive_performance_test("std::pow                             : ", max_test_time, test_pow_std);
  run_adaptive_performance_test("pow int 10                           : ", max_test_time, test_pow_int);
  run_adaptive_performance_test("pow int4v 10                         : ", max_test_time, test_pow_int4v);
  run_adaptive_performance_test("pow int4v^int 10                     : ", max_test_time, test_pow_int4v_int);
  run_adaptive_performance_test("pow int4v^int4v 10                   : ", max_test_time, test_pow_int4v_int4v);
  run_adaptive_performance_test("std::pow int 10                      : ", max_test_time, test_pow_int_std);
  run_adaptive_performance_test("sse::pow<100>(float4v)               : ", max_test_time, test_pow_int_template);
  run_adaptive_performance_test("sse::pow<10>(int4v)                  : ", max_test_time, test_pow_int4v_int_template);

	run_adaptive_performance_test("sse::frexp                           : ", max_test_time, test_frexp_sse);
  run_adaptive_performance_test("std::frexp                           : ", max_test_time, test_frexp_std);
  run_adaptive_performance_test("sse::ldexp                           : ", max_test_time, test_ldexp_sse);
  run_adaptive_performance_test("std::ldexp                           : ", max_test_time, test_ldexp_std);

#if 1
  index_manip *tester = new index_manip;
  coord_normalize_tester = tester;
  tester->nx_ = 22234;
  tester->ny_ = 33456;
  tester->inv_xinc_ = 0.02;
  tester->inv_yinc_ = 0.02;
  tester->xmin_ = 23456.7;
  tester->ymin_ = 43256.7;

  tester->init();
  run_adaptive_performance_test("Index normalize - index_manip                   : ", max_test_time, index_normalize_test);

  struct empty_index_manip : public index_manip_base {
    bool compute(double&, double&) const { return true; }
  };
  coord_normalize_tester = new empty_index_manip;
  run_adaptive_performance_test("Empty index normalize - empty_index_manip       : ", max_test_time, index_normalize_test);

  struct indirect_index_manip : public index_manip_base {
    const index_manip_base *next_;
    indirect_index_manip(const index_manip_base *p) : next_(p) {}
    bool compute(double& x, double& y) const { return next_->compute(x,y); }
  };

  coord_normalize_tester = new indirect_index_manip(coord_normalize_tester);
  run_adaptive_performance_test("Virtual x2 index normalize                      : ", max_test_time, index_normalize_test);
  coord_normalize_tester = new indirect_index_manip(coord_normalize_tester);
  run_adaptive_performance_test("Virtual x3 index normalize                      : ", max_test_time, index_normalize_test);
  coord_normalize_tester = new indirect_index_manip(coord_normalize_tester);
  run_adaptive_performance_test("Virtual x4 index normalize                      : ", max_test_time, index_normalize_test);
  coord_normalize_tester = new indirect_index_manip(coord_normalize_tester);
  run_adaptive_performance_test("Virtual x5 index normalize                      : ", max_test_time, index_normalize_test);
  coord_normalize_tester = new indirect_index_manip(coord_normalize_tester);
  run_adaptive_performance_test("Virtual x6 index normalize                      : ", max_test_time, index_normalize_test);
  coord_normalize_tester = new indirect_index_manip(coord_normalize_tester);
  run_adaptive_performance_test("Virtual x7 index normalize                      : ", max_test_time, index_normalize_test);
  coord_normalize_tester = new indirect_index_manip(coord_normalize_tester);
  run_adaptive_performance_test("Virtual x8 index normalize                      : ", max_test_time, index_normalize_test);

  {
    bilinear_index_manip *bitester = new bilinear_index_manip;
    bitester->nx_ = 50;
    bitester->ny_ = 50;
    bitester->inv_xinc_ = 1.0/20000.0;;
    bitester->inv_yinc_ = bitester->inv_xinc_;
    bitester->xmin_ = 23456.7;
    bitester->ymin_ = 43256.7;
    bitester->knots = (double *)real_sample_data;

    bitester->init();
    coord_normalize_tester = bitester;
    run_adaptive_performance_test("Real bilinear interp - bilinear_index_manip   : ", max_test_time, index_normalize_test);
  }
  {
    bilinear_index_manip_sse ssetester;
    ssetester.nx_ = 50;
    ssetester.ny_ = 50;
    ssetester.inv_xinc_ = 1.0/20000.0;;
    ssetester.inv_yinc_ = ssetester.inv_xinc_;
    ssetester.xmin_ = 23456.7;
    ssetester.ymin_ = 43256.7;
    ssetester.knots_sse = (oyrke::algorithm::sse::double2v *)real_sample_data;

    ssetester.init();
    coord_normalize_tester = &ssetester;
    run_adaptive_performance_test("Real bilinear interp SSE - bilinear_index_manip_sse : ", max_test_time, index_normalize_test);
  }
  {
    bilinear_index_manip_ssenan ssetester; // = new bilinear_index_manip_ssenan;
    ssetester.nx_ = 50;
    ssetester.ny_ = 50;
    ssetester.inv_xinc_ = 1.0/20000.0;;
    ssetester.inv_yinc_ = ssetester.inv_xinc_;
    ssetester.xmin_ = 23456.7;
    ssetester.ymin_ = 43256.7;
    ssetester.knots_sse = (oyrke::algorithm::sse::double2v *)real_sample_data;

    ssetester.init();
    coord_normalize_tester = &ssetester;
    run_adaptive_performance_test("Real bilinear interp SSENaN - bilinear_index_manip_ssenan: ", max_test_time, index_normalize_test);
    run_adaptive_performance_test("Real bilinear interp SSE[2] - bilinear_index_manip_ssenan: ", max_test_time, index_normalize_pos_test);

    double x=real_sample_data[0], y=real_sample_data[1];
    ssetester.knots_sse[ssetester.buffer_index(x, y)] = oyrke::algorithm::sse::double2v(std::numeric_limits<double>::quiet_NaN()).value();
    run_adaptive_performance_test("Real bilinear SSENaN - bilinear_index_manip_ssenan(NaN)  : ", max_test_time, index_normalize_test);
  }
#endif
}


