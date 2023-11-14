#include <Utility/tiny_algo.h>
#include <functional>

namespace oyrke { namespace algorithm { namespace test {

  /* test                           int         double
  -------------------------------------------------------------------------------------------
  o find                            cr
  o find_if                         cr
  o count                           cr
  o equal                           cr
  o lexicographical_compare         cr
  o lexicographical_compare_3way    cr
  o hash                            cr
  o copy                            cr
  o fill                            cr
  o reverse                         cr
  o transform                       cr
  o transform bi                    cr
  o swap                            cr
  o min                             cr
  o min_element
  o max                             cr
  o max_element
  o accumulate                      cr
  o accumulate func                 cr
  o inner_product                   cr
  o inner_product func              cr
  o partial_sum                     cr
  o partial_sum func                cr  fail
  o iota                            cr




  */

	template <typename T, size_t N>
	size_t array_size(T (&)[N]) { return N; }

// local variable is initialized but not referenced
#pragma warning(disable: 4189)
  
void
tiny_algo_test() {
  typedef tiny_algo<3>   algo_3;
  int idata[3];
  int idata_2[3];

  size_t size = array_size(idata);
  algo_3::fill(idata_2, 5);
  algo_3::iota(idata, 2);  // fill w/ 2, 3, 4

	int y = eval_polynomial(idata, 2);

  int pos = algo_3::find(idata, 3); // pos = 1
  int no_pos = algo_3::find(idata, -9999);
  
  pos = algo_3::find_if(idata, [](auto x) { return x % 2; });

  size_t n = algo_3::count(idata, 3);
  bool equal = algo_3::equal(idata, idata);
  equal = algo_3::equal(idata, idata_2);

  size_t hashcode = algo_3::hash(idata);

  bool   less = algo_3::lexicographical_compare(idata, idata_2);
  less = algo_3::lexicographical_compare(idata, idata);
  int    cmp  = algo_3::lexicographical_compare_3way(idata, idata_2);
  cmp  = algo_3::lexicographical_compare_3way(idata, idata);

  algo_3::reverse(idata);
  algo_3::reverse(idata);

  int tmp[3];
  algo_3::copy(idata, tmp);

  algo_3::swap(tmp, idata_2);
  algo_3::swap(tmp, idata_2);

  double ddata[3];
  algo_3::transform(idata, ddata, [](double x) { return x + 3.0f; });
  algo_3::transform(idata, idata_2, ddata, std::plus<double>());

  int sum = algo_3::accumulate(idata);
  int product = algo_3::accumulate(idata, std::multiplies<int>(), 1);
  sum = algo_3::inner_product(idata, idata_2, 0);
  // sum_2 same as sum
  int sum_2 = algo_3::inner_product(idata, idata_2, std::plus<int>(), std::multiplies<int>());
  
  // swap plus and multiply, i.e. (x0+y0)*(x1+y1)*(x2+y2);
  sum_2 = algo_3::inner_product(idata, idata_2, std::multiplies<int>(), std::plus<int>(), 13);


  algo_3::partial_sum(idata, tmp);
  algo_3::partial_sum(idata, tmp, std::multiplies<int>());

  int low  = algo_3::min(idata);
  int high = algo_3::max(idata);

  //float fdata[3];
  //float_3_algo::iota(fdata, 4.0f);
  //hashcode = float_3_algo::hash(fdata);
}





}}}
