#include <boost/static_assert.hpp>

#include <algorithm>
#include <limits>

/*
y[i] = w[i] * x[i];
==> adt_TinySTL<float, 3>::transform(w, x, y, multiplies<float>());

y[i] = x[i];
==> adt_TinySTL::copy(x, y);

y[i] = c;
==> adt_TinySTL::fill(y, c);
*/


namespace oyrke { namespace algorithm {

namespace test {
  void tiny_algo_test();
}

template <int WordSize> 
struct hasher_constants_by_size {
  static const size_t seed        =        123u;
  static const size_t multiplier  = 1592763461u;
  static const size_t additive    =         13u;
};

// 
template <> 
struct hasher_constants_by_size<1> {
  static const size_t seed        =         13u;
  static const size_t multiplier  =         91u;
  static const size_t additive    =         13u;
};

template <> 
struct hasher_constants_by_size<2> {
  static const size_t seed        =         67u;
  static const size_t multiplier  =      27647u;
  static const size_t additive    =         13u;
};

template <typename T>
struct hasher_constants : public hasher_constants_by_size<sizeof(T)> { };

struct hasher_constants_2 {
  static const size_t seed        =        119u;
  static const size_t multiplier  = 1592763461u;
  static const size_t mult_2      =  910047527u;
  static const size_t invert      = 0x945da3b2u;
};


template <class T, int N>
class tiny_algo {
  BOOST_STATIC_ASSERT(0 < N);

  typedef tiny_algo<T, N-1> next_t;
  typedef tiny_algo<T, N-2> next_2_t;

public:
  static const int endpos = N;

  /*
  ** Initialization
  */

  static __forceinline void iota(T *dest, T value) {
    *dest = value;
    next_t::iota(dest+1, value+1);
  }

  static __forceinline T* copy(const T *src, T *dest) {
    copy_impl(src, dest);
    return dest+N;
  }


  static __forceinline void fill(T *dest, T value) {
    *dest = value;
    next_t::fill(1+dest, value);
  }

  /*
  Searching
  */
  static __forceinline int find(const T* src, T value) {  // position of 1st value
    return *src == value ? 0 : 1+next_t::find(src+1, value);
  }

  template <typename UnaryPred>
  static __forceinline int find_if(const T* src, const UnaryPred& pred) {  // position of 1st value
    return pred(*src) ? 0 : 1+next_t::find_if(src+1, pred);
  }

  static __forceinline size_t count(const T *src, T value) {
    return ((*src == value) ? 1 : 0) + next_t::count(src+1, value);
  }



  /*
  Comparison
  */
  static __forceinline bool equal(const T* src1, const T *src2) {
    return (*src1 == *src2) && next_t::equal(src1+1, src2+1);
  }

  static __forceinline bool lexicographical_compare(const T* src1, const T *src2) {
    return (*src1 < *src2)
      || (*src1 == *src2) && next_t::lexicographical_compare(src1+1, src2+1);
  }

  static __forceinline int lexicographical_compare_3way(const T* src1, const T *src2) {
    int result = *src1 < *src2
               ? -1
               : *src1 > *src2 ? 1 : next_t::lexicographical_compare_3way(1+src1, 1+src2);
    return result;
  }


  //typedef hasher_constants<T> magic_numbers;
  typedef hasher_constants_2 magic_numbers;

  static __forceinline /*inline*/ size_t hash(const T *src, size_t seed=magic_numbers::seed) {
    BOOST_STATIC_ASSERT(std::numeric_limits<T>::is_specialized);
    BOOST_STATIC_ASSERT(std::numeric_limits<T>::is_integer);

    //seed = magic_numbers::mult_1 * size_t(*src) + magic_numbers::mult_2 * (*src ^ seed);
    //seed = magic_numbers::multiplier * (size_t(*src) + seed) + magic_numbers::additive;
    seed = magic_numbers::multiplier * ((size_t(*src) ^ magic_numbers::invert) + seed);
    return next_t::hash(src+1, seed);
  }

  //! Polynomial evaluation by Horners rule
  static __forceinline T polynomial(const T *coeffs, T x) {
    return next_t::polynomial(coeffs+1, x) * x + coeffs[0];
  }


  static __forceinline void swap(T *lhs, T *rhs) {
    std::swap(*lhs, *rhs);
    next_t::swap(1+lhs, 1+rhs);
  }


  // todo: add specialization reverse<1>(T *) { } 
  static __forceinline void reverse(T *src) {
    enum { actual_N = N > 1 ? N : 0};
    typedef tiny_algo<T, actual_N> reverse_impl_t;

    reverse_impl_t::reverse_impl(src);
  }

  template <typename U, typename UnaryOp>
  static __forceinline void transform(const T *src, U *dest, const UnaryOp& op) {
    *dest = op(*src);
    next_t::transform(1+src, 1+dest, op);
  }

  template <typename U, typename BinaryOp>
  static __forceinline void transform(const T *src1, const T *src2, U *dest, const BinaryOp& op) {
    *dest = op(*src1, *src2);
    next_t::transform(1+src1, 1+src2, 1+dest, op);
  }

  // What if N==0?
  static __forceinline const T& min(const T *src) {
    return next_t::min_impl(*src, 1+src);
  }

  static __forceinline const T& max(const T *src) {
    return next_t::max_impl(*src, 1+src);
  }

  static inline int min_element(const T* src) {
  }

  static inline int max_element(const T* src) {
  }


  // const T& median(const T *src);
  // const T& nth_element(const T *src);

  static __forceinline T accumulate(const T *src, T init = T()) {
    return init + next_t::accumulate(1+src, *src);
  }

  template <class BinaryOp>
  static __forceinline typename BinaryOp::result_type
  accumulate(const T *src, const BinaryOp& op, T init = T()) {
    return op(init, next_t::accumulate(1+src, op, *src));
  }

  // same as
  // inner_product(src1, src2, multiplies<T>(), plus<T>());
  template <typename U, typename V>
  static __forceinline V inner_product(const T *src1, const U *src2, V init = V()) {
    return init +  next_t::inner_product(1+src1, 1+src2, *src1 * *src2); // TODO
  }

  template <class BinaryOp1, class BinaryOp2>
  static __forceinline T inner_product(const T *src1, const T *src2,
                                const BinaryOp1& op1, const BinaryOp2& op2, T init = T()) {
    return op1(init,
               next_t::inner_product(1+src1, 1+src2, op1, op2, op2(*src1, *src2)));
  }


  static __forceinline void partial_sum(const T *src, T *dest) {
    *dest = *src;
    next_t::partial_sum_impl(1+src, 1+dest, *src);
  }

  template <class BinaryOp>
  static __forceinline void partial_sum(const T *src, T *dest, const BinaryOp& op) {
    *dest = *src;
    next_t::partial_sum_impl(1+src, 1+dest, *src, op);
  }

public:
  static __forceinline void copy_impl(const T *src, T *dest) {
    *dest = *src;
    next_t::copy_impl(src+1, dest+1);
  }

  static __forceinline const T& min_impl(const T& x, const T* src) {
    return next_t::min_impl(std::min(x, *src), 1+src);
  }

  static __forceinline const T& max_impl(const T& x, const T* src) {
    return next_t::max_impl(std::max(x, *src), 1+src);
  }

  static __forceinline void reverse_impl(T *src) {
    std::swap(*src, *(src+N-1));
    next_2_t::reverse(1+src);
  }

  static __forceinline void partial_sum_impl(const T *src, T *dest, T init) {
    init += *src;
    *dest = init;
    next_t::partial_sum_impl(1+src, 1+dest, init);
  }

  template <class BinaryOp>
  static __forceinline void partial_sum_impl(const T *src, T *dest, T init, const BinaryOp& op) {
    init  = op(init, *src);
    *dest = init;
    next_t::partial_sum_impl(1+src, 1+dest, init, op);
  }
};



// specialize for 0
template <class T>
class tiny_algo<T, 0> {
public:
  static const int endpos = 0;

  static __forceinline int    find(const T*, T)                                  { return 0; }
  template <typename UnaryPred>
    static __forceinline int find_if(const T*, const UnaryPred&)                 { return 0; }
  static __forceinline size_t count(const T *, T)                                { return 0; }
  static __forceinline bool   equal(const T*, const T *)                         { return true; }
  static __forceinline bool   lexicographical_compare(const T*, const T *)       { return false; }
  static __forceinline int    lexicographical_compare_3way(const T*, const T *)  { return 0; }
  static __forceinline size_t hash(const T *, size_t seed)                       { return seed; }
  static __forceinline T polynomial(const T *coeffs, T x)                 { return coeffs[0]; }
  static __forceinline void   copy(const T *, T *)                               { }
  static __forceinline void   fill(T *, T)                                       { }

  // todo: add specialization reverse<1>(T *) { } 
  static __forceinline void reverse(T *)                                         { }

  template <typename U, typename UnaryOp>
    static __forceinline void transform(const T *, U *, const UnaryOp&)          { }
  
  template <typename U, typename BinaryOp>
    static __forceinline void transform(const T *, const T *, U *, const BinaryOp&) { }

  // also no-op for N==1
  static __forceinline void swap(T *, T *)                                       { }

  // What if N==0?  alternative compile_assert, make special case for N==1
  //{ return *src; }

  static __forceinline void min(const T *)                                       { }
  static __forceinline void max(const T *)                                       { }
  static __forceinline int min_element(const T*)                                 { return 0; }
  static __forceinline int max_element(const T*)                                 { return 0; }

  static __forceinline T accumulate(const T *, T init = T())                     { return init; }

  // tricky, identity value depends on op.  E.g. if multiply, then
  // return value should be T(1).
  template <class BinaryOp>
  static __forceinline T accumulate(const T *, const BinaryOp&, T init = T())    { return init; } 

  template <typename U, typename V>
  static __forceinline V inner_product(const T *, const U *, V init = V())       { return init; }

    // tricky return value, see accumulate.
  template <class BinaryOp1, class BinaryOp2>
  static __forceinline T inner_product(const T *, const T *,
                                const BinaryOp1&, const BinaryOp2&, T init = T()) { return init; }


  static __forceinline void partial_sum(const T *, T *)                          { }

  template <class BinaryOp>
    static __forceinline void partial_sum(const T *, T *, const BinaryOp&)       { }

  static __forceinline void iota(T *, T)                                         { }

public:
  static __forceinline void copy_impl(const T *, T *)                            { }
  static __forceinline const T& min_impl(const T& x, const T*)                   { return x; }
  static __forceinline const T& max_impl(const T& x, const T*)                   { return x; }
  static __forceinline void reverse_impl(T*)   { }
  static __forceinline void partial_sum_impl(const T *, T *, T) { }
  template <class BinaryOp>
    static __forceinline void partial_sum_impl(const T *, T *, T, const BinaryOp&) { }
};

// special cases for N==1

}
}
