#pragma once
#include <functional>

namespace oyrke { namespace algorithm {

template <typename T>
struct tautology {
    bool operator()(const T&) const { return true; }
};


template <typename T>
struct contradiction  {
    bool operator()(const T&) const { return false; }
};



template <typename T>
inline const T&
min3(const T& a, const T& b, const T& c) {
  return min(min(a,b), c);
}

template <typename T, typename Compare>
inline const T&
min3(const T& a, const T& b, const T& c, Compare comp) {
  return min(min(a,b,comp), c, comp);
}

template <typename T>
inline const T&
min4(const T& a, const T& b, const T& c, const T& d) {
  return min(min(a,b), min(c,d));
}

template <typename T, typename Compare>
inline const T&
min4(const T& a, const T& b, const T& c, const T& d, Compare comp) {
  return min(min(a,b,comp), min(c,d,comp), comp);
}

template <typename T>
inline const T&
max3(const T& a, const T& b, const T& c) {
  return max(max(a,b), c);
}

template <typename T, typename Compare>
inline const T&
max3(const T& a, const T& b, const T& c, Compare comp) {
  return max(max(a,b,comp), comp);
}

template <typename T>
inline const T&
max4(const T& a, const T& b, const T& c, const T& d) {
  return max(max(a,b), max(c,d));
}

template <typename T, typename Compare>
inline const T&
max4(const T& a, const T& b, const T& c, const T& d, Compare comp) {
  return max(max(a,b,comp), max(c,d,comp), comp);
}



template <typename T>
inline const T&
clip(const T& val, const T& lo, const T& hi) {
  // Require( lo <= hi )
  return val < lo
       ? lo
       : (hi < val) ? hi : val;
}



template <typename T, typename Compare>
inline const T&
clip(const T& val, const T& lo, const T& hi, Compare comp) {
  // Require( !comp(hi, lo) )
  return comp(val, lo)
       ? lo
       : (comp(hi,val) ? hi : val);
}


template <typename T>
inline bool
in_interval(const T& val, const T& lo, const T& hi) {
  // Require( lo <= hi )
  bool not_ok = val < lo ||  hi < val;
  return !not_ok;
}

template <typename T, typename Compare>
inline bool
in_interval(const T& val, const T& lo, const T& hi, Compare comp) {
  // Require( !comp(hi, lo) )
  bool not_ok = comp(val, lo) || comp(hi, val);
  return !not_ok;
}


// return sign of x.  Same as signum_interval(x, 0, 0).
template <typename T>
inline int
signum(const T& x) {
	return x > T(0)
		 ? 1
		 : x < T(0)  ? -1 : 0;
}

template <typename T>
inline int
signum_interval(const T& x, const T& lo, const T& hi) {
	return x < lo
		? -1
		: x > hi  ? 1 : 0;
}

template <typename T, typename Compare>
inline int
signum_interval(const T& x, const T& lo, const T& hi, Compare comp) {
  // Require( !comp(hi, lo) )
	return comp(x, lo) 
		? -1
		: comp(hi, x)  ? 1 : 0;
}



template <typename T>
inline T
sqr(const T& x) {
  return x * x;
}


template <typename T>
inline T
cube(const T& x) {
  return x * x * x;
}


template <typename Scalar>
inline int
nearest_int(Scalar x) {
  x += x > 0.0 ? 0.5 : -0.5;
  return int(x);
}


template <typename T>
inline void
order(T& a, T& b) {
  if (b < a) {
    swap(a,b);
  }
}


template <typename T, typename Compare>
inline void
order(T& a, T& b, Compare comp) {
  if (comp(b, a)) {
    swap(a,b);
  }
}



}}
