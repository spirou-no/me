#pragma once

#include <Utility/algobase.h>
#include <Utility/iterator.h>

//#include <chs_adt/adt_primitives.hh>
#define EXTENT_HAS_OSTREAM

#ifdef EXTENT_HAS_OSTREAM
#include <iostream>    // for operator<<
#endif

#include <algorithm>
#include <limits>
#include <assert.h>

//#include <Utilitychs_adt/adt_algobase.hh>
//#include <chs_adt/adt_iterator.hh>

#define Require assert

namespace oyrke { namespace basic {

template <class T>        class basic_extent;
template <class T, int N> class basic_kd_extent;

typedef basic_extent<int>           extent;
typedef basic_kd_extent<int, 2>        extent_2d;
typedef basic_kd_extent<int, 3>     extent_3d;
typedef basic_kd_extent<int, 4>     extent_4d;
typedef basic_extent<double>        real_extent;
typedef basic_kd_extent<double, 2>  real_extent_2d;
typedef basic_kd_extent<double, 3>  real_extent_3d;
typedef basic_kd_extent<double, 4>  real_extent_4d;


/* \brief
Vector type container with regular contents

Conceptually, an extent is a (possibly empty) sequence of integers, 
in either increasing or decreasing order, with a constant increment.
An basic_extent can be used anywhere a const vector<int> can be used,
i.e. the const mfs of vector are also available in basic_extent.

E = { e0, e1, e2, ..., en }

where e0 is the first element, en is the last, the count is (n+1),
and the increment is (e1-e0). The _position_ of element ei is i, and
we know that 0<= pos<n for all elements in E. If the count is 
zero, the increment and the elements e0 and en is undefined.

\note Some original mfs are kept for backwards compatibility
(-> end of public section).  Do not use these mfs in new code.
\note When changing code from old to new interface, remove the old mfs by
defining ADT_EXTENT_NO_OBSOLETE.
\note Declare some mfs inline even in class scope.
Otherwise, the SGI compiler complains that front, back and op[] are
redefined inline after use!!! - xy 20.05.98 
*/
template <class T>
class basic_extent {
public:
  typedef basic_extent<T>              self_type;
  typedef T                            value_type;
  typedef void                         pointer;    // something...
  typedef value_type                   const_reference;
  typedef const_reference              reference;
/*
  TODO iterator
  typedef adt_RndIterator<const self_type,
    value_type,
    pointer,
    const_reference> const_iterator;
  typedef const_iterator               iterator;
*/
  typedef ::ptrdiff_t									difference_type;
  typedef ::size_t										size_type;  // should be unsigned

  // dummy iterators 
  typedef random_index_iterator<const basic_extent, T, T> const_iterator;
  typedef const_iterator                               iterator;
  typedef std::reverse_iterator<iterator>              reverse_iterator;
  typedef reverse_iterator                             const_reverse_iterator;

/*
TODO iterator
#if defined(__STL_CLASS_PARTIAL_SPECIALIZATION) || defined(LINUX)
  // Since Solaris & sgi compilers doesn't support partial class
  // specializations yet, the iterator_traits based reverse_iterator
  // template class can't be used.
  typedef __STD::reverse_iterator<iterator>  const_reverse_iterator;
#else
  typedef __STD::reverse_iterator<iterator, value_type, const_reference,
    iterator::difference_type>  const_reverse_iterator;
#endif

  typedef const_reverse_iterator      reverse_iterator;
*/
  //
  // Creators
  //

  //! copy other, but make sure new extent has stride > 0
  static inline basic_extent<T> 
  copy_increasing(const self_type& other);

  //! create extent, reverse if needed so incr > 0
  static inline basic_extent<T>
  create_increasing(value_type beg, size_type count, value_type incr);

  //! create an empty extent
  static inline basic_extent<T> create_empty();

  //! Default constructor, creates an empty extent
  basic_extent() : begin_(0), size_(0), stride_(1) { }

  //! Create an extent over given range
  inline basic_extent(value_type begin, size_type sz, value_type inc);

#if 0
  basic_extent(const std::slice& s);
  // copy from standard slice.  Note that slice memebers are unsigned.
#endif

  //! Copy extent
  basic_extent(const self_type &rhs);

  //! Replace my extent with \c rhs
  self_type& operator=(const self_type &rhs);

  //! Return copy of me running in the opposite direction
  self_type    reverse() const;

  //! 2 extents are equal if they are empty, size() is 1 with equal front(),
  //! front(), stride() and() size() are equal.
  bool operator==(const self_type &rhs) const;

  bool operator!=(const self_type& rhs) const { 
    return !(*this == rhs);
  }

  //
  // extent "capacity"
  //

  //! Return current size of extent
  inline size_type size() const     { return size_; }

  //! Alias for \ref size.
  size_type count() const    { return size(); }

  //! Return \c true if extent is empty
  bool empty() const    { return size() == 0; }

  //! Extent size only limited by size's value type
  size_type max_size() const { return std::numeric_limits<size_type>::max(); }

  //! 1st value in extent. \pre !empty()
  value_type front() const;

  //! last value in extent, \pre !empty()
  value_type back() const;

  //! n-th element in extent.
  //! \pre 0 <= n < size()
  inline value_type operator[](size_t n) const;
  value_type get_at(size_t n) const { return (*this)[n]; }

  //
  // normal & reverse iterators
  //
  const_iterator         begin()  const { return iterator(this, 0); } 
  const_iterator         end()    const { return iterator(this, size_); }
  const_reverse_iterator rbegin() const { return reverse_iterator(end()); }
  const_reverse_iterator rend()   const { return reverse_iterator(begin()); }

  //! Returns iterator to req'd element, or end() if \c val not in extent
  iterator find(value_type val) const;

  //! Return iterator to nearest pos in extent, or end() if extent is empty
  iterator nearest(double x) const;

  //! Return true iff val is part of extent, i.e. find(val) != end()
  bool contains(value_type val) const;

  //! Current increment between extent elements
  value_type stride() const { return stride_; }

  //! Returns the begin value of the extent 
  value_type origin() const { return begin_;}

  //
  // modifiers
  //

  //! Set new start value for extent
  void set_front(value_type n_begin) { begin_ = n_begin; }

  //! resizes extent
  //! \pre n_size >= 0
  //! \note iterators are invalidated if new size > previous size
  inline void resize(size_type n_size);

  //! Set new stride/increment for extent
  // \pre n_stride != 0
  void set_stride(value_type n_stride);

#ifdef EXTENT_HAS_OSTREAM
  //! print me to ostream in style begin:stride:end (i.e. matlab style)
  std::ostream& print(std::ostream&) const;
#endif

private:
  value_type     begin_;
  value_type     stride_;
  size_type      size_;
};


/*
// These 3 guys are used 1 place each...  Should be free functions
bool isSubsetOf(const self_type &other) const ;
// chr_bfa/gen_synth.cc

bool intersectsWith(const self_type &other) const ;
// chr_lmf/lmf_dialog2.cc

int  closestValue(double x) const ;
// value_type snap(double x) const;
// chr_bfa/cubefunc.cc
*/



//
// basic_extent implementation of trivial mfs
//
template <class T>
inline basic_extent<T>
basic_extent<T>::create_empty() {
  return basic_extent<T>(T(0), 0, T(1));
}


template <class T>
inline basic_extent<T>::basic_extent(value_type begin, size_type sz, value_type inc)
  : begin_(begin), size_(sz), stride_(inc)
{
  Require(sz >= 0);
  Require(stride_ != 0);
}

template <class T>
inline basic_extent<T>::basic_extent(const self_type &rhs)
  : begin_(rhs.begin_), size_(rhs.size_), stride_(rhs.stride_)
{
}

#if 0
template <class T>
inline basic_extent::basic_extent(const std::slice& s)
  : begin_(s.start()), size_(s.size()), stride_(s.stride())
{
  // nothing.  s.size() is unsigned so no need to check size>=0
}
#endif

template <class T>
inline basic_extent<T>&
basic_extent<T>::operator=(const self_type &rhs) {
  new(this) basic_extent<T>(rhs);
  return *this;
}



template <class T>
inline
basic_extent<T>
basic_extent<T>::create_increasing(value_type beg, size_type count, value_type incr) {
  if (incr < 0) {
    beg += incr * (count - 1);
    incr = -incr;
  }
  return basic_extent(beg, count, incr);
}



template <class T>
inline
basic_extent<T>
basic_extent<T>::copy_increasing(const self_type& rhs) {
  return rhs.stride() > 0 ? rhs : rhs.reverse();
}



template <class T>
inline T
basic_extent<T>::front() const {
  Require(!empty());
  return begin_;
}


template <class T>
inline T
basic_extent<T>::back() const {
  Require(!empty());
  return (*this)[size()-1];
}

template <class T>
inline T
basic_extent<T>::operator[](size_t pos) const {
  Require(oyrke::algorithm::in_interval(pos, static_cast<size_type>(0), size()-1));
  return static_cast<T>(begin_ + pos * stride_);
}

template <class T>
inline void
basic_extent<T>::set_stride(value_type inc) {
  Require(inc != 0);
  stride_ = inc;
}


template <class T>
inline void
basic_extent<T>::resize(size_type s) {
  Require(s >= 0);
  size_ = s;
}


/*
** 2 extents are equal if one of the following is true
** o both extents are empty
** o both have count==1 and equal begin values
** o they have the equal count, begin and stride values.
*/
template <class T>
bool
basic_extent<T>::operator==(const self_type &rhs) const {
  bool eq = false;

  if (size() == rhs.size()) {
    eq = size()  == 0             
      || size()  == 1           && front()  == rhs.front()
      || front() == rhs.front() && stride() == rhs.stride(); 
  }

  return eq;
}



template <class T>
typename basic_extent<T>::iterator
basic_extent<T>::find(value_type val) const {
  size_type pos = (val - begin_) / stride_;
  value_type t = static_cast<value_type>(begin_ + pos * stride_);

  return t == val && !empty() && oyrke::algorithm::in_interval(pos, static_cast<size_type>(0), size()-1)
       ? begin() + pos
       : end();
}



template <class T>
bool
basic_extent<T>::contains(value_type val) const {
  return find(val) != end();
}


template <class T>
typename basic_extent<T>::iterator
basic_extent<T>::nearest(double x) const {
  basic_extent<T>::iterator result;

  if (empty()) {
    result = end();
  }
  else {
    double pos = (x - begin_) / stride_;
    size_type ipos   = oyrke::algorithm::nearest_int(pos);

    ipos       = oyrke::algorithm::clip(ipos, static_cast<size_type>(0), size()-1);
    result     = begin() + ipos;
  }

  return result;
}


template <class T>
basic_extent<T>
basic_extent<T>::reverse() const {
  return empty()
    ? basic_extent<T>(begin_, 0, -stride())
    : basic_extent<T>(back(), size(), -stride());
}

#ifdef EXTENT_HAS_OSTREAM
// Outputs an extent 4:32:2 meaning an extent running
// from 4 to 32 inclusive in increments of 2.
// Notation picked from matlab.
template <class T>
std::ostream&
basic_extent<T>::print(std::ostream& os) const {
  if (empty()) {
    os << "<empty>";
  }
  else {
    os << front() << ":" << stride() << ":" << back();
  }

  return os;
}

template <class T>
inline std::ostream&
operator<<(std::ostream& os, const basic_extent<T>& ext) {
  return ext.print(os);
}
#endif


basic_extent<int>
make_intersection(const basic_extent<int> &lhs, 
                  const basic_extent<int> &rhs);

bool
is_subset_of(const basic_extent<int> &lhs, 
             const basic_extent<int> &rhs);

bool
is_intersecting(const basic_extent<int> &lhs, 
                const basic_extent<int> &rhs);

// Give location of child in master
// Result is empty if child is not a subset of master (or either empty)
basic_extent<int>
make_index_extent(const basic_extent<int>& child,
                  const basic_extent<int>& master);

// Make a couple of transformation routines to help us out
template <class T>
double
extent_index_from_value(const basic_extent<T>& in, double val) {
  double res = (val - in.origin())/in.stride();
  return res;
}

template <class T>
double
extent_value_from_index(const basic_extent<T>& in, double index){
  double res = in.origin() + index * in.stride();
  return res;
}

template <class T, class U>
double
transform_extent_value(
  const basic_extent<T>& from, 
  const basic_extent<U>& to,
  double from_value
){
  double ind = extent_index_from_value(from, from_value);
  double to_val = extent_value_from_index(to, ind);
  return to_val;
}

/*!
Transform an extent by index transformation

The result extent has the same relationship to \c to_master as \c from_subset 
has to \c from_master.
It is assumed that the transformed extent is on the same grid as the
source extent from_master, i.e. if the from_subset and from_master extents where 
infinite, they would have a non-empty intersection.

\note The \c from_subset extent can be larger than \c from_master.
\note in.stride / src.stride() is assumed to be an integer.  I.e. for real
extents, the ratio is rounded to the nearest integer.
*/
template <class T, class U>
basic_extent<U>
transform_extent(
  const basic_extent<T>& from_subset, 
  const basic_extent<T>& from_master,
  const basic_extent<U>& to_master
) {
  int ratio = nearest_int(from_subset.stride() / from_master.stride());
  U loc_stride = to_master.stride() * ratio;

  // index is really an int
  double index = extent_index_from_value(from_master, from_subset.origin());
  U beg = to_master.origin() + nearest_int(index) * to_master.stride();
  return basic_extent<U>(beg, from_subset.size(), loc_stride);
}

/*******************************************************************************
** basic_kd_extent<T, N>
**
** Give an object with array behavior of a size N of basic_extents<T>
**
**
*******************************************************************************/

template <class T, int N>
class basic_kd_extent {
public:
  enum { X, Y, Z, W};
  enum { n_dims = N };

  typedef basic_kd_extent<T, N>    self_type;
  //   typedef basic_kd_extent<T, (int)less_one> minor_type;
  typedef basic_extent<T>         value_type;


  basic_kd_extent();
  basic_kd_extent(const value_type in[], int i=N); //  N-1 <= i <= N
  basic_kd_extent(const self_type& rhs);

  // assignment ok
  self_type& operator=(const self_type& rhs);

  bool operator==(const self_type& rhs) const;
  bool operator!=(const self_type& rhs) const {
    return !(*this == rhs);
  }

  self_type intersect(const self_type& rhs) const;
  self_type& self_intersect(const self_type& rhs);

  value_type& operator[](int i);
  const value_type& operator[](int i) const;

  bool empty() const;  // All extents are empty
  int  size() const { return N; }
  size_t kd_size() const ; //ext[0].size()*ext[1].size()..*ext[N-1].size()

  const value_type * data() const { return arr_; };


protected:
  int min_index() const { return 0; }
  int max_index() const { return N-1;}


private:
  value_type arr_[N];
};

/********************* IMPLEMENTATION ******************************/

template <class T, int N>
basic_kd_extent<T,N>::basic_kd_extent() {
  std::fill_n(arr_, N, value_type());
}

template <class T, int N>
inline
basic_kd_extent<T,N>::basic_kd_extent(const self_type& rhs) {
  std::copy(rhs.arr_, rhs.arr_ + N, arr_);
}

template <class T, int N>
inline basic_kd_extent<T,N>&
basic_kd_extent<T,N>::operator=(const self_type& rhs) {
  if (this != &rhs) {
    std::copy(rhs.arr_,rhs.arr_ + N, arr_);
  }
  return *this;
}



#if 0
template <class T, int N>
inline
basic_kd_extent<T,N>::basic_kd_extent(const minor_type& rhs)
{
  for( int i = 0; i < less_one ; i++){
    arr_[i] = rhs.arr_[i];
  }
  arr_[N-1] = basic_extent<T>();
}
#endif

template <class T, int N>
basic_kd_extent<T,N>::basic_kd_extent(const value_type  rhs[],
                                          int num)
{
  Require( num >= 0 && num <= N );
  std::copy(rhs,rhs + num, arr_);
  std::fill_n(arr_+num, N-num, value_type());
}


template <class T, int N>
bool
basic_kd_extent<T,N>::operator==(const self_type& rhs) const {
  return std::equal(arr_, arr_+N, rhs.arr_);
}


template <class T, int N>
basic_extent<T>&
basic_kd_extent<T,N>::operator[](int i) {
  Require(in_interval(i, min_index(), max_index()));
  return arr_[i];
}


template <class T, int N>
const basic_extent<T>&
basic_kd_extent<T,N>::operator[](int i) const {
  Require( in_interval(i, min_index(), max_index() ) );
  return arr_[i];
}


template <class T, int N>
bool
basic_kd_extent<T,N>::empty() const {
  bool is_empty = arr_[0].empty();
  for(int i = 1; (i < N ) && is_empty; i++){
    is_empty = arr_[i].empty();
  }
  return is_empty;
}


template <class T, int N>
size_t
basic_kd_extent<T,N>::kd_size() const {
  size_t loc_size = 1;
  for( int i = 0; i < N ; i++){
    loc_size *= arr_[i].size();
  }
  return loc_size;
}

template <class T, int N>
basic_kd_extent<T,N>&
basic_kd_extent<T,N>::self_intersect(const basic_kd_extent<T,N>& rhs) {
  std::transform(arr_, arr_+N, rhs.arr_, arr_, make_intersection);
  return *this;
}


template <class T, int N>
basic_kd_extent<T,N>
basic_kd_extent<T,N>::intersect(const basic_kd_extent<T,N>& rhs) const {
  basic_kd_extent<T,N> tmp(*this);
  tmp.self_intersect(rhs);

  return tmp;
}


}}