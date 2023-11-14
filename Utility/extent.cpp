#include <Utility/extent.h>
#include <Utility/algorithm.h>
#include <algorithm>
#include <utility>
#include <assert.h>

namespace {
  using namespace oyrke::basic;

  int gcd(int u, int v) {
    if (v > u) {
      std::swap(u, v);
    }

    while (v != 0) {
      u = u % v;
      std::swap(u, v);
    }

    return 0 <= u ? u : -u;
  }

  basic_extent<int>
  check_contains(const basic_extent<int> &lhs, const basic_extent<int> &rhs) {
    // At least one has size 1
    Require((lhs.size() == 1) || (rhs.size() == 1));

    basic_extent<int> res = basic_extent<int>::create_empty();
    if (( lhs.size() == 1 ) && rhs.size() == 1){
      if (lhs.front() == rhs.front()) {
        res = lhs;
      }
    }
    else if (lhs.size() == 1) {
      if (rhs.contains(lhs.front())) {
        res = lhs;
      }
    }
    else {
      if (lhs.contains(rhs.front())) {
        res = rhs;
      }
    }
    return res;
  }
}


namespace oyrke { namespace basic {

  basic_extent<int>
  make_intersection(const basic_extent<int> &lhs, const basic_extent<int> &rhs) {
    // Must do some more magic for float extents
    // compile_assert(numeric_limits<adt_Extent::value_type>::is_integer);

    basic_extent<int> res = basic_extent<int>::create_empty();
    if (lhs.empty() || rhs.empty()){
      return res;
    }

    if (lhs.size() == 1 || rhs.size() == 1){
      return check_contains(lhs, rhs);
    }

    // Create temp extents, make sure they have positive strides and
    // that t2 starts to the right of t1
    basic_extent<int> t1 = basic_extent<int>::copy_increasing(lhs);
    basic_extent<int> t2 = basic_extent<int>::copy_increasing(rhs);

    if (t2.front() < t1.front()) {
      std::swap(t1, t2);
    }

    basic_extent<int>::value_type n_inc = t1.stride() * t2.stride();
    n_inc /= gcd(t1.stride(), t2.stride());
    // search limited by ratio new stride / extent stride
    basic_extent<int>::iterator t1_begin = t1.nearest(t2.front());
    basic_extent<int>::iterator t1_end   = 1+t1_begin + n_inc/t1.stride();
    basic_extent<int>::iterator t2_end   = 1+t2.begin() + n_inc/t2.stride();

    // make sure we don't fall off the end.
    t1_end = std::min(t1_end, t1.end());
    t2_end = std::min(t2_end, t2.end());

    typedef std::pair<basic_extent<int>::iterator, basic_extent<int>::iterator> ext_IP;
    ext_IP IP = oyrke::algorithm::first_set_intersection(t1_begin, t1_end, t2.begin(), t2_end);

    // if t1 and t2 lie on the same "grid", and t2 starts within t1,
    // it now refers to that position.  Otherwise, it == t1.end().

    if (IP.first != t1_end) {
      basic_extent<int>::value_type n_beg = *IP.first;
      basic_extent<int>::value_type n_end = std::min(t1.back(), t2.back());

      res = basic_extent<int>(n_beg, 1 + (n_end - n_beg)/n_inc, n_inc);
      if (lhs.stride() < 0) {
        res = res.reverse();
      }
    }

    return res;
  }



  bool
  is_subset_of(
    const basic_extent<int> &lhs, 
    const basic_extent<int> &rhs
  ){
    return make_intersection(lhs, rhs).count() == lhs.count();
  }





  // Create an extent describing the indexes of a subset
  //   i.e  result = adt_IndexExtent(child, master)
  //        satisfies master[result[i]] == child[i]
  //                  and result.size() == child.size()
  // Returns an empty extent if child is not a subset of master or
  // any of the args is empty.
	basic_extent<int>
	make_index_extent(
		const basic_extent<int>& child,
		const basic_extent<int>& master
	) {
	extent::difference_type start = 0;
	extent::value_type inc   = 0;
	extent::size_type  size  = 0;

	basic_extent<int> res = basic_extent<int>::create_empty();
	if (!(child.empty() || master.empty())) {
		if ((child.size() == 1 ) || ( master.size() == 1 )) {
			if (child.size() == master.size()) {
				if (child.front() == master.front()) {
					start = 0;
					size  = 1;
					inc   = 1;
				}
			}
			else { 
				if (child.size() == 1 && (master.size() > 1)  && master.contains(child.front())) {
					start = master.find(child.front()) - master.begin();
					size  = child.size();
					inc   = child.stride() % master.stride() == 0 
						  ? child.stride() / master.stride() 
						  : 1;
				}
				// Else the master is smaller than the child
			}
		}
		else {
			if (child.stride() % master.stride() == 0 &&
				master.contains(child.front()) && 
				master.contains(child.back())) {
					start = master.find(child.front()) - master.begin();
					size  = child.size();
					inc   = child.stride() / master.stride();
			}
		}
	}
	return extent(static_cast<extent::value_type>(start), size, inc);
  }

  //
  // Return true iff lhs and rhs overlaps
  //
  bool
  is_intersecting(const basic_extent<int> &lhs, const basic_extent<int> &rhs) {
    basic_extent<int> t = make_intersection(lhs, rhs);
    return !t.empty();
  }

}}
