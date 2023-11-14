#pragma once

#include <Utility/algobase.h>

//#include <boost/random.hpp>

#include <algorithm>
#include <functional>
#include <iterator>
#include <utility>

namespace oyrke { 
    namespace algorithm {


namespace detail {
    template <typename T1, typename T2>
    struct compare_helper {
        bool operator()(const T1& lhs, const T2& rhs) const { return lhs < rhs; }
        bool operator()(const T2& lhs, const T1& rhs) const { return lhs < rhs; }
    };

    template <typename T>
    struct compare_helper<T, T> : public std::less<T> { };
    
    template <typename InIter, typename T>
    struct compare_picker {
        typedef compare_helper<typename std::iterator_traits<InIter>::value_type, T> comparer_type;
    };
}

template <class InputIterator1, class InputIterator2, class BinaryPredicate>
std::pair<InputIterator1, InputIterator2>
first_set_intersection(
  InputIterator1			first1, 
  InputIterator1			last1,
  InputIterator2			first2, 
  InputIterator2			last2,
  const BinaryPredicate&	compare
) {
  // __stl_debug_check(__check_range(first1, last1));
  // __stl_debug_check(__check_range(first2, last2));

  while (first1 != last1 && first2 != last2) {
    if (compare(*first1, *first2)) {
      ++first1;
    }
    else if (compare(*first2, *first1)) {
      ++first2;
    }
    else {
      break;
    }
  }

  return first1 != last1 && first2 != last2
    ? std::make_pair(first1, first2)
    : std::make_pair(last1, last2);
}


//! Returns the 1st common element in the 2 input sequences.  The result is the
// same as the first element output from set_intersection() with the same args.
// If no common element was found, the 2 "last" iterators are returned.
// 
template <class InputIterator1, class InputIterator2>
std::pair<InputIterator1, InputIterator2>
first_set_intersection(
  InputIterator1 first1, 
  InputIterator1 last1,
  InputIterator2 first2, 
  InputIterator2 last2
) {
  // return first_set_intersection(first1, last1, first2, last2, 
  // __stl_debug_check(__check_range(first1, last1));
  // __stl_debug_check(__check_range(first2, last2));
    typedef typename std::iterator_traits<InputIterator1>::value_type vt1;
    typedef typename std::iterator_traits<InputIterator2>::value_type vt2;
    
    return first_set_intersection(first1, last1, first2, last2,
                                  detail::compare_helper<vt1, vt2>());
}


template<typename InIter, typename UnaryOp, typename Pred> 
inline 
UnaryOp 
for_each_if(
    InIter			first, 
    InIter			last, 
    const UnaryOp&	func, 
    const Pred&		pred
) {
    for (; first != last; ++first) {
        if (pred(*first)) {
            func(*first);
        }
    }
    return func;
}


// TODO remove, its part of standard
template<typename FwdIter, typename OutIter, typename Pred> 
inline 
OutIter 
copy_if(
    FwdIter			first, 
    FwdIter			last, 
    OutIter			out, 
    const Pred&		pred
) {
    for (; first != last; ++first) {
        if (pred(*first)) {
            *out = *first;
            ++first;
            ++out;
        }
    }
    return out;
}



template<typename BiIter_1, typename BiIter_2, typename Pred> 
inline
BiIter_2 
copy_backward_if(
    BiIter_1		first, 
    BiIter_1		last, 
    BiIter_2		out_last, 
    const Pred&		pred
) {
    for (; first != last; --last) {
        if (pred(*last)) {
            --last;
            --out_last;
            *out_last = *last;
        }
    }
    return out_last;
}


// execute func whilst pred is true.  Return iterator to 1st non-processed iterator, or end
// if all elements were processed.
template<typename FwdIter, typename UnaryOp, typename Pred> 
inline 
FwdIter 
while_do(
    FwdIter		first, 
    FwdIter		last, 
    UnaryOp		func, 
    Pred		  pred
) {
    while (first != last && pred(*first)) {
        func(*first);
        ++first;
    }
    return first;
}



template<typename FwdIter, typename UnaryOp, typename Pred> 
inline 
FwdIter 
do_while(
    FwdIter		first, 
    FwdIter		last, 
    UnaryOp		func, 
    Pred		pred
) {
    if (first != last) {
        do {
            func(*first);
            ++first;
        } while (first != last && pred(*first));
    }
    return first;
}



template<typename FwdIter, typename OutIter, typename UnaryFunc, typename Pred> 
inline
OutIter 
transform_if(
    FwdIter		first, 
    FwdIter		last, 
    OutIter		out, 
    UnaryFunc	func, 
    Pred		pred
) {
    for (; first != last; ++first) {
        if (pred(*first)) {
            *out = func(*first);
            ++out;
        }
    }
    return out;
}



template<typename FwdIter_1, typename FwdIter_2, typename OutIter, typename BinaryFunc, typename Pred> 
inline
OutIter 
transform_if(
    FwdIter_1			first_1, 
    FwdIter_1			last_1, 
    FwdIter_2			first_2, 
    OutIter				out, 
    const BinaryFunc&	func, 
    const Pred&			pred
) {
    for (; first_1 != last_1; ++first_1, ++first_2) {
        if (pred(*first_1, *first_2)) {
            *out = func(*first_1, *first_2);
            ++out;
        }
    }
    return out;
}



template <typename FwdIter, typename BinaryPred, typename Pred>
inline
FwdIter 
min_element_if(
    FwdIter				first, 
    FwdIter				last, 
    const BinaryPred&	comp, 
    const Pred&			pred
) {
    bool    accept_first = true;
    FwdIter result       = last;

    for (; first != last; ++first) {
        if (pred(*first)) {
            if (accept_first || comp(*first, *result)) {
                result = first;
            }
            accept_first = false;
        }
    }
    return result;
}



template <typename FwdIter, typename Pred>
inline
FwdIter 
min_element_if(
    FwdIter		first, 
    FwdIter		last, 
    const Pred& pred
) {
    typedef typename std::iterator_traits<FwdIter>::value_type value_type;
    return min_element_if(first, last, std::less<value_type>(), pred);
}



template <typename FwdIter, typename BinaryPred, typename Pred>
inline
FwdIter 
max_element_if(
    FwdIter				first, 
    FwdIter				last, 
    const BinaryPred&	comp, 
    const Pred&			pred
) {

    bool    accept_first = true;
    FwdIter result       = last;

    for (; first != last; ++first) {
        if (pred(*first)) {
            if (accept_first || comp(*result, *first)) {
                result = first;
            }
            accept_first = false;
        }
    }
    return result;
}



template <typename FwdIter, typename Pred>
inline
FwdIter 
max_element_if(
    FwdIter		first, 
    FwdIter		last, 
    const Pred& pred
) {
    typedef typename std::iterator_traits<FwdIter>::value_type value_type;
    return max_element_if(first, last, std::greater<value_type>(), pred);
}


// create a contains function that returns true if
// the container contains the given value
//
template <typename InputIterator, typename T>
bool 
contains(
    InputIterator	first, 
    InputIterator	last, 
    const T&		value
) {
    InputIterator result = find(first, last, value);
    return result != last;
}


template <typename InputIterator, typename UnaryPred>
bool 
contains_if(
    InputIterator		first, 
    InputIterator		last, 
    const UnaryPred&	pred
) {
    InputIterator result = find_if(first, last, pred);
    return result != last;
}

namespace detail {
    template <typename BiIter, typename T, typename BinaryPredicate>
    std::pair<BiIter, BiIter>
    binary_expand_left(
        BiIter      first, 
        BiIter      hint_first, 
        typename std::iterator_traits<BiIter>::difference_type step, 
        const T&    value, 
        BinaryPredicate compare
    ) {
        typedef typename std::iterator_traits<BiIter>::difference_type diff_t;

        diff_t max_d = std::distance(first, hint_first);
        BiIter hint_last;
        bool   range_exhausted = false;
        bool   value_is_smaller = false;
        step = std::max(step, diff_t(1));

        do {
            step = std::min(step, max_d);
            hint_last = hint_first;
            std::advance(hint_first, -step);
            range_exhausted = step == max_d;
            value_is_smaller = compare(value, *hint_first);
            step = 2 * step;
        } while (value_is_smaller && !range_exhausted);

        return value_is_smaller 
             ? std::make_pair(first, first)  
             : std::make_pair(hint_first, hint_last);
    }


    template <typename BiIter, typename T, typename BinaryPredicate>
    std::pair<BiIter, BiIter>
    binary_expand_right(
        BiIter              hint_first, 
        BiIter              last, 
        typename std::iterator_traits<BiIter>::difference_type step, 
        const T&            value, 
        BinaryPredicate     compare
    ) {
        typedef typename std::iterator_traits<BiIter>::difference_type diff_t;
        
        diff_t max_d = max_d = std::distance(hint_first, last);
        BiIter hint_last = hint_first;
        bool   range_exhausted  = false;
        bool   value_is_smaller = false;
        step = std::max(step, diff_t(1));

        do {
            step = std::min(step, max_d);
            hint_first = hint_last;
            std::advance(hint_last, step);
            range_exhausted = step == max_d;
            value_is_smaller = !range_exhausted && compare(value, *hint_last);
            step = 2 * step;
        } while (!value_is_smaller && !range_exhausted);

        return value_is_smaller 
             ? std::make_pair(last, last)  
             : std::make_pair(hint_first, hint_last);
    }



    // returned range has p.first <= value < p.last
    template <typename BiIter, typename T, typename BinaryPredicate>
    std::pair<BiIter, BiIter>
    binary_expand(
        BiIter			first,
        BiIter			last,
        const T&		value,
        BiIter			hint,
        BinaryPredicate compare
    ) {

        std::pair<BiIter, BiIter> range = hint != last
                                        ? compare(value, *hint)
                                          ? detail::binary_expand_left(first, hint, 1, value, compare)
                                          : detail::binary_expand_right(hint, last, 1, value, compare)
                                        : std::make_pair(first, last);
        return range;
    }

        
    // returned range has p.first <= value < p.last
    template <typename BiIter, typename T, typename BinaryPredicate>
    std::pair<BiIter, BiIter>
    binary_expand(
        BiIter			first,
        BiIter			last,
        const T&		value,
        BiIter			hint_first,
        BiIter			hint_last,
        BinaryPredicate compare
    ) {
        
        if (hint_first == last) {
            return std::make_pair(first, last);
        }
        
        std::pair<BiIter, BiIter> range;
        if (compare(value, *hint_first)) {
            // someplace first..hint_first
            range = detail::binary_expand_left(first, hint_first, 
                                               2 * std::distance(hint_first, hint_last), 
                                               value, compare);
        }
        else if (hint_last != last) {
            range = compare(value, *hint_last)
                  ? std::make_pair(hint_first, hint_last)
                  : detail::binary_expand_right(hint_last, last, 
                                                2 * std::distance(hint_first, hint_last),  
                                                value, compare);
        }
        else {
            range = std::make_pair(hint_first, hint_last); // remember hint_last == last
        }
        return range;
    }
}



template <typename BiIter, typename T, typename BinaryPredicate>
BiIter
lower_bound_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint,
    BinaryPredicate compare
) {
    std::pair<BiIter, BiIter> range = detail::binary_expand(first, last, value, hint, compare);
    return std::lower_bound(range.first, range.second, value, compare);
}


template <typename BiIter, typename T>
BiIter
lower_bound_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint
) {
    return lower_bound_hinted(first, last, value, hint, 
                              detail::compare_picker<BiIter, T>::comparer_type());
}

template <typename BiIter, typename T, typename BinaryPredicate>
BiIter
lower_bound_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint_first,
    BiIter			hint_last,
    BinaryPredicate compare
) {
    std::pair<BiIter, BiIter> range = binary_expand(first, last, value, hint_first, hint_last, compare);
    return std::lower_bound(range.first, range.second, value, compare);
}

template <typename BiIter, typename T>
BiIter
lower_bound_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint_first,
    BiIter			hint_last
) {
    return lower_bound_hinted(first, last, value, hint_first, hint_last, 
                              detail::compare_picker<BiIter, T>::comparer_type());
}



template <typename BiIter, typename T, typename BinaryPredicate>
BiIter
upper_bound_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint,
    BinaryPredicate compare
) {
    std::pair<BiIter, BiIter> range = detail::binary_expand(first, last, value, hint, compare);
    return std::upper_bound(range.first, range.second, value, compare);
}


template <typename BiIter, typename T>
BiIter
upper_bound_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			seed
) {
    return upper_bound_hinted(first, last, value, seed, 
                              detail::compare_picker<BiIter, T>::comparer_type());
}

template <typename BiIter, typename T, typename BinaryPredicate>
BiIter
upper_bound_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint_first,
    BiIter			hint_last,
    BinaryPredicate compare
) {
    std::pair<BiIter, BiIter> range = binary_expand(first, last, value, hint_first, hint_last, compare);
    return std::upper_bound(range.first, range.second, value, compare);
}

template <typename BiIter, typename T>
BiIter
upper_bound_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint_first,
    BiIter			hint_last
) {
    return lower_bound_hinted(first, last, value, hint_first, hint_last, 
                              detail::compare_picker<BiIter, T>::comparer_type());
}



template <typename BiIter, typename T, typename BinaryPredicate>
std::pair<BiIter, BiIter>
equal_range_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint,
    BinaryPredicate compare
) {
    std::pair<BiIter, BiIter> range = detail::binary_expand(first, last, value, hint, compare);
    return std::equal_range(range.first, range.second, value, compare);
}


template <typename BiIter, typename T>
std::pair<BiIter, BiIter>
equal_range_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			seed
) {
    return equal_range_hinted(first, last, value, seed, 
                              detail::compare_picker<BiIter, T>::comparer_type());
}

template <typename BiIter, typename T, typename BinaryPredicate>
std::pair<BiIter, BiIter>
equal_range_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint_first,
    BiIter			hint_last,
    BinaryPredicate compare
) {
    std::pair<BiIter, BiIter> range = binary_expand(first, last, value, hint_first, hint_last, compare);
    return std::equal_range(range.first, range.second, value, compare);
}

template <typename BiIter, typename T>
std::pair<BiIter, BiIter>
equal_range_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint_first,
    BiIter			hint_last
) {
    return equal_range_hinted(first, last, value, hint_first, hint_last, 
                              detail::compare_picker<BiIter, T>::comparer_type());
}





template <typename BiIter, typename T, typename BinaryPredicate>
bool
binary_search_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint,
    BinaryPredicate compare
) {
    std::pair<BiIter, BiIter> range = detail::binary_expand(first, last, value, hint, compare);
    return std::binary_search(range.first, range.second, value, compare);
}


template <typename BiIter, typename T>
bool
binary_search_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			seed
) {
    return binary_search_hinted(first, last, value, seed, 
                                detail::compare_picker<BiIter, T>::comparer_type());
}

template <typename BiIter, typename T, typename BinaryPredicate>
bool
binary_search_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint_first,
    BiIter			hint_last,
    BinaryPredicate compare
) {
    std::pair<BiIter, BiIter> range = binary_expand(first, last, value, hint_first, hint_last, compare);
    return std::binary_search(range.first, range.second, value, compare);
}

template <typename BiIter, typename T>
bool
binary_search_hinted(
    BiIter			first,
    BiIter			last,
    const T&		value,
    BiIter			hint_first,
    BiIter			hint_last
) {
    return binary_search_hinted(first, last, value, hint_first, hint_last, 
                                detail::compare_picker<BiIter, T>::comparer_type());
}




/* Interface for random number generator
class RandomNumberGenerator : std::unary_function<size_t, size_t> {
    size_t operator()(size_t max);
}
*/

namespace detail {

    template <typename InIt, typename UnaryPred, typename RandomNumberGenerator>
    InIt 
    random_1_if(
        const InIt&         first, 
        const InIt&         last, 
        const UnaryPred&    pred,
        RandomNumberGenerator& rng
    ) {
        typedef typename std::iterator_traits<InIt>::difference_type diff_type;
        InIt result       = last;
        diff_type current = 0;

        for (InIt it = first; it != last; ++it) {
            if (pred(*it)) {
                diff_type index = rng(1+current);  
                bool accept = index == (current/4);  // any number in the range 0..current
                if (accept) {
                    result = it;
                }
                ++current;
            }
        }
        return result;
    }


    template <typename InIt, typename RandomNumberGenerator>
    inline InIt
    random_1_dispatch(
        const InIt& first, 
        const InIt& last, 
        RandomNumberGenerator& rng,
        std::input_iterator_tag
    ) {
        typedef typename std::iterator_traits<InIt>::value_type vt;
        return detail::random_1_if(first, last, tautology<vt>(), rng);
    }



    //! Pick a random element, specialization for random access iterators
    template <typename FwdIt, typename RandomNumberGenerator>
    FwdIt
    random_1_dispatch(
        FwdIt   first, 
        FwdIt   last, 
        RandomNumberGenerator& rng,
        std::forward_iterator_tag
    ) {
        // NOTE: As is, this algorithm will work on all iterators,but will require
        // 2 passes for non-random access iterators (use advance instead of +)
        // if the input size was passed as a parameter, then a single pass would do.
        typedef typename std::iterator_traits<FwdIt>::difference_type diff_type;
        diff_type max = std::distance(first, last);
        // diff_type index = last-first; //std::distance(first, last);
        diff_type index = rng(max);  // call random number generator
        std::advance(first, index);
        return first;
    }



    /*! Pick N random elements from input range.

    Since we don't know the size of the input range, the output must be rewritable, 
    i.e. random, bidir or forward iterator.
    Probabilities are adjusted based on how many inputs have been tested.
    Complexity: O(N) if out is random access, O(N*count) otherwise.
    */
    template <typename InIt, typename FwdIt, typename UnaryPred, typename RandomNumberGenerator>
    FwdIt
    random_n_if(
        InIt                    first, 
        InIt                    last, 
        size_t                  count, 
        FwdIt                   out, 
        const UnaryPred&        pred,
        RandomNumberGenerator&  rng
    ) {
        size_t current_input = 0;
        FwdIt  out_start  = out;

        for (InIt it = first; it != last; ++it) {
            if (pred(*it)) {
                if (current_input < count) {
                    *out = *it;
                    ++out;
                }
                else {
                    size_t test = rng(1+current_input);
                    if (test < count) {
                        // got count already, replace an arbitrary element in output
                        size_t replace_ix = rng(count);
                        FwdIt  tmp = out_start;
                        std::advance(tmp, replace_ix);
                        *tmp = *it;
                    }
                }

                ++current_input;
            }
        }

        return out;
    }


    template <typename InIt, typename FwdIt, typename RandomNumberGenerator>
    inline FwdIt
    random_n_dispatch(
        const InIt&             first, 
        const InIt&             last, 
        size_t                  count, 
        FwdIt                   out, 
        RandomNumberGenerator&  rng,
        std::input_iterator_tag /* input iterator tag */
    ) {
        typedef typename std::iterator_traits<InIt>::value_type vt;
        // assert out is forward or "better"
        return detail::random_n_if(first, last, count, out, tautology<vt>(), rng);
    }



    /*! Pick N random elements from input range.

    Since input range is given as a random access iterator range, output can be any
    kind of writable iterator, even a stream or inserter iterator.
    Probabilities are adjusted based on how many inputs and outputs are remaining.
    */
    template <typename InIt, typename OutIt, typename RandomNumberGenerator>
    OutIt
    random_n_from_m(
        const InIt&             first, 
        size_t                  input_size, 
        size_t                  count, 
        OutIt                   out,
        RandomNumberGenerator&  rng
    ) {
        typedef typename std::iterator_traits<InIt>::difference_type diff_type;
        size_t remaining_output = count;
        size_t remaining_input  = input_size; // std::distance(first, last);

        for (InIt it = first; remaining_input != 0 && remaining_output != 0; --remaining_input) {
            size_t test = rng(remaining_input);
            bool accept = test < remaining_output;
            if (accept) {
                --remaining_output;
                *out = *it;
                ++out;
            }
        }

        return out;
    }

    template <typename FwdIt, typename OutIt, typename RandomNumberGenerator>
    inline 
    OutIt
    random_n_dispatch(
        const FwdIt&            first, 
        const FwdIt&            last, 
        size_t                  count, 
        OutIt                   out,
        RandomNumberGenerator&  rng,
        std::forward_iterator_tag /* input_range_tag */
    ) {
        size_t input_size = std::distance(first, last);
        return random_n_from_m(first, input_size, count, out, rng);
    }
    
    
    // wrapper for system rand functions to fit into boost/TR1 random number framework
    class system_engine_rng {
    public:
        typedef int value_type;
        
        system_engine_rng() {}
        value_type operator()() { return ::rand(); }
        value_type min() const { return 0; }
        value_type max() const { return RAND_MAX; }
    };

    template <typename IntType>
    class system_rng {
    public:
        system_rng() : base_rng_() {} //, rng_(base_rng_) {}
        IntType operator()(IntType x) {
          const int RANDOM_BITS = 12;	// minimum random bits from rand() is 15, use 12 to be conservative
          const int RANDOM_MAX = (1U << RANDOM_BITS) - 1;

          IntType Rm = RANDOM_MAX;
          IntType Rn = ::rand() & RANDOM_MAX;
          while (Rm < x && Rm != ~0UL) {
            Rm = Rm << RANDOM_BITS | RANDOM_MAX;
            Rn = Rn << RANDOM_BITS | (::rand() & RANDOM_MAX);	// build random value
          }

          return Rn % x;
          //return rng_(x);
        }
    
    private:
        // no copy or assignment
        system_rng(const system_rng&);
        void operator=(const system_rng&);
        
        system_engine_rng base_rng_;
        //boost::random_number_generator<detail::system_engine_rng, IntType> rng_;
    };
}



/*! Return a random element in the given range.

Algorithm complexity is O(1) for random access iterator, O(N) otherwise.
*/
template <typename InIt, typename RandomNumberGenerator>
inline 
InIt
random_1(
    const InIt&             first, 
    const InIt&             last,
    RandomNumberGenerator&  rng
) {
    return detail::random_1_dispatch(first, last, rng, typename std::iterator_traits<InIt>::iterator_category());
}



template <typename InIt>
inline 
InIt
random_1(
    const InIt& first, 
    const InIt& last
) {
    detail::system_rng<typename std::iterator_traits<InIt>::difference_type> rng;
    return random_1(first, last, rng);
}


/*! Return a random element in the range satisfying the predicate.

Algorithm is O(N) for any iterator type.
*/
template <typename InIt, typename UnaryPred, typename RandomNumberGenerator>
inline InIt
random_1_if(
    const InIt&             first, 
    const InIt&             last, 
    const UnaryPred&        pred,
    RandomNumberGenerator&  rng
) {
    return detail::random_1_if(first, last, pred, rng);
}

template <typename InIt, typename UnaryPred>
inline InIt
random_1_if(
    const InIt&         first, 
    const InIt&         last, 
    const UnaryPred&    pred
) {
    detail::system_rng<typename std::iterator_traits<InIt>::difference_type> rng;
    return oyrke::algorithm::random_1_if(first, last, pred, rng);
}




/*! Pick N random elements from given range.

One of InIt or OutIt must be (at least) a forward iterator.
Complexity: O(N) if InIt or OutIt is random access
O(N*count) otherwise
*/
template <typename InIt, typename OutIt, typename RandomNumberGenerator>
inline 
OutIt
random_n(
    const InIt&     first, 
    const InIt&     last, 
    size_t          count, 
    OutIt           out,
    RandomNumberGenerator &rng
) {
    typedef typename std::iterator_traits<InIt>::iterator_category category;
    return detail::random_n_dispatch(first, last, count, out, rng, category());
}


template <typename InIt, typename OutIt>
inline 
OutIt
random_n(
    const InIt&     first, 
    const InIt&     last, 
    size_t          count, 
    OutIt           out
) {
    detail::system_rng<typename std::iterator_traits<InIt>::difference_type> rng;
    return random_n(first, last, count, out, rng);
}


/*! Pick N random elements from range, where range is given by first element and
count.

Same as the iterator range version, but main difference is relaxed requirements to input and
output iterators.  Input can be input_iterator, and the output can be an output_iterator
at the same time.
*/
template <typename InIt, typename OutIt, typename RandomNumberGenerator>
inline OutIt
random_n_of_m(
    const InIt& first, 
    size_t input_size, 
    size_t count, 
    OutIt out,
    RandomNumberGenerator &rng
) {
    return detail::random_n_from_m(first, input_size, count, out, rng);
}


template <typename InIt, typename OutIt>
inline OutIt
random_n_of_m(
    const InIt& first, 
    size_t input_size, 
    size_t count, 
    OutIt out
) {
    detail::system_rng<typename std::iterator_traits<InIt>::difference_type> rng;
    return random_n_from_m(first, input_size, count, out, rng);
}


/*! Pick N random elements from sequence.

Complexity. N is last-first.
output is random access: O(N)
output is bidir/fwd    : O(N*count)
*/
template <typename InIt, typename FwdIt_out, typename UnaryPred, typename RandomNumberGenerator>
inline 
FwdIt_out
random_n_if(
    const InIt& first, 
    const InIt& last, 
    size_t count, 
    FwdIt_out out, 
    const UnaryPred& pred,
    RandomNumberGenerator &rng
) {
    return detail::random_n_if(first, last, count, out, pred, rng);
}

template <typename InIt, typename FwdIt_out, typename UnaryPred>
inline 
FwdIt_out
random_n_if(
    const InIt& first, 
    const InIt& last, 
    size_t count, 
    FwdIt_out out, 
    const UnaryPred& pred
) {
    detail::system_rng<typename std::iterator_traits<InIt>::difference_type> rng;
    return oyrke::algorithm::random_n_if(first, last, count, out, pred, rng);
}


}}