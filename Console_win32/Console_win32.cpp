// Console_win32.cpp : Defines the entry point for the console application.
//
//#undef _SECURE_SCL
//#define _SECURE_SCL 1
//#define _HAS_ITERATOR_DEBUGGING 1

#include "stdafx.h"
#include <utility/bit_algo.h>
#include <utility/stopwatch.h>
#include <utility/algobase.h>
#include <utility/performance_testing.h>
#include <utility/sprintf.h>

#include <algorithm>
#include <iterator>
#include <limits>
#include <stdio.h>
#include <stdarg.h>
#include <time.h>

#include <vector>
#include <utility/const_array.h>


extern void test_interpolate_performance();
extern void test_interpolate();
extern void test_hash_performance();
extern void test_basic_performance();
extern void test_sse();

/*
Test to compute union of N elements in a vector.
o elements are int, aggregate vector<int>
o elements are vector<int>, aggregate vector<int>

test impact of
o convert single element to aggregate_type (when value_type != aggregate_type
o copy single element
o use boost::reference_wrapper<T> to avoid actual copy
o use iterator version which holds reference to iterator.
*/

struct counted_int {
    int value;
    static int count_default;
    static int count_copy;
    static int count_assign;
    static int count_destruct;
    counted_int() : value(0) { ++count_default; }
    counted_int(const counted_int& rhs) : value(rhs.value) { ++count_copy; }
    counted_int& operator=(const counted_int& rhs) { value = rhs.value; ++count_assign; }
    ~counted_int() { ++count_destruct; }
    
    bool operator==(const counted_int& rhs) const { return value == rhs.value; }
    bool operator!=(const counted_int& rhs) const { return value != rhs.value; }
    bool operator< (const counted_int& rhs) const { return value < rhs.value; }
};


struct compute_traits {
    typedef std::vector<counted_int> aggregate_type;
    void combine(const counted_int& lhs, const counted_int& rhs, aggregate_type& result) {
        result.push_back(std::min(lhs, rhs));
        if (lhs != rhs) {
            result.push_back(std::max(lhs, rhs));
        }
    }
    
    void combine(const aggregate_type& lhs, const counted_int& rhs, aggregate_type& result) {
        oyrke::basic::const_array<counted_int> tmp(1, rhs);
        std::set_union(lhs.begin(), lhs.end(), tmp.begin(), tmp.end(), std::back_inserter(result));
        
        // copy lhs while *lhs < rhs
        // ignore lhs while *lhs == rhs
        // copy rhs
        // copy rest of lhs.
        
        // or wrap the single value rhs into a container type which provides iterators.
    }
    void combine(const aggregate_type& lhs, const aggregate_type& rhs, aggregate_type& result) {
        std::set_union(lhs.begin(), lhs.end(), rhs.begin(), rhs.end(), std::back_inserter(result));
    
    }
    
    void assign(const counted_int& source, aggregate_type& target) {
        target.push_back(source);
    }
    void assign(const aggregate_type& source, aggregate_type& target) {
        target = source;
    }
    void move(aggregate_type& source, aggregate_type& target) {
        target.swap(source);
    }
    
};

// local variable is initialized but not referenced
#pragma warning(push)
#pragma warning(disable: 4189)


#include <utility/tiny_algo.h>

void
test_tiny_algo() {
  oyrke::utility::stopwatch timer;
  oyrke::algorithm::test::tiny_algo_test();
  oyrke::utility::stopwatch::tick_t ticks = timer.elapsed_ticks();
  double time = timer.elapsed_time();
}

#pragma warning(pop)


class bit_statistics {
public:
  enum { word_size = std::numeric_limits<size_t>::digits };

  bit_statistics();
  void operator()(size_t x);

private:
  size_t n_samples_;
  size_t n_bits_;
  size_t bit_count_[1+word_size];
  size_t bit_field_[word_size];
};


bit_statistics::bit_statistics()
: n_samples_(0) {
  std::fill_n(bit_count_, int(1+word_size), 0);
  std::fill_n(bit_field_, int(word_size), 0);
}

void 
bit_statistics::operator()(size_t x) {
  ++n_samples_;
  size_t nbits = oyrke::algorithm::bit_count(x);
  n_bits_ += nbits;
  ++bit_count_[nbits];
  int bit_no = 0;
  while (x != 0) {
    if ((x & 0x1) != 0) {
      ++bit_field_[bit_no];
    }
    ++bit_no;
    x >>= 1;
  }
}

extern void random_match_test(int, int, bool);

static void
gather_hasher_statistics() {
  srand((unsigned)time(0));

  bit_statistics stats;
  int loops = 100000000;
  int max_rand = 1000;
  for (int i = 0; i < loops; ++i) {
    int pos[3];
    pos[0] = int(max_rand * double(rand()) / RAND_MAX);
    pos[1] = int(max_rand * double(rand()) / RAND_MAX);
    pos[3] = int(max_rand * double(rand()) / RAND_MAX);
    size_t hashcode = oyrke::algorithm::tiny_algo<3>::hash(pos);
    stats(hashcode);
  }
}

namespace oyvind {
  template <typename T, size_t N>
  struct foo {
    //enum { Size = N };
    static const int Size = N;
    T data_[N];

    T* begin() { return data_; }
    T* end() { return data_ + N; }
  };

  struct balle : public foo<float, 2> {
  };

}


#include <vector>
#include <map>
static void debugmap() {
  oyvind::foo<int, 5> koko;
  koko.data_[0] = 1;
  koko.data_[1] = 2;
  koko.data_[2] = 6;
  koko.data_[3] = 24;
  koko.data_[4] = 120;

  std::string xxx("hei");
  oyvind::balle kook;
  std::vector<float> bla(4, -2.0f);
  std::map<int, std::string> foo;
  foo[9]="I";
  foo[1]="A";
  foo[5]="E";
  foo[26]="Z";
}

int _tmain(int argc, _TCHAR* argv[]) {
//  oyrke::algorithm::test::bit_algo_test();
//  debugmap();
  test_sse();
	test_tiny_algo();
  test_interpolate();
  oyrke::utility::stopwatch timer;
  test_basic_performance();
 // test_hash_performance();
  test_interpolate_performance();
  //gather_hasher_statistics();
  
  ::scanf("\n");

#if 0
  srand(unsigned(oyrke::utility::stopwatch::ticks_now()));
  int nteams = argc > 1 ? atoi(argv[1]) : 12;
  int nreps  = argc > 2 ? atoi(argv[2]) : 2;
  bool print = argc > 3;

  random_match_test(nteams, nreps, print);
  oyrke::utility::stopwatch::tick_t now = timer.elapsed_ticks();
  timer.restart();
  now = timer.elapsed_ticks();
#endif

  return 0;
}


