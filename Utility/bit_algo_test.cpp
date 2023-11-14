#include <Utility/bit_algo.h>
#include <iostream>
#include <assert.h>

namespace oyrke { namespace algorithm { namespace test {


template <class Bits>
void
test_one() {
  //Bits koko;

  std::cout.width(10);
  std::cout << unsigned(Bits::value) << ":  bit count = ";

  std::cout.width(2);
  std::cout << Bits::bit_count << " MSB = ";

  std::cout.width(2);
  std::cout << Bits::msb << " LSB = ";

  std::cout.width(2);
  std::cout << Bits::lsb
    << " power of 2? " << (Bits::is_power_2 ? "yes" : "no ")
    << " prev power 2 exp= ";
  std::cout.width(2);
  std::cout << Bits::prev_power_2_exp << " next power 2 exp= ";
  std::cout.width(2);
  std::cout << Bits::next_power_2_exp << " prev power 2= ";
  std::cout.width(10);
  std::cout << unsigned(Bits::prev_power_2) << " next power 2= ";
  std::cout.width(10);
  std::cout << unsigned(Bits::next_power_2) << std::endl;


  assert(Bits::bit_count          == bit_count(Bits::value));
  assert(Bits::msb                == most_significant_bit(Bits::value));
  assert(Bits::lsb                == least_significant_bit(Bits::value));
  assert(Bits::is_power_2         == is_power_2(Bits::value));
  assert(Bits::prev_power_2_exp   == prev_power_2_exp(Bits::value));
  assert(Bits::next_power_2_exp   == next_power_2_exp(Bits::value));
  assert(Bits::prev_power_2       == prev_power_2(Bits::value));
  assert(Bits::next_power_2       == next_power_2(Bits::value));
}

void
bit_algo_test() {
  enum { max_test = (1u<<31) -1  };
  // enum { max_test = USHRT_MAX  };
  //enum { max_test = UINT_MAX };

  test_one<bit_algo<0> >();
  test_one<bit_algo<1> >();
  test_one<bit_algo<2> >();
  test_one<bit_algo<7> >();
  test_one<bit_algo<63> >();
  test_one<bit_algo<64> >();
  test_one<bit_algo<122> >();
  test_one<bit_algo<1024> >();
  test_one<bit_algo<1204> >();
  test_one<bit_algo<12040> >();
  test_one<bit_algo<max_test/2-1> >();
  test_one<bit_algo<max_test/2> >();
  test_one<bit_algo<max_test/2+1> >();
  test_one<bit_algo<max_test-2> >();
  test_one<bit_algo<max_test-1> >();
  test_one<bit_algo<max_test> >();
}

}}}