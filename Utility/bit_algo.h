#pragma once

#include <stddef.h>         // size_t

namespace oyrke { namespace algorithm {

namespace test {
  void bit_algo_test();
}

template <size_t N>
class bit_algo {
  typedef bit_algo<N/2>    next_t;

  enum {
    is_bit0_set       = (N & 0x1) == 1,
  };
public:
  static const size_t value = N;

  enum States {
    is_zero           = N == 0,
//    value             = N,

    bit_count         = (is_bit0_set ? 1 : 0) + next_t::bit_count,
    msb               = 1+next_t::msb,
    lsb               = !is_bit0_set ? 1+next_t::lsb : 0,

    is_power_2        = bit_count == 1,
    // next_power_2_exp= next_m1_t::msb + 1,
    prev_power_2_exp  = msb,
    next_power_2_exp  = is_power_2 ? prev_power_2_exp : 1+prev_power_2_exp,

    prev_power_2      = 1u << prev_power_2_exp,
    next_power_2      = 1u << next_power_2_exp
  };
};

// template <>
template<>
class bit_algo<0u> {
public:
  enum {
    value             = 0,
    is_zero           = true,
    is_bit0_set       = false,

    bit_count         = 0,
    msb               = -1,
    lsb               = -1,

    is_power_2        = false,
    next_power_2_exp  = -1,
    prev_power_2_exp  = -1,
    next_power_2      = 0,
    prev_power_2      = 0
  };
};


int     most_significant_bit(size_t x);
int     least_significant_bit(size_t x);
size_t  bit_count(size_t x);
bool    is_power_2(size_t x);
int     next_power_2_exp(size_t x);
int     prev_power_2_exp(size_t x);
size_t  next_power_2(size_t x);
size_t  prev_power_2(size_t x);

}
}


