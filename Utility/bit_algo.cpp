#include <Utility/bit_algo.h>

namespace oyrke { namespace algorithm {

int
most_significant_bit(size_t x) {
  int msb = -1;
  while (x != 0) {
    ++msb;
    x >>= 1;
  }

  return msb;
}


int
least_significant_bit(size_t x) {
  int lsb = x == 0 ? -1 : 0;
  while ((x != 0) && (x & 0x1) == 0) {
    ++lsb;
    x >>= 1;
  }

  return lsb;
}

size_t
bit_count(size_t x) {
  size_t count = 0;
  while (x != 0) {
    if ((x & 0x1) != 0) {
      ++count;
    }
    x >>= 1;
  }
  return count;
}


bool
is_power_2(size_t x) {
  return bit_count(x) == 1;
}


int next_power_2_exp(size_t x) {
  return x != 0 ? 1 + most_significant_bit(x-1) : -1;
}

int prev_power_2_exp(size_t x) {
  return most_significant_bit(x);
}

size_t
next_power_2(size_t x) {
  return x != 0 ? (1u << next_power_2_exp(x)) : 0;
}


size_t
prev_power_2(size_t x) {
  return x != 0 ? (1u << prev_power_2_exp(x)) : 0;
}

}}
