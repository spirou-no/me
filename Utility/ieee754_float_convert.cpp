/* Utility to convert between float representations that follow the IEEE-754 
convention of sign, a biased exponent and a mantissa.

Useful to convert between e.g. standard float formats and half-floats or compressed formats.
*/

#include <assert.h>

namespace oyrke { namespace algorithm { 

	class ieee754_float_descriptor {
		bool has_sign_;
		int exponent_bits_;
		int exponent_bias_;
		int mantissa_bits_;

	public:
		ieee754_float_descriptor(int exp_bits, int mantissa_bits);
		ieee754_float_descriptor(int exp_bits, int mantissa_bits, int exp_base, bool has_sign);

		bool has_sign() const { return has_sign_; }
		int exponent_bits() const { return exponent_bits_; }
		int exponent_bias() const { return exponent_bias_; }
		int mantissa_bits() const { return mantissa_bits_; }
	};

	ieee754_float_descriptor::ieee754_float_descriptor(int exp_bits, int mantissa_bits, int exp_bias, bool has_sign)
		: has_sign_(has_sign), exponent_bits_(exp_bits), exponent_bias_(exp_bias), mantissa_bits_(mantissa_bits)
	{}


	class ieee754_float_converter {
		ieee754_float_descriptor from_;
		ieee754_float_descriptor to_;

		typedef char from_type;
		typedef char result_type;

	public:
		ieee754_float_converter(const ieee754_float_descriptor& from, ieee754_float_descriptor& to);
		result_type operator()(const from_type& x) const;

	private:
		static void unpack_float(const ieee754_float_descriptor& desc, const from_type& x, bool &sign, int& expnt, int& mantissa);
		static result_type pack_float(const ieee754_float_descriptor& desc, bool sign, int expnt, int mantissa);
	};


	ieee754_float_converter::ieee754_float_converter(
		const ieee754_float_descriptor& from, 
		ieee754_float_descriptor& to)
		: from_(from), to_(to) {}

	void
	ieee754_float_converter::unpack_float(
		const ieee754_float_descriptor& desc, 
		const ieee754_float_converter::from_type& x, 
		bool &sign, 
		int& expnt, 
		int& mantissa
	) {
	}

	ieee754_float_converter::result_type
	ieee754_float_converter::pack_float(
		const ieee754_float_descriptor& desc, 
		bool sign, 
		int expnt, 
		int mantissa
	) {
		assert(desc.has_sign() || !desc.has_sign() && !sign);
		result_type t = 0;
		if (desc.has_sign()) {
			t &= (sign ? 1 : 0);
		}
		t <<= desc.exponent_bits();
		if (expnt > 0 && expnt < 0xfff) {
			expnt -= desc.exponent_bias();
		}

		t &= expnt;
		t <<= desc.mantissa_bits();
		t &= mantissa;

		return t;
	}

	ieee754_float_converter::result_type 
	ieee754_float_converter::operator()(const ieee754_float_converter::from_type& x) const {
		bool from_sign;
		int from_exp;
		int from_mantissa;

		bool to_sign = false;
		int to_exp = 0;
		int to_mantissa = 0;

		unpack_float(from_, x, from_sign, from_exp, from_mantissa);
		// check inf or nan or zero
		// ...
		if (from_exp == 0xffff) {
			// check inf or nan or zero
		}
		else if (from_sign && !to_.has_sign()) {
			to_sign = false;
			to_exp = 0;
			to_mantissa = 0;
		}
		else {

			// subtract 1 and 2 respectively since 0x0..0 is reserved for zero and subnormals
			// and 0xf..f is reserved for infinity and NaN
			int to_exp_min = 1 - to_.exponent_bias();
			int to_exp_max = (1 << to_.exponent_bits()) - 2 - to_.exponent_bias();
			to_sign = from_sign;

			// check for underflow
			if (from_exp < to_exp_min) {
				// can be represented as subnormal?

				if (from_exp < to_exp_min + to_.mantissa_bits()) {
					// too small for subnormal, result is 0. copy sign from input.
					to_exp = 0;
					to_mantissa = 0;
				}
				else {
					to_exp = 0;  // subnormal
					to_mantissa = from_mantissa >> (to_exp_min - from_exp);
				}
			}
			else if (from_exp > to_exp_max) {
				// overflow --> infinity
				to_exp = 0xffffff;
				to_mantissa = 0;
			}
			else {
				// normal case
				to_exp = from_exp;
				to_mantissa = from_mantissa;
			}
		}

		result_type y = pack_float(to_, to_sign, to_exp, to_mantissa);
		return y;
	}


}}