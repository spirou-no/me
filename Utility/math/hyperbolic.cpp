#include <utility/math.h>

#if 0
namespace oyrke { namespace algorithm { namespace sse { namespace math {

	float4v 
	cosh(const float4v& x) {
		float4v ex = exp(x);
		return 0.5f * (ex + 1.0f/ex);
	}

	float4v sinh(const float4v& x) {
		float4v ex = exp(x);
		return 0.5f * (ex - 1.0f/ex);
	}

	float4v tanh(const float4v& x) {
		// see cephes
		float4v ex = exp(x);
		float4v inv_ex = 1.0f/ex;
		return (ex - inv_ex) / (ex + inv_ex);  // sinh(x)/cosh(x)
	}

	float4v acosh(const float4v& x);																										// arc hyperbolic cosine
	float4v asinh(const float4v& x);																										// arc hyperbolic sine
	float4v atanh(const float4v& x);																										// arc hyperbolic tangent
}}}}

#endif