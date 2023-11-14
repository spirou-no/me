
/* Bessel function of noninteger order
 */
#include <utility/math.h>
#include <utility/sse.h>
#include <utility/math/constants.h>
#include <cmath>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

    float4v yv(float v, const float4v& x) {

        float4v y = select(isnan(x), x, y);  // if both v and x are NaN, x wins
        mask4v openarg = !isnan(x);

        float v_int = std::floorf(v);
        bool intarg = v == v_int;

        if (intarg) {
            float4v y_int = yn(static_cast<int>(v_int), x);
            y = select(openarg, y_int, y);
        }
        else {
            float4v sine, cosine;
            sincos(float4v(PI*v), sine, cosine);
            float4v y_float = (cosine * jv(v, x) - jv(-v, x)) / sine;
            y = select(openarg, y_float, y);
	    }

        return y;
    }
}}}}

