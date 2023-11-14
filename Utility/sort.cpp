/*
Sort elements in a vector, ascending or descending order.
*/

#include <utility/math.h>
#include <utility/sse.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

    float4v sort_4_asc(const float4v& x) {
        float4v p = shuffle<0,0,2,2>(x,x);
        float4v P = shuffle<1,1,3,3>(x,x);

        float4v r = min(p, P);
        float4v R = max(p, P);

        //r = shuffle<0,0,2,2>(r, r); -- not needed
        R = shuffle<3,3,1,1>(R, R);

        float4v s = min(r, R);
        float4v S = max(r, R);

        float4v st = shuffle<0,0,3,3>(s, S);
        float4v St = shuffle<3,3,0,0>(s, S);

        float4v t = min(st, St);
        float4v T = max(st, St);

        static const __m128i mask = {0xffffffff, 0x0, 0xffffffff, 0x0};
        //static mask4v mask = __m128i({0xffffffff, 0x0, 0xffffffff, 0x0});
        float4v y = select(mask4v(mask), t, T);
        return y;
    }
    
    // sort_{2,3,4}_{asc,desc}
}}}}