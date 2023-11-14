/*
Multiply row or column major matrix with vector
	float4v row_matrix_multiply(const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3, const float4v& x);
	float4v col_matrix_multiply(const float4v& c0, const float4v& c1, const float4v& c2, const float4v& c3, const float4v& x);

*/
#include <utility/math.h>

namespace oyrke { namespace algorithm { namespace sse { namespace math {

	float4v 
	row_matrix_multiply(
		const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3, 
		const float4v& x
	) {
		float4v y = sum_row(r0*x, r1*x, r2*x, r3*x);
		return y;
	}

	float4v 
	row_matrix_multiply_perspective(
		const float4v& r0, const float4v& r1, const float4v& r2, const float4v& r3, 
		const float4v& x
	) {
        float4v w = r3*x;
		float4v y = sum_row(r0*x, r1*x, r2*x, w);
        y /= w;
		return y;
	}

	
	
    //! Matrix multiply.  
    // Matrix is column ordered.  
    // x, y, z, w are all expanded, i.e. all 4 elements in x are the same.
	float4v 
	col_matrix_multiply(
		const float4v& c0, const float4v& c1, const float4v& c2, const float4v& c3, 
		const float4v& x,  const float4v& y,  const float4v& z,  const float4v& w
	) {
        float4v out = c0 * x   +   c1 * y   + c2 * z   + c3 * w;
		return out;
	}

	float4v 
	col_matrix_multiply(
		const float4v& c0, const float4v& c1, const float4v& c2, const float4v& c3, 
		const float4v& x
	) {
		//float4v y = c0*x.at<0>() + c1*x.at<1>() + c2*x.at<2>() + c3*x.at<3>();
        float4v y = col_matrix_multiply(c0, c1, c2, c3, 
                                        float4v::create_replicated<0>(x), 
                                        float4v::create_replicated<1>(x), 
                                        float4v::create_replicated<2>(x), 
                                        float4v::create_replicated<3>(x));
		return y;
	}

	void 
	matrix_transpose(float4v& rc0, float4v& rc1, float4v& rc2, float4v& rc3) {
		_MM_TRANSPOSE4_PS(rc0.value(), rc1.value(), rc2.value(), rc3.value());
	}


}}}}