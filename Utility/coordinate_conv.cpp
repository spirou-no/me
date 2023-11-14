   /*
    factor 	Prime factors.
    gcd	Greatest common divisor.            int4v gcd(int4v, int4v)
    isprime	True for prime numbers.         mask4v isprime(int4v)
    lcm	Least common multiple.              int4v lcm(int4v, int4v)
    nchoosek All combinations of n elements taken k at a time.  ???
    perms	All possible permutations.
    primes	Generate list of prime numbers.
    rat	Rational approximation.             void rat(float4v, eps, int4v& int4v&)
    rats	Rational output.

Coordinate System Transforms  Function
	Purpose
    cart2pol Transform Cartesian coordinates to polar.      void cart2pol(float4v x, float4v y, float4v& r, float4v& phi)
                                                            float4v cart2pol(float4v xy)
    cart2sph Transform Cartesian coordinates to spherical.  void cart2sph(float4v x,y,x, float4v& r, phi, theta)
                                                            float4v cart2sph(float4v xyz)
    pol2cart Transform polar coordinates to Cartesian.      void pol2cart(r phi, x y)
                                                            float4v pol2cart(float4v r_phi)
    sph2cart Transform spherical coordinates to Cartesian.  void sph2cart(r phi theta, x y z)
                                                            float4v sph2cart(float4v r_phi_theta)
*/
#include <utility/math.h>


namespace oyrke { namespace algorithm { namespace sse { namespace math {


    void 
    cartesian2polar(const float4v& x, const float4v& y, float4v& r, float4v& phi) {
        r = sqrt(x*x + y*y);  // hypot(x, y)
        phi = atan2(y, x);
    }

    void 
    polar2cartesian(const float4v& r, const float4v& phi, float4v& x, float4v& y) {
        float4v s, c;
        sincos(phi, s, c);
        x = r * c;
        y = r * s;
    }

    void 
    cartesian2spherical(const float4v& x,const float4v& y,const float4v& z, float4v& r, float4v& phi, float4v& theta) {
        float4v r_xy = sqrt(x*x + y*y);  // hypot(x, y);
        theta = atan2(z, r_xy);
        phi   = atan2(y, x);
        r     = sqrt(r_xy*r_xy + z*z);      // hypot(r_xy, z)
    }


    void 
    spherical2cartesian(const float4v& r,const float4v& phi,const float4v& theta, float4v& x, float4v& y, float4v& z) {
        float4v sphi, cphi, stheta, ctheta;
        sincos(phi, sphi, cphi);
        sincos(theta, stheta, ctheta);

        x = r * cphi * ctheta;
        y = r * sphi * ctheta;
        z = r * stheta;
    }
}}}}