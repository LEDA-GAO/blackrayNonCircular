#include "def.h"
long double omegakNC(long double r)
{
    long double omega;
    long double th = Pi/2;
    std::array<std::array<long double, 4>, 4> met_rder;
    met_rder = metric_rderivativesNC(r, th);
    //compute the orbital velocity for the material in the accretion disk 
    omega = ( -met_rder[0][3] + sqrt( (met_rder[0][3]*met_rder[0][3]) - (met_rder[0][0]*met_rder[3][3]) ))/(met_rder[3][3]);
    return omega;

}