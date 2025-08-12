#include "def.h"
long double omegak(long double r)
{
    long double omega;
    long double th = Pi/2;
    long double met_rder[4][4];
    metric_rderivatives(r, th, met_rder);
    
    omega = ( -met_rder[0][3] + sqrt( (met_rder[0][3]*met_rder[0][3]) - (met_rder[0][0]*met_rder[3][3]) ))/(met_rder[3][3]);
    return omega;

}