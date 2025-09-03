#include "def.h"
long double veffpotential(long double r)
{
    
    long double Omega,udt,const1,const2;
    long double th = Pi/2.0;
    std::array<std::array<long double, 4>, 4> metNC;
    Omega = omegakNC(r);
    metNC = metricNC(r,th);
    udt = 1.0/sqrt(-metNC[0][0] - (2.0*metNC[0][3]*Omega) - (metNC[3][3]*Omega*Omega));
    const1 = udt*(metNC[0][0]+metNC[0][3]*Omega);
    const2 = udt*(metNC[0][3]+metNC[3][3]*Omega);

    return const1;

}