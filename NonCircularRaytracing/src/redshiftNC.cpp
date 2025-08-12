#include "def.h"
long double redshiftNC(long double r, long double kt0, long double kuf, long double krf, long double kpf)
{
    long double Omega,udt,g,ku,kp;
    long double th = Pi/2.0;
    std::array<std::array<long double, 4>, 4> metNC;
    Omega = omegakNC(r);
    metNC = metricNC(r,th);
    udt = 1.0/sqrt(-metNC[0][0] - (2.0*metNC[0][3]*Omega) - (metNC[3][3]*Omega*Omega));
    ku = metNC[0][0]*kuf+metNC[0][1]*krf+metNC[0][3]*kpf;
    kp = metNC[3][0]*kuf+metNC[3][1]*krf+metNC[3][3]*kpf;
    g = -kt0/(udt*ku+Omega*udt*kp);
    return g;

}
