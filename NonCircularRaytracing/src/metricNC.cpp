#include "def.h"
std::array<std::array<long double, 4>, 4> metricNC(long double r, long double th)
{
    long double t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15, t16;
    long double t17, t18, t19, t22, t23;
    long double guu, gthth, gpp, gur, grp, gup; /*metric components with lowered indicies*/
    std::array<std::array<long double, 4>, 4> mn;
    t2 = pow(r,2);
    t3 = pow(spin,2);
    t5 = cos(th);
    t6 = pow(t5,2);
    t7 = t3*t6;
    t8 = t2 + t7;
    t10 = betaNC/2.;
    t9 = 1/t8;
    t11 = pow(3,t10);
    t12 = pow(4,betaNC);
    t13 = pow(lNP,4);
    t14 = pow(t8,-3);
    t15 = t13*t14;
    t16 = pow(t15,t10);
    t17 = t11*t12*t16;
    t18 = 1 + t17;
    t19 = 1/t18;
    t22 = sin(th);
    t23 = pow(t22,2);
    guu = -1 + 2*r*t19*t9;
    gthth = t8;
    gpp = t23*(pow(t2 + t3,2) - t23*t3*(-2*r*t19 + t2 + t3))*t9;
    gur = 2;
    grp = -2*spin*t23;
    gup = -4*r*spin*t19*t23*t9;

    
    mn[0][0] = guu;
    mn[0][1] = gur/2.0;
    mn[1][0] = mn[0][1];
    mn[0][3] = gup/2.0;
    mn[3][0] = mn[0][3];
    mn[1][3] = grp/2.0;
    mn[3][1] = mn[1][3];
    mn[2][2] = gthth;
    mn[3][3] = gpp;
    return mn;
}
