#include "def.h"
std::array<std::array<long double, 4>, 4> metric_rderivativesNC(long double r, long double th)
{
    long double t2,t3,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t28,t29,t31,t34,t35,t36;
    long double dguudr,dgppdr,dgupdr;
    std::array<std::array<long double, 4>, 4> dmn;
    t9 = pow(lNP,4);
    t10 = pow(r,2);
    t6 = betaNC/2.;
    t11 = pow(spin,2);
    t12 = cos(th);
    t13 = pow(t12,2);
    t14 = t11*t13;
    t15 = t10 + t14;
    t16 = pow(t15,-3);
    t17 = t16*t9;
    t21 = pow(3,t6);
    t22 = pow(4,betaNC);
    t23 = pow(t17,t6);
    t24 = t21*t22*t23;
    t25 = 1 + t24;
    t29 = 1/t25;
    t31 = 1/t15;
    t2 = 2*betaNC;
    t3 = 1 + t2;
    t5 = pow(2,t3);
    t7 = 1 + t6;
    t8 = pow(3,t7);
    t18 = -1 + t6;
    t19 = pow(t17,t18);
    t26 = pow(t25,-2);
    t34 = sin(th);
    t35 = pow(t34,2);
    t28 = pow(t15,-2);
    t36 = t10 + t11;
    t20 = pow(t15,-5);
    dguudr = -4*t10*t28*t29 + 2*t29*t31 + betaNC*t10*t19*t20*t26*t5*t8*t9;
    dgppdr = -2*r*t28*t35*(-(t11*(t10 + t11 - 2*r*t29)*t35) + pow(t36,2)) + t31*t35*(4*r*t36 - t11*t35*(2*r - 2*t29 - (betaNC*t10*t19*t26*t5*t8*t9)/pow(t15,4)));
    dgupdr = 4*spin*t10*t28*t29*t35 - 2*spin*t29*t31*t35 - spin*betaNC*t10*t19*t20*t26*t35*t5*t8*t9;
    dmn[0][0] = dguudr;
    dmn[3][3] = dgppdr;
    dmn[0][3] = dgupdr;
    dmn[3][0] = dmn[0][3];
    return dmn;
}