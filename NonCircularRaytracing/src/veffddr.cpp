#include "def.h"
long double veffddr(long double r, long double h)
{
    //find the second derivative for the effective potential at radius r 
    //h is the step size 
    long double veffd2;
    veffd2 = (veffpotential(r+h)-veffpotential(r-h))/(2.0*h);

    return veffd2;

}