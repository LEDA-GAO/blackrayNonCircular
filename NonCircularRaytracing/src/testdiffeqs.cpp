#include "def.h"

VariablePair testdiffeqs(long double xobs, long double yobs, long double iobs)
{
    long double dobs, xobs2, yobs2, fact1, fact2, r02, r0, th0, phi0, s0, s02, kphi0;
    long double kt0, b, kr0, kth0, r, th, phi, kr, kth, fact3;
    long double b_v1, kt0_v1, kr0_v1, kth0_v1, kphi0_v1;
    long double c1,c2;
      long double met[4][4];
      array<long double, 5> diffs, vars; 
      std::array<std::array<long double, 4>, 4> metNC;
      const int max_prec = std::numeric_limits<long double>::max_digits10;
    dobs = 1.0e+5;    /* distance of the observer */
  xobs2 = xobs*xobs;
  yobs2 = yobs*yobs;
  
  fact1 = yobs*sin(iobs) + dobs*cos(iobs);
  fact2 = dobs*sin(iobs) - yobs*cos(iobs);
  
  r02 = xobs2 + yobs2 + dobs*dobs;
  
  r0 = sqrt(r02);
  th0 = acos(fact1/r0);
  phi0 = atan2(xobs,fact2);
  s0  = sin(th0);
  s02 = s0*s0;
			
  kr0 = dobs/r0;
  kth0 = -(cos(iobs) - dobs*fact1/r02)/sqrt(r02-fact1*fact1);
  kphi0 = -xobs*sin(iobs)/(xobs2+fact2*fact2);

  //test the old initial conditions 
   metric(r0, th0, met);
   //This fact3 is E 
    fact3 = sqrt(met[0][3]*met[0][3]*kphi0*kphi0-met[0][0]*(met[1][1]*kr0*kr0+met[2][2]*kth0*kth0+met[3][3]*kphi0*kphi0));
  
  kt0 = -(met[0][3]*kphi0+fact3)/met[0][0];

  // b = Lz/E. In Kerr metric, Eq(13)(14) in note
  
  b = -(met[3][3]*kphi0+met[0][3]*kt0)/(met[0][0]*kt0+met[0][3]*kphi0);
//find the constant of motion for the horizon penetrating metric (ingoing Kerr coordinate)
    metNC = metricNC(r0,th0);
  c1 = metNC[0][0]*(kt0+(r02+spin*spin)*kr0/(r02-2.0*r0+spin*spin))+metNC[0][1]*kr0;
  cout<<"c1 "<<std::setprecision(max_prec) <<c1<<endl;
  c2 = metNC[0][3]*(kt0+(r02+spin*spin)*kr0/(r02-2.0*r0+spin*spin))+metNC[1][3]*kr0+metNC[3][3]*(kphi0+spin*kr0/(r02-2.0*r0+spin*spin));
  cout<<"c2 "<<std::setprecision(max_prec) <<c2<<endl;

  c1 /= fact3;
  c2 /= fact3;


  
  kr0_v1 = kr0/fact3;
  kth0_v1 = kth0/fact3;
  kphi0_v1 = kphi0/fact3;
//below we assume an symptotically flat spherical coordinate result. 
  //  kt0 = sqrt(kr0*kr0+r02*kth0*kth0+r02*s02*kphi0*kphi0);

  // b = Lz/E. In Kerr metric, Eq(13)(14) in note. We can try to calculate the b from a flat spacetime. 
  
  //b = r02*s02*kphi0/kt0;
  
  //It is to change lambda' = -E lambda. In this case, kt0 = dt/(-E dlambda) = 1. Therefore, the time step h is -1, a negative value. fact3 = -E
  kr0 /= fact3;
  kth0 /= fact3;
  kphi0 /= fact3;
  //Q3: Why we can set the mass to 1 and only input the spin of the BH? 
			
  /* ----- carter constant ----- */
  //There is no carter constant for this metric 
  //c02 = 1. - s02;
			
  //carter = yobs2 - spin2*c02 + xobs2*c02;
  //carter = sqrt(carter);
			
  /* ----- solve geodesic equations ----- */
			
  r = r0;
  th = th0;
  phi = phi0;
	
  //ku = kt0 + kr0; // transform to the momentum in the ingoing Kerr coordinate
  kr = kr0;
  kth = kth0;

    vars[0] = r;
    vars[1] = th;
    vars[2] = phi;
    vars[3] = kr0;
    vars[4] = kth0;


  //diffs =   diffeqs(b,vars);
  diffs = diffeqsNC(c1,c2, vars);
    return {vars, diffs};
}