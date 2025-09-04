#include "def.h"
#include <functional>
#include <stdexcept>

//try to use only double precision because the horizon only needs precision about 10e-3
//double horizon_delta(double r, double theta){
//    
//    double delta;
//    delta = r*r + spin*spin - 2.0*r/(1+pow((48*pow(lNP,4)/pow(r,6)),beta/2.0));
//    return delta;
//}

double horizon_delta(double r, double theta) {
    double cos_theta2 = cos(theta)*cos(theta);
    double r2 = r * r;
    double spin2 = spin*spin;
    double delta;

    //this is the delta = r^2+a^2-2*r*m(r,theta)
    delta = r2 + spin2 - 2.0*r/(1+pow((48*pow(lNP,4)/pow(r2+spin2*cos_theta2,3)),betaNC/2.0));
    
    // Check for valid expression under square root
    /* if (delta >= 0) {
        return 0.0;
    }
    */

    return delta;
}




const double EPS = 1e-12;

double find_root(
    const std::function<double(double)>& func,
    double x0, double x1, 
    double tol, 
    double bracket_size, int max_iter
) {
    
    // bisection (robust)
    double a = x0;
    double b = x1;
    double fa = func(x0);
    double fb = func(x1);
    
    if (fa * fb > 0) {
        // Expand bracket if needed
        double expand = 1.0;
        while (fa * fb > 0 && expand < 10.0) {
            a = x0 - bracket_size * expand;
            b = x0 + bracket_size * expand;
            fa = func(a);
            fb = func(b);
            expand += 0.5;
        }
        if (fa * fb > 0) {
            throw std::runtime_error("Cannot find horizon");
        }
    }
    
    for (int i = 0; i < max_iter; ++i) {
        double c = (a + b) / 2.0;
        double fc = func(c);
        
        if (std::abs(fc) < tol || std::abs(b - a) < EPS) {
            return c;
        }
        
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    
    throw std::runtime_error("Root finder failed to converge");
}


/*
long double find_horizon(long double r1, long double r2, int max_iter, long double tol){
    //find the horizon for the new metric
    long double rlow, rhigh, rH;
    rlow = r1;
    rhigh = r2;
    while(horizon_delta(rlow)*horizon_delta(rhigh)>=0){
        rlow/=2; 
        rhigh *=2;
    };
    
    for(auto i=1;i<=max_iter;i++){
        rH = (rlow+rhigh)/2;
        if(abs(horizon_delta(rH))<tol) return rH;
        if(horizon_delta(rlow)*horizon_delta(rH)<0 ){
            rhigh = rH;
        }
        else{
            rlow = rH;
        }
    };
    return rH;
}
*/



