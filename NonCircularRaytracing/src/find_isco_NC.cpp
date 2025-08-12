#include "def.h"
long double find_isco_NC(long double r1, long double r2, long double h0, long double tol, int max_iter, long double hmin)
{
    //r1 and r2 are the lower and upper bounds for the root finder
    //The input h0 is the step size to decide the accuracy of the second derivative. 
    //Try to find the ISCO located at where the second derivative of the effective potential is 0 
    long double rEH, rlow, rhigh, rmax,h, rISCO;
    h = h0;
    rmax = 1000.0;
    rEH = 1.0+sqrt(1.0-spin*spin);

    while(h>hmin){
        rlow = r1;
        rhigh = r2;
        while(veffddr(rlow,h)*veffddr(rhigh,h)>=0){
            if(rlow>rEH) rlow/=2; 
            if(rhigh<rmax) rhigh *=2;
            if(rlow<=rEH && rhigh>=rmax){
                cout<<"cannot find ISCO!"<<endl;
                return 0.0;
            }
        };
    
        for(auto i=1;i<=max_iter;i++){
            rISCO = (rlow+rhigh)/2;
            if(abs(veffddr(rISCO,h))<tol) return rISCO;
            if(veffddr(rlow,h)*veffddr(rISCO,h)<0 ){
                rhigh = rISCO;
            }
            else{
                rlow = rISCO;
            }
        };
        h/=2;
    };
    cout<<"cannot find the ISCO radius satisfying the tolerance "<<tol<<endl;
    cout<<"second derivative of Veff = "<<abs(veffddr(rISCO,h))<<endl;
    return rISCO;
}