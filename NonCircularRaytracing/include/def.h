#ifndef _DEF_H
#define _DEF_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <array>
#include <iomanip>      // For setprecision
#include <limits>       // For numeric_limits
#include <functional>
#include <stdexcept>

using namespace std;

#define imax 400

const long double Pi = acos(-1.0);

extern long double xobs, yobs;
extern long double epsi3, a13, a22, a52;
extern long double spin;
extern long double iobs_deg;
extern long double beta, lNP;

/*-----------------------------------------------------------*/
struct VariablePair {
    std::array<long double, 5> vars;
    std::array<long double, 5> diffs; 
};

struct IntegrationOut {
    std::array<long double, 3> varsOut;
    int stop_integration; 
};

const size_t HorizonN = 5001; // or your preferred size

template<size_t N>
struct HorizonOut {
    std::array<double, N> H;
    std::array<double, N> theta;
    // std::array<double, N> dHdth; // add if needed
};

// Declare global horizon variable
extern HorizonOut<HorizonN> global_horizon;

//void raytrace(long double xobs, long double yobs, long double iobs, long double xin, long double disk_length_combined ,long double traces[], int& stop_integration);
IntegrationOut raytracingNC(long double xobs, long double yobs, long double iobs, long double rin, long double disk_length_combined);
std::array<long double, 5> diffeqs(long double b, std::array<long double, 7> vars);
void redshift(long double r, long double ktkp, long double& gg);
//void redshift_polish_doughnut(long double r, long double th, long double l ,long double ktkp, long double& gg);
void intersection(long double x_1, long double y_1, long double z_1, long double x_2, long double y_2, long double z_2, long double x_d[]);
void metric(long double z1, long double z2, long double mn[][4]);
void metric_rderivatives(long double z1, long double z2, long double dmn[][4]);
void find_isco(long double z1, long double& isco);
VariablePair testdiffeqs(long double xobs, long double yobs, long double iobs);
std::array<long double, 5> diffeqsNC(long double c1, long double c2, std::array<long double, 5> vars);
std::array<std::array<long double, 4>, 4> metricNC(long double r, long double th);
std::array<std::array<long double, 4>, 4> metric_rderivativesNC(long double r, long double th);
long double omegak(long double r);
long double omegakNC(long double r);
long double veffpotential(long double r);
long double veffddr(long double r, long double h);
long double find_isco_NC(long double r1, long double r2, long double h0, long double tol, int max_iter, long double hmin);
long double redshiftNC(long double r, long double kt0, long double kuf, long double krf, long double kpf);
double horizon_delta(double r, double theta) ;
//long double find_horizon(long double r1, long double r2, int max_iter, long double tol);
double find_root(const std::function<double(double)>& func, double x0, double x1, double tol = 1e-6, double bracket_size = 0.1, int max_iter = 500);
long double three_point_interpolation(long double x, long double x0, long double x1, long double x2, long double f0, long double f1, long double f2);

class HighPrecisionRootFinder {
public:
    using Function = std::function<double(double)>;

    HighPrecisionRootFinder(Function f, double tol = 1e-30, int max_iter = 100)
        : f_(f), tol_(tol), max_iter_(max_iter) {}

    double newton_raphson(double x0) {
        double x = x0;
        double fx = f_(x);
        try{
        for (int i = 0; i < max_iter_; ++i) {
            if (std::abs(fx) < tol_) return x;
            double h = optimal_step_size(x);
            double dfdx = derivative(x, h);
            if (std::abs(dfdx) < tol_ * (std::abs(x) + 1)) {
                return brent(x - 0.1 * h, x + 0.1 * h);
            }
            double dx = fx / dfdx;
            double x_new = x - dx;
            if (std::abs(dx) < tol_ * (std::abs(x) + 1)) {
                return x_new;
            }
            x = x_new;
            fx = f_(x);
        }
        return x;
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error in Newton-Raphson method: " << e.what() << std::endl;
            return NAN; // Return the last computed value
        }
    }

    double brent(double a, double b) {
        double fa = f_(a);
        double fb = f_(b);
        if (fa * fb >= 0) {
            throw std::invalid_argument("Root not bracketed");
        }
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
        double c = a, fc = fa;
        bool mflag = true;
        double d = 0;
        for (int i = 0; i < max_iter_; ++i) {
            if (std::abs(fb) < tol_ || std::abs(b - a) < tol_) {
                return b;
            }
            double s;
            if (fa != fc && fb != fc) {
                s = a * fb * fc / ((fa - fb) * (fa - fc))
                   + b * fa * fc / ((fb - fa) * (fb - fc))
                   + c * fa * fb / ((fc - fa) * (fc - fb));
            } else {
                s = b - fb * (b - a) / (fb - fa);
            }
            double cond1 = (s - (3*a + b)/4) * (s - b);
            double cond2 = mflag && std::abs(s - b) >= std::abs(b - c)/2;
            double cond3 = !mflag && std::abs(s - b) >= std::abs(c - d)/2;
            double cond4 = mflag && std::abs(b - c) < tol_;
            double cond5 = !mflag && std::abs(c - d) < tol_;
            if (cond1 >= 0 || cond2 || cond3 || cond4 || cond5) {
                s = (a + b) / 2;
                mflag = true;
            } else {
                mflag = false;
            }
            double fs = f_(s);
            d = c;
            c = b;
            fc = fb;
            if (fa * fs < 0) {
                b = s;
                fb = fs;
            } else {
                a = s;
                fa = fs;
            }
            if (std::abs(fa) < std::abs(fb)) {
                std::swap(a, b);
                std::swap(fa, fb);
            }
        }
        return b;
    }

private:
    Function f_;
    double tol_;
    int max_iter_;

    double optimal_step_size(double x) const {
        constexpr double eps = std::numeric_limits<double>::epsilon();
        double scale = std::max(std::abs(x), 1.0);
        return std::sqrt(eps) * scale;
    }

    double derivative(double x, double h) const {
        double h2 = h * 2;
        return (-f_(x + h2) + 8*f_(x + h) - 8*f_(x - h) + f_(x - h2)) / (12 * h);
    }
};

template<size_t N>
HorizonOut<N> find_horizon(double dtheta){
    std::array<double, N> H;
    std::array<double, N> dHdth;
    std::array<double, N> theta;
    int Npts = H.size();
    if (Npts < 1) throw std::runtime_error("Npts must be >= 1");
    theta[0] = Pi/2;
    auto findInitialH = [&](double x) {
        return horizon_delta(x, theta[0]);
    };
    double rHkerr = 1 + sqrt(1 - spin*spin);

    double initial_guess = rHkerr;
    constexpr double tol = 1e-10;
    constexpr int max_iter = 100;

    HighPrecisionRootFinder root_finder(findInitialH, tol, max_iter);
    double rinitial = root_finder.newton_raphson(initial_guess);

    H[0] = rinitial;
    dHdth[0] = 0;
    for (int i = 1; i < Npts; ++i) {
        theta[i] = theta[i-1]+dtheta;
        auto findH = [&](double x) {
            double delta = horizon_delta(x, theta[i]);
            if (delta < 0) {
                return -x + H[i-1] + (dHdth[i-1] - sqrt(-delta)) * dtheta / 2.0;
            }
            else{
                return -x + H[i-1] + (dHdth[i-1]) * dtheta / 2.0; 
            }
        };

        // Use the previous H and dHdth to find the next H
        //H[i] = find_root(findH, H[i-1]-0.1, H[i-1] + dHdth[i-1] * dtheta+0.1, 1e-6, 0.1, 500); --- IGNORE ---
        HighPrecisionRootFinder root_finder(findH, tol, max_iter);
        H[i] = root_finder.newton_raphson(H[i-1] + dHdth[i-1] * dtheta);
        if(horizon_delta(H[i], theta[i]) < 0){
            dHdth[i] = - sqrt(-horizon_delta(H[i], theta[i]));
        }
        else{
            dHdth[i] = 0; 
        }
    }
    HorizonOut<N> horizon;
    horizon.H = H;
    //horizon.dHdth = dHdth;
    horizon.theta = theta;
    return horizon;
}


//void polish_doughnut(long double r, long double theta, long double phi ,long double angm_disk, long double spin, long double& w_current);
//void emission_angle_rth(long double r, long double th, long double l, long double a, long double kr, long double kphi ,long double kth , long double lambda , long double& em_angle);



#endif
