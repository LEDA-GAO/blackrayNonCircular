#include "def.h"
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>

// Define the global horizon variable
HorizonOut<HorizonN> global_horizon;

// Find rHorizon for exact match
double get_rHorizon_for_spin(double spin, const std::vector<double>& spin_vals, const std::vector<double>& rHorizon_vals) {
    for (size_t i = 0; i < spin_vals.size(); ++i) {
        if (std::abs(spin_vals[i] - spin) < 1e-8) { // tolerance for floating point
            return rHorizon_vals[i];
        }
    }
    throw std::runtime_error("Spin value not found in spin_vals");
}

// Or: Find rHorizon for nearest spin
double get_rHorizon_nearest(double spin, const std::vector<double>& spin_vals, const std::vector<double>& rHorizon_vals) {
    if (spin_vals.empty()) throw std::runtime_error("spin_vals is empty");
    size_t best = 0;
    double min_diff = std::abs(spin_vals[0] - spin);
    for (size_t i = 1; i < spin_vals.size(); ++i) {
        double diff = std::abs(spin_vals[i] - spin);
        if (diff < min_diff) {
            min_diff = diff;
            best = i;
        }
    }
    return rHorizon_vals[best];
}

int main(int argc, char *argv[])
{
    long double alpha;
    long double iobs;
    long double robs, pobs;
    long double robs_i, robs_f, rstep, pstep;
    long double xin, xout;
    long double isco, iscoNC,risco;
    long double omega,omegaNC,gfactor;
    long double rHkerr;
    long double rHNC, rH;
    long double isco_initial_guess;
    long double pp, qq, fr;
    double r;
    double h = 1e-7;
    double horizon_limit;
    //long double traced[5];    
    int stop_integration_condition = 0;
    int photon_index = 0, index, i;
    std::array<long double, imax> E_obs{}, N_obs{}, fphi{}; 

    const int max_prec = std::numeric_limits<long double>::max_digits10;
    char filename_o[128];
    char filename_o2[128];
    //std::ofstream outFile("data/isco_lNP_a_E_v2.dat");


    FILE *foutput;
    FILE *foutput_coord;
    
    /* ----- Set free parameters ----- */
    //We use unit c=G=M = 1 and set the mass of the BH scaled to be 1
    spin = atof(argv[1]);
    iobs_deg = atof(argv[2]);
    beta = atof(argv[3]);
    lNP = atof(argv[4]);
    rstep = atof(argv[5]);
    pstep = atof(argv[6]);
    alpha = atof(argv[7]);
    //beta = 1.0;
    //lNP = 1.0e-30;
    
    iobs = Pi/180*iobs_deg;
    //const size_t Npts = 5001; // Number of points for horizon calculation
    double dtheta = -Pi/2/(HorizonN-1); // Step size for theta in horizon calculation
    global_horizon = find_horizon<HorizonN>(dtheta);

    // Read rHorizonLimit.dat into arrays
    std::vector<double> spin_vals, rHorizon_vals;
    std::ifstream rHorizon_file("data/rHorizonLimit.dat");
    if (!rHorizon_file) {
        std::cerr << "Error: Could not open data/rHorizonLimit.dat" << std::endl;
        exit(1);
    }
    std::string line;
    while (std::getline(rHorizon_file, line)) {
        std::istringstream iss(line);
        double spin_val, rHorizon_val;
        if (iss >> spin_val >> rHorizon_val) {
            spin_vals.push_back(spin_val);
            rHorizon_vals.push_back(rHorizon_val);
        }
    }
    rHorizon_file.close();

    /*
    std::ofstream horizon_file("data/horizon_points.dat");
    std::ofstream horizon_file2("data/horizon_points_2.dat");
    
    if (horizon_file.is_open()) {
        for (size_t i = 0; i < HorizonN; ++i) {
            horizon_file << std::setprecision(max_prec)<<global_horizon.theta[i] << " " << global_horizon.H[i]  << "\n";
        }
        horizon_file.close();
    } else {
        std::cerr << "Error opening horizon_points.dat for writing." << std::endl;
    }   
    
    if (horizon_file2.is_open()) {
    for(double theta1 = 0.0; theta1 <= Pi/2; theta1 += 0.01) {
        //three point interpolation 
        index = (int) floor(-(Pi/2.0 - theta1) / dtheta); 
        if (index <= 0) {
            index = 1;
            //rH = horizon.H[0]+ (horizon.H[1] - horizon.H[0]) * (theta1 - horizon.theta[0]) / (horizon.theta[1] - horizon.theta[0]);
        }
        else if (index >= HorizonN-1) {
            index = HorizonN - 2;
        }
        rH = three_point_interpolation(theta1,global_horizon.theta[index-1], global_horizon.theta[index], global_horizon.theta[index+1],
                                       global_horizon.H[index-1], global_horizon.H[index], global_horizon.H[index+1]);
        // Can you write this rH to a separate file?
        horizon_file2 << std::setprecision(max_prec) << index << " "<<theta1 << " " << rH << "\n"; 
    }
        horizon_file2.close();
    } else {
        std::cerr << "Error opening horizon_points.dat for writing." << std::endl;
    }   
    */
    
    

    /* ----- Set inner and outer radius of the disk ----- */
    // find the isco in the new coordinate system 
    //find_isco(15.0, isco);
    //cout<<"isco v1 "<<isco<<endl;
    //isco = find_isco_NC(2, 40, 1e-5, 1e-8, 1000, 1e-7);
    //cout<<"isco NC "<<isco<<endl;
     // Initial step size for the second derivative
    auto findIsco = [&](double x) {
        return veffddr(x, h);
    };
    constexpr double tol = 1e-10;
    constexpr int max_iter = 100;

    /*
    std::ofstream isco_file("data/isco_file.dat");
    if (isco_file.is_open()) {
        for(double r = global_horizon.H[0]+0.1; r <= 100; r += 0.01) {
            isco_file << std::setprecision(max_prec) <<r<<"  "<< findIsco(r) << "\n";
        }   
        isco_file.close();
    } else {
        std::cerr << "Error opening isco_file.dat for writing." << std::endl;
    }
    */

    
    
    for(r = 10 ; r >= global_horizon.H[0]+0.1; r -= 0.001) {
        if(findIsco(r)*findIsco(r-0.001) < 0) {
            isco_initial_guess = r;
            break;
        }
    }   
    

    HighPrecisionRootFinder root_finder_isco(findIsco, tol, max_iter);
    isco = root_finder_isco.newton_raphson(isco_initial_guess);
    
   
    //Find ISCO for different spin and lNP
    /*
    outFile << std::setprecision(max_prec);
    for(spin=0.6;spin<=0.99;spin+=0.01){
        horizon_limit = get_rHorizon_for_spin(spin, spin_vals, rHorizon_vals);
        for(lNP=1e-40;lNP<horizon_limit-0.0001;lNP+=0.01){
            global_horizon = find_horizon<HorizonN>(dtheta);
            for(r = 10 ; r >= global_horizon.H[0]+0.1; r -= 0.001) {
                if(findIsco(r)*findIsco(r-0.001) < 0) {
                    isco_initial_guess = r;
                    break;
                }
            }
            HighPrecisionRootFinder root_finder_isco(findIsco, tol, max_iter);
            isco = root_finder_isco.newton_raphson(isco_initial_guess); 
                      
            outFile << lNP << " " << spin  <<" "<<isco<<'\n';
        }
    }
    
    outFile.close();
    */
    xin = isco;
    xout = 500;

    double rstep2 = (rstep - 1)/rstep;

    /* ----- Set model for the spectral line ----- */
	
	double E_line = 6.4;   /* energy rest of the line in keV */
	double N_0    = 1.0;   /* normalization */
	// pstep  = 2*Pi/720;
	
	E_obs[0] = 0.0125000002;     /* minimum photon energy detected by the observer; in keV */
	N_obs[0] = 0;
	for (i = 1; i <= imax - 1; i++) {
		E_obs[i] = E_obs[i - 1] + 0.025;
		N_obs[i] = 0;
	}
	
	/*Iron line output file*/
	//sprintf(filename_o,"iron_a%.03f.epsilon_r%.02f.epsilon_t%.02f.i%.02f.dat",spin,epsi3,iobs_deg);
    // sprintf(filename_o,"ironline_data/iron_a%.05Le.i%.02Le.e_%.02Le.a13_%.02Le.a22_%.02Le.a52_%.02Le.dat",spin,iobs_deg,epsi3,a13,a22,a52);
    snprintf(filename_o, sizeof(filename_o), "ironline_data/iron_a_%.05Lf_i_%.05Lf_beta_%.05Lf_lNP_%.05Lf.dat",spin,iobs_deg,beta,lNP);
    
    /* ----- Set computational parameters ----- */
    robs_i = 1;
    robs_f = 500;
    //robs_f = 3;

   snprintf(filename_o2, sizeof(filename_o2), "data/photons_data_a%.05Lf_i_%.05Lf_beta_%.05Lf_lNP_%.05Lf.dat",
             spin, iobs_deg, beta, lNP);
    
    foutput_coord = fopen(filename_o2,"w");
    if (foutput_coord == NULL) {
        perror("Error opening output file");
        exit(1);
    }
    
    /* ----- Photon tracing loop ----- */
    for (robs = robs_i; robs < robs_f; robs = robs*rstep) {
        for (pobs = 0; pobs < 2*Pi - 0.5*pstep; pobs = pobs + pstep) {
        //for (pobs = 0; pobs < Pi/3; pobs = pobs + pstep) {
            xobs = robs*cos(pobs);
            yobs = robs*sin(pobs);
            
            
            auto traceOut = raytracingNC(xobs, yobs, iobs, xin, xout);
            //fprintf(foutput_coord, "%d %d %Lf %Lf %Lf %Lf %Lf\n", 
            //            photon_index, stop_integration_condition, xobs, yobs, traced[0], traced[3], traced[1]);
            //photon_index++;
            stop_integration_condition = traceOut.stop_integration;
            
            if (stop_integration_condition == 1) { 
                //printf("Writing photon %d to %s\n", photon_index, filename_o2);
                fprintf(foutput_coord, "%d %d %Lf %Lf %Lf %Lf %Lf\n", 
                        photon_index, stop_integration_condition, xobs, yobs, traceOut.varsOut[0], traceOut.varsOut[1], traceOut.varsOut[2]);
                //fprintf(foutput_coord, "%d %Lf  %Lf %Lf %Lf %Lf\n", 
                //        photon_index, xobs, yobs, traced[0], traced[3], traced[1]);
                photon_index++;
                gfactor = traceOut.varsOut[1];
				pp = gfactor*E_line;
				/* --- integration - part 1 --- */
				
				for (i = 0; i <= imax - 2; i++) {
					if (E_obs[i] < pp && E_obs[i + 1] > pp) {
						
						qq = gfactor*gfactor*gfactor*gfactor;
						qq = qq*pow(traceOut.varsOut[0],alpha);
						
						fphi[i] = fphi[i] + qq;
						
					}
				}		
            } 
        }

        /* --- integration - part 2 --- */
		
		for (i = 0; i <= imax - 1; i++) {			
			fr = robs*robs*fphi[i]*rstep2;		
			N_obs[i] = N_obs[i] + fr;			
		}			
    }


    /* --- print spectrum --- */
	
	foutput = fopen(filename_o,"w");	
    if (foutput == NULL) {
        perror("Error opening output iron file");
        exit(1);
    }
	long double N_tot  = 0.0;
	
	for (i = 0; i <= imax - 1; i++) {
		N_obs[i] = N_0*N_obs[i]/E_obs[i];
		N_tot = N_tot + N_obs[i];
	}

	for (i = 0; i <= imax - 1; i++) {
		fprintf(foutput,"%Lf %.10Lf\n",E_obs[i],(N_obs[i]/N_tot));
	}
	
	fclose(foutput);
    
    fclose(foutput_coord);
    return 0;
}



