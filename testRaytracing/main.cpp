#include "def.h"

int main(int argc, char *argv[])
{
    long double alpha;
    long double iobs;
    long double robs, pobs;
    long double robs_i, robs_f, rstep, pstep;
    long double xin, xout;
    long double isco;
    long double traced[5];    
    int stop_integration_condition = 0;
    int photon_index = 0;
    char filename_o2[128];
    
    FILE *foutput_coord;
    
    /* ----- Set free parameters ----- */
    spin = atof(argv[1]);
    iobs_deg = atof(argv[2]);
    a13 = atof(argv[3]);
    a22 = atof(argv[4]);
    a52 = atof(argv[5]);
    epsi3 = atof(argv[6]);
    alpha = atof(argv[7]);
    rstep = atof(argv[8]);
    pstep = atof(argv[9]);
    
    iobs = Pi/180*iobs_deg;
    
    /* ----- Set inner and outer radius of the disk ----- */
    find_isco(15.0, isco);
    xin = isco;
    xout = 500;
    
    /* ----- Set computational parameters ----- */
    robs_i = 1;
    //robs_f = 400;
    robs_f = 200;
    snprintf(filename_o2, sizeof(filename_o2), "data/photons_data_a%.05Lf_i_%.05Lf_e_%.05Lf_a13_%.05Lf_a22_%.05Lf_a52_%.05Lf.dat",
             spin, iobs_deg, epsi3, a13, a22, a52);
    
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
            
            raytrace(xobs, yobs, iobs, xin, xout, traced, stop_integration_condition);
            //fprintf(foutput_coord, "%d %d %Lf %Lf %Lf %Lf %Lf\n", 
            //            photon_index, stop_integration_condition, xobs, yobs, traced[0], traced[3], traced[1]);
            //photon_index++;
            
            if (stop_integration_condition == 1) { 
                //printf("Writing photon %d to %s\n", photon_index, filename_o2);
                fprintf(foutput_coord, "%d %d %Lf %Lf %Lf %Lf %Lf\n", 
                        photon_index, stop_integration_condition, xobs, yobs, traced[0], traced[3], traced[1]);
                //fprintf(foutput_coord, "%d %Lf  %Lf %Lf %Lf %Lf\n", 
                //        photon_index, xobs, yobs, traced[0], traced[3], traced[1]);
                photon_index++;
            } 
        }	
    }
    
    fclose(foutput_coord);
    return 0;
}
