import os
import sys
import numpy as np

# ----------------------------------------------------------------------------------------------------
# INPUT PARAMETERS
# ----------------------------------------------------------------------------------------------------

#path_to_xillver_file = "xillver-a-Ec5.fits"
gamma = 1.8
afe = 1
logxi = 1.
ecut = 300


spin = 0.8
incl = 30
beta = 2.0
lNP = 0.32
alpha = -3

Nr = 2000
Nph = 500
Rmax = 500
rstep  = np.e**((np.log(Rmax)-np.log(np.cos(np.pi*incl/180)))/Nr)
pstep  = 2*np.pi/Nph

# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------

rt_command = "./testRayTracing %f %f %f %f %f %f %f"%(spin, incl, beta, lNP, rstep, pstep, alpha)
# os.system("cd raytracer/")
print("Start of ray-tracing part ...")
os.system(rt_command)
print("End of ray-tracing part")