"""
This example creates a waveguide, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os
from datetime import datetime

from waveguideSolver_funcs import *

# Define geometry paramaters 
#waveguide taper cell number
WN = 9
#taper cell number (left mirror region)
TN = 8
#mirror cell number (left region) 
MN_L = 20
#mirror cell number (right region)
MN_R = 3
#lattice constant
a = 2.801277586646125e-7
#hole diameter prefactor 
d = 0.642241102375429
#beam width prefactor
w = 1.692961161046466
#taper prefactor (for the defect region)
t = 0.84
#taper prefactor (for the waveguide region)
t_wvg = 0.9
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 15e-6
# The target resonance frequency, in Hz
# 1280nm = 234.212857812500e12
# 916nm = 327.3e12
target_frequency = 327.3e12
#the prefactor associated with the weaker mirror region
prefactor_mirror_R = 0.9
#the refractive index associated with the material 
n_f = 2.6
# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations 
FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
# engine = LumericalEngine(mesh_accuracy=5, hide=False, lumerical_path=FDTDloc, working_path="./fsps")
engine = LumericalEngine(mesh_accuracy=5, hide=True, lumerical_path=FDTDloc, save_fsp=False)
amin = a*t
a_R = a*prefactor_mirror_R

#testing the tapered lattice constants 
a_taper_R = cubic_tapering(a,amin,TN)
a_taper_L = a_taper_R[::-1]
a_tapering = np.concatenate((a_taper_L, a_taper_R), axis=None)

#simulate the band gap 
for i in a_tapering:
    i_nm = i*1e9
    print("a: %f nm" % (i_nm))
    sim_bandGap(i,d,w,h0,n_f,engine)
# #simulate the band structure 
# band_structure(a,d,w,h0,n_f,engine)

# #Optimize the unit cell 
# p0=[a,d,w]
# bnd = ((2.50e-07,3.00e-07),(0,1),(1,None))
# # unitCellOptimization_SiC(p0) #debugging
# popt = scipy.optimize.minimize(unitCellOptimization_SiC,p0,method='Nelder-Mead',bounds=bnd)
