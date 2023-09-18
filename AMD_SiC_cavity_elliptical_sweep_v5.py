"""
This example creates a sweep the parameters of the unit cells for the waveguide portion of the cavity to find a good initial guess
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

#lattice constant
a = 263.3e-09
#hole diameter in the x direction 
hx = 72.40e-09 
#hole diameter in the y direction 
hy = 127.5e-09
#beam width prefactor
w0 = 468.6e-09
#taper prefactor (for the defect region)
t = 0.818
#taper prefactor (for the waveguide region)
t_wvg = 0.852
#beam height (set by epi-layer thickness)
h0 = 500e-9
# cavity beam length
l = 20e-6
# The target resonance frequency, in Hz
# 916nm = 327.3e12
target_frequency = 327.3e12
target_wavelength = 9.16e-07
#the prefactor associated with the weaker mirror region
prefactor_mirror_R = 0.965
#the refractive index associated with the material 
n_f = 2.6
#the minimum lattice constant in the waveguide region 
amin_wvg = t_wvg*a
# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations 
# FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
FDTDloc="C:/Program Files/Lumerical/v221/" # note: this is specified to be run on Feynman 
# engine = LumericalEngine(mesh_accuracy=5, hide=False, lumerical_path=FDTDloc, working_path="./fsps")
engine = LumericalEngine(mesh_accuracy=4, hide=True, lumerical_path=FDTDloc, save_fsp=False)
#the minimum lattice constant in the tapering region
amin = a*t
#the minimum radius prefactor we are tapering to 
d_min = 0.437
#the left mirror cell number 
MN_L = 10 
#the right mirror cell number 
MN_R = 3
#the number of taper unit cells 
TN = 6
TN_L = 6 
TN_R = 6
#set the center of the device
centerCell = MN_L+TN_L-1
#the number of cells in the waveguide region
WN = 3
#filename we are saving the data under 
file_name = "OptimizeListFull_Overcoupled_cavity_500nm_testRun_2.csv"

# debugging 
#hx_min = d_min*hx
#hy_min = d_min*hy
#sim_ellipticalCavity_v2(a,hx,hx_min,hy,hy_min,t,t_wvg,target_frequency,file_name,MN_L,MN_R,TN,w0=w0,engine=engine)

# sweeping the unit cell parameters


# # sweeping the taper prefactor 
# t_list = np.linspace(0.4,1,20)
# for t in t_list:
#     param = [t]
#     sweep_tapering_elliptical_cavity(param)

# sweeping the beam width  
# w_list = np.linspace(4.00e-7,8.00e-7,20)
# for w in w_list:
#     param = [w]
#     w_nm = w*1e9
#     print("the beam width: %f nm" %(w_nm))
#     sweep_beamWidth_ellipticalCavity_v2(param)

# optimize the beam width  
# p0 = [w0]
# bnd = [(None,1e-6)]
# popt = scipy.optimize.minimize(sweep_beamWidth_ellipticalCavity_v2,p0,method='Nelder-Mead')
# sweep_beamWidth_ellipticalCavity_v2(param)


# # optimization algorithm (only the beam width)
# p0 = [w0]
# bnd = [(None,None)]
# popt = scipy.optimize.minimize(sweep_beamWidth_ellipticalCavity_v2,p0,method='Nelder-Mead')

# # sweeping the cell numbers 
# TN_list = [2,3,4,5,6,7,8]
# MN_L_list = [1,2,3,4,5,6,7,8,9,10]
# for TN in TN_list:
#     for MN_L in MN_L_list:
#         param = [MN_L,TN]
#         sweep_cellNum_ellipticalCavity(param)