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
a = 2.838218812324136e-7
#hole diameter in the x direction
hx = 7.160169206987993e-08
#hole diameter in the y direction
hy = 1.652696864561149e-07
#beam width 
w0 = 5.005507792174242e-07
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 10e-6
# The target resonance frequency, in Hz
# 916nm = 327.3e12
target_frequency = 327.3e12
#the prefactor associated with the weaker mirror region
prefactor_mirror_R = 1
#the refractive index associated with the material 
n_f = 2.6
#the minimum lattice constant in the waveguide region 
amin_wvg = t_wvg*a

# sweep through hy and record the bandgap and midband value 
hy_list = np.linspace(hy/2,hy,20)
for hy_sweep in hy_list:
    hy_nm = hy_sweep*1e9
    print("hy: %f nm =============" %(hy_nm))
    p0 = [a,hx,hy_sweep,w0,h0]
    unitCellOptimization_SiC_elliptical(p0)

# # do a low resolution sweep over the desired parameter range (LOOPING CODE/Mirror region) ############
# # a_list = np.linspace(2.7e-07,2.9e-07,10)
# a_list = np.linspace(2.50e-7,3.00e-7,10)
# dx_list = np.linspace(0.5,0.9,5)
# dy_list = np.linspace(0.5,1.75,5)
# for a in a_list:
#     for dx in dx_list:
#         hx = dx*a/2
#         for dy in dy_list:
#             hy = dy*a/2
#             w0 = w*a
#             p0 = [a,hx,hy,w0]
#             hx_nm = hx*1e9
#             hy_nm = hy*1e9
#             a_nm = a*1e9
#             print("a: %f nm hx: %f nm hy: %f nm" %(a_nm, hx_nm, hy_nm))
#             unitCellOptimization_SiC_elliptical(p0)
            
# # optimizing for the mirror unit cells (SWEEPING CODE) ###################
# a_0 = a
# hx_0 = hx*0.75
# hy_0 = hy*0.75
# w0 = w*a_0
# p0 = [a_0,hx_0,hy_0,w0,h0]
# bnd = ((2.5e-07,3.0e-07),(None,None),(None,None),(None,None),(None,None))
# popt = scipy.optimize.minimize(unitCellOptimization_SiC_elliptical,p0,method='Nelder-Mead')

# # optimizing for the waveguide unit cells
# p0 = [amin_wvg]
# bnd = ((0.5*a,a))
# popt = scipy.optimize.minimize(unitCellOptimization_SiC_waveguide_elliptical,p0,method='Nelder-Mead')

# # do a low resolution sweep over the desired parameter range (LOOPING CODE/defect region) ############
# a_list = np.linspace(2.7e-07,2.9e-07,5)
# hx_list = np.linspace(5e-8,1e-7,5)
# hy_list = np.linspace(1e-7,2e-7,5)
# w0_list = np.linspace(4e-7,8e-7,5)
# for a in a_list:
#     for hx in hx_list:
#         for hy in hy_list:
#             for w0 in w0_list:
#                 w0 = w*a
#                 a_def = a*t
#                 p0 = [a_def,hx,hy,w0]
#                 print("tapering prefactor: %f" %(t))
#                 unitCellOptimization_SiC_elliptical(p0)

# testing the algorithm
# unitCellOptimization_SiC_elliptical(p0)

# # Sweep the beam width 
# w0_list = np.linspace(4e-7,8e-7,10)
# for w0 in w0_list:
#     p0 = [a,hx,hy,w0,h0]
#     print("beam width: %f" %(w0))
#     unitCellOptimization_SiC_elliptical(p0)