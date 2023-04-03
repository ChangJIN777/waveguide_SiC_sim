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
a = 2.775047787137548e-07
#hole diameter in the x direction
hx = 7.076170823568784e-08
hx /= 2
#hole diameter in the y direction
hy = 1.730259002115936e-07
hy /= 2
#beam width 
w0 = 4.699560661981287e-7
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

# # sweep through hy and record the bandgap and midband value 
# hy_list = np.linspace(hy/2,hy,20)
# for hy_sweep in hy_list:
#     hy_nm = hy_sweep*1e9
#     print("hy: %f nm " %(hy_nm))
#     p0 = [a,hx,hy_sweep,w0,h0]
#     unitCellOptimization_SiC_elliptical(p0)
# debugging
# p0 = [a,hx,hy,w0,h0]
# unitCellOptimization_SiC_elliptical(p0)

# # do a low resolution sweep over the desired parameter range (LOOPING CODE/Mirror region) ############
# # a_list = np.linspace(2.7e-07,2.9e-07,10)
# # a_list = np.linspace(2.50e-7,3.00e-7,10)
# hx_list = np.linspace(5e-08,1e-07,10)
# hy_list = np.linspace(5e-08,2e-07,20)
# for hx in hx_list:
#     for hy in hy_list:
#         p0 = [a,hx,hy,w0,h0]
#         hx_nm = hx*1e9
#         hy_nm = hy*1e9
#         a_nm = a*1e9
#         print("hx: %f nm hy: %f nm a: %f nm" %(hx_nm, hy_nm, a_nm))
#         unitCellOptimization_SiC_elliptical(p0)
            
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

# Sweep the lattice constants
a_list = np.linspace(2.00e-7,3.00e-7,20)
for a in a_list:
    p0 = [a,hx,hy,w0,h0]
    a_nm = a*1e09
    print("the lattice constant: %f nm" %(a_nm))
    unitCellOptimization_SiC_elliptical(p0)