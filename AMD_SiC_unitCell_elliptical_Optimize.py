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
# a = 2.801277586646125e-7
a = 2.737678125000000e-07
#hole diameter prefactor 1
d1 = 0.67
#hole diameter prefactor 2
# d2 = 0.77
d2 = 1.2
#beam width prefactor
w = 1.75
#taper prefactor (for the defect region)
t = 0.8
#taper prefactor (for the waveguide region)
t_wvg = 0.875
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 15e-6
# The target resonance frequency, in Hz
# 916nm = 327.3e12
target_frequency = 327.3e12
#the prefactor associated with the weaker mirror region
prefactor_mirror_R = 1
#the refractive index associated with the material 
n_f = 2.6
#the minimum lattice constant in the waveguide region 
amin_wvg = t_wvg*a

# # do a low resolution sweep over the desired parameter range (LOOPING CODE) ############
# a_list = np.linspace(2.7e-07,2.9e-07,10)
# dx_list = np.linspace(0.2,0.6,10)
# dy_list = np.linspace(0.5,w,10)
# for a in a_list:
#     for dx in dx_list:
#         hx = dx*a/2
#         for dy in dx_list:
#             hx = dy*a/2
#             w0 = w*a
#             p0 = [a,hx,hy,w0]
#             unitCellOptimization_SiC_elliptical(p0)

# optimizing for the mirror unit cells (SWEEPING CODE) ###################
a_0 = 2.70e-7
hx_0 = 8.10e-08
hy_0 = 7.00e-8
w0 = w*a_0
p0 = [a_0,hx_0,hy_0,w0]
bnd = ((2.5e-07,3.0e-07),(None,a_0),(None,a_0),(None,None))
popt = scipy.optimize.minimize(unitCellOptimization_SiC_elliptical,p0,method='Nelder-Mead')

# # optimizing for the waveguide unit cells
# p0 = [amin_wvg]
# bnd = ((0.5*a,a))
# popt = scipy.optimize.minimize(unitCellOptimization_SiC_waveguide_elliptical,p0,method='Nelder-Mead')


# testing the algorithm
# unitCellOptimization_SiC_elliptical(p0)