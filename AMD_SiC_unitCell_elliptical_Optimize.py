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
w = 1.69
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

# optimizing for the mirror unit cells
p0 = [a,d1,d2]
bnd = ((2.0e-07,3.0e-07),(0.1,0.8),(0.67,w))
# popt = scipy.optimize.minimize(unitCellOptimization_SiC_elliptical,p0,method='Nelder-Mead')

# # optimizing for the waveguide unit cells
# p0 = [amin_wvg]
# bnd = ((0.5*a,a))
# popt = scipy.optimize.minimize(unitCellOptimization_SiC_waveguide_elliptical,p0,method='Nelder-Mead')


# testing the algorithm
unitCellOptimization_SiC_elliptical(p0)