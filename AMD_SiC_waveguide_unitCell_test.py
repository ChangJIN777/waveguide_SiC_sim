"""
This script simulates the band structures of the waveguide unit cells to make sure that it does not pick up higher order modes 
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
a = 263.3e-09
#hole diameter prefactor 
d = 0.64
#beam width prefactor
w = 1.69
#taper prefactor (for the defect region)
t = 0.8
#taper prefactor (for the waveguide region)
t_wvg = 0.75
#beam height (set by epi-layer thickness)
h0 = 500e-9
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
#hole diameter prefactor (in the waveguide region)
d_wvg = 0.64

band_structure(amin_wvg,d_wvg,w,h0,n_f)
