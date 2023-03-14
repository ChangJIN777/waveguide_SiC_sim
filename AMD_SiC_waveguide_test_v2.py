"""
This code creates only the waveguide portion of the cavity and simulate its band structure to get a good initial guess for the optimization code. 
"""

from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os
import scipy.constants
import datetime as datetime
from waveguideSolver_funcs import *


# Define geometry paramaters 
#waveguide taper cell number
WN = 6 
#lattice constant
a = 280e-9
#hole diameter prefactor 
d = 0.64
#beam width prefactor
w = 1.69
#taper prefactor (for the defect region)
t = 0.8
#taper prefactor (for the waveguide region)
t_wvg = 1
#taper prefactor for the hole size of the waveguide region 
d_wvg_min = 0.2
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 15e-6
# The target resonance frequency, in Hz
# 916nm = 327.3 THz 
target_frequency = 327.3e12
target_wavelength = 9.16e-07


# Define geometry dependencies
#beam width
w0 = w*a
# Radius of the air holes in the cells
r0 = (d*a)/2
#min tapered lattice (waveguide)
amin_wvg = t_wvg*a 
#min taper hole diameter (waveguide)
rmin_wvg = d_wvg_min*a/2
#radius taper rate (for the waveguide region)
r_wvg_tr = (r0-rmin_wvg)/WN
#refractive index of the material we are trying to simulate (SiC = 2.6)
n_f = 2.6

# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
# FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
FDTDloc='C:/Program Files/Lumerical/v221/' # for running on the local desktop
engine = LumericalEngine(mesh_accuracy=4, hide=True, lumerical_path=FDTDloc, save_fsp=False)

# build the tapered waveguide unit cells 
r_wvg = np.linspace(rmin_wvg,r0,WN)

# # simulate the band gap 
# for i in r_wvg:
#     i_nm = i*1e9
#     d = i/a
#     print("radius: %f nm" % (i_nm))
#     sim_bandGap(a,d,w,h0,n_f,engine)

#simulate the band structure 
d_BS = d_wvg_min #choose which cell to simulate
target_frequency = 327.3e12 #Hz 
freq_span = 50e12 #Hz
fmin = target_frequency-freq_span
fmax = target_frequency+freq_span
f_grating = 5000
kmin = 0.2
kmax = 0.5
knum = 5
band_structure(a,d_BS,w,h0,n_f,engine)