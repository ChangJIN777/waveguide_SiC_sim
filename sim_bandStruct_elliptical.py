from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os
import csv
import struct
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import style
from waveguideSolver_funcs import *

#lattice constant
a = 2.708900299422566e-7
#the tapering prefactor 
t = 0.6
#hole diameter prefactor 1
hx = 6.936330787756172e-8
#hole diameter prefactor 2
hy = 1.639998046215848e-7
#beam width prefactor
w = 1.75 
w0 = 5.102578564047907e-7
#refractive index of the dielectric material 
n_f = 2.6
#the number of taper cells
TN = 8
#define the frequency range 
target_frequency = 327.3e12 #Hz 
freq_span = 50e12 #Hz
fmin = target_frequency-freq_span
fmax = target_frequency+freq_span
f_grating = 500
amin = t*a
h0 = 250e-9
a_taper = cubic_tapering(a,amin,taperNum=TN) 
a_taper = a_taper[::-1]
# simulate the band gap associated with the unit cell 
# sim_bandGap_elliptical(a,d1,d2)

# # simulate the band structure associated with the unit cell (mirror region)
# band_structure_elliptical(a,hx,hy,w0)
# sim_bandGap_elliptical(a,hx,hy,w0)

# # simulate the band structure associated with the unit cell (defect region)
t_list = np.linspace(0,1,20)
for t in t_list:
    a_def = a*t
    # band_structure_elliptical(a_def,hx,hy,w0)
    print("t = %f" % (t))
    sim_bandGap_elliptical(a_def,hx,hy,w0)
