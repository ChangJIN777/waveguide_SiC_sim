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
a = 2.80e-7
#the tapering prefactor 
t = 0.8
#hole diameter prefactor 1
d1 = 0.64
#hole diameter prefactor 2
d2 = 0.8
#beam width prefactor
w = 1.69
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
sim_bandGap_elliptical(a,d1,d2)
# band_structure(a_taper[-1],d,w,h0,n_f,engine=engine)
# bandStructSim(a,d,w,n_f,fmin,fmax,f_grating)
