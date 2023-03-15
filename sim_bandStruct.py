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

# simulating the unit cell band structures
def bandStructSim(a,d,w,n_f,fmin,fmax,f_grating,h0=250e-9):
    """this function simulate the band structure of a unit cell with parameters given by the inputs
        a: the lattice constant 
        d: hole diameter prefactor 
        w: beam width prefactor
        h0: beam height 
        n_f: the refractive index of the dielectric material 
    """
    #beam width
    w0 = w*a
    # Radius of the air holes in the cells
    r0 = (d*a)/2
    FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
    engine = LumericalEngine(mesh_accuracy=4, hide=True, lumerical_path=FDTDloc, working_path="./fsps", save_fsp=False)

    # building the unit cell 
    # the sim material is set to be SiC with refractive index = 2.6 
    cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    simulate_unit_cell = UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)
    
    simulate_unit_cell.simulate("bandstructure", freqs=(fmin, fmax, f_grating))
    
    return 

#lattice constant
a = 2.80e-7
#hole diameter prefactor 
d = 0.64
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
amin = d*a
h0 = 250e-9
a_taper = cubic_tapering(a,amin,taperNum=TN) 
a_taper = a_taper[::-1]
band_structure(a_taper[-1],d,w,h0,n_f,engine=engine)
# bandStructSim(a,d,w,n_f,fmin,fmax,f_grating)
