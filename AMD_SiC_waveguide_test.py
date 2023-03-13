"""
This example creates a waveguide, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os

from waveguideSolver_funcs import *

# Define geometry paramaters 
#waveguide taper cell number
WN = 9
#taper cell number (left mirror region)
TN = 8
#mirror cell number (left region) 
MN_L = 20
#mirror cell number (right region)
MN_R = 3
#lattice constant
a = 3.348909692268754e-07
#hole diameter prefactor 
d = 0.64
#beam width prefactor
w = 1.75
#taper prefactor (for the defect region)
t = 0.84
#taper prefactor (for the waveguide region)
t_wvg = 0.5
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 15e-6
# The target resonance frequency, in Hz
# 1280nm = 234.212857812500e12
# 916nm = 327.3e12
target_frequency = 327.3e12
#the prefactor associated with the weaker mirror region
prefactor_mirror_R = 1
#the refractive index associated with the material 
n_f = 2.6
# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations 
FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
engine = LumericalEngine(mesh_accuracy=5, hide=False, lumerical_path=FDTDloc, working_path="./fsps")
# engine = LumericalEngine(mesh_accuracy=5, hide=True, lumerical_path=FDTDloc, save_fsp=False)
#the minimum lattice constant in the tapering region
amin = a*t
    
#build the left mirror cell region 
mirror_cells_left = buildMirrorRegion(a,d,w,h0,n_f,MN_L,engine)

#build the right mirror cell region 
a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
mirror_cells_right = buildMirrorRegion(a_R,d,w,h0,n_f,MN_R,engine)

#building cubic tapered cell region
taper_cells = buildTaperRegion(a,a_R,amin,d,w,h0,n_f,TN,engine)

#set the center of the device
centerCell = MN_L+TN-1

#adding one sided cubic tapered waveguide region to the cavity
waveguide_cells_R = buildWaveguideRegion_right(a_R,d,w,t_wvg,h0,WN,n_f,engine)

cavity = Cavity1D(
unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right + waveguide_cells_R ,
structures=[ BoxStructure(Vec3(0), Vec3(l, w*a, h0), DielectricMaterial(n_f, order=2, color="red")) ],
center_cell=centerCell,
center_shift=0,
engine=engine
)

# By setting the save path here, the cavity will save itself after each simulation to this file
cavity.save("cavity_testing.obj")

#define mesh size (use 12nm for accuracy, currently set to 50nm)
man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 20e-9, dy=None, dz=None)

# simulating the resonance and the Q =================================================
# r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
#                     analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,4,10))
r1 = cavity.simulate("resonance", target_freq=target_frequency, 
                    analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,8,10))

# Print the reults and plot the electric field profiles
print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
r1["freq"], r1["vmode"],
1/(1/r1["qxmin"] + 1/r1["qxmax"]),
1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
))
# r1["xyprofile"].show()
# r1["yzprofile"].show()

cavity = Cavity1D(load_path="cavity_testing.obj",engine=engine)
Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
Vmode = r1["vmode"]
F = r1["freq"]

Q = 1/((1/Qsc) + (1/Qwvg))
P = (Q*Qsc) / (Vmode*Vmode)
print("Q: %f, P: %f" % ( Q, P))

r1 = cavity.get_results("resonance")[-1]
print(r1['res']["xyprofile"].max_loc())
print(r1['res']["yzprofile"].max_loc())
# r1["sess_res"].show()
# # ======================================================================================

# # evaluate the quasipotential
r2 = cavity.simulate("quasipotential", target_freq=target_frequency)
r2.show()

# file = open("OptimizeList.txt","a") 
# file.write("\n" + str(a) + " " + str(Q) + " " + str(Vmode)+ " " + str(F) + "\n") 
# file.close()

