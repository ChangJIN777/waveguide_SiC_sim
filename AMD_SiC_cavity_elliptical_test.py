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
a = 2.811543831305853e-7
#hole diameter in the x direction
hx = 7.002993766800964e-8
#hole diameter in the y direction
hy = 1.680746149495162e-7
#beam width prefactor
# w = 1.75 
w = 1.75
w0 = 4.945210709892193e-7
# w0 = 4.699560661981287e-7
#taper prefactor (for the defect region)
t = 0.8153
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 10e-6
# The target resonance frequency, in Hz
# 916nm = 327.3e12
target_frequency = 327.3e12
#the prefactor associated with the weaker mirror region
prefactor_mirror_R = 1
#taper prefactor (for the waveguide region)
t_wvg = 0.75
#the refractive index associated with the material 
n_f = 2.6
#the minimum lattice constant in the waveguide region 
amin_wvg = t_wvg*a
# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations 
FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
engine = LumericalEngine(mesh_accuracy=5, hide=False, lumerical_path=FDTDloc, working_path="./fsps")
# engine = LumericalEngine(mesh_accuracy=5, hide=True, lumerical_path=FDTDloc, save_fsp=False)
#the minimum lattice constant in the tapering region
amin = a*t
#the minimum radius prefactor we are tapering to 
d_min = 0.3
#the left mirror cell number 
MN_L = 10
#the right mirror cell number 
MN_R = 3
#the number of taper unit cells 
TN = 6
#set the center of the device
centerCell = MN_L+TN-1 

#build the left mirror cell region 
mirror_cells_left = buildMirrorRegion_elliptical(a,hx,hy,MN_L,w0,h0,n_f,engine)

#build the right mirror cell region 
a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
mirror_cells_right = buildMirrorRegion_elliptical(a_R,hx,hy,MN_R,w0,h0,n_f,engine)

#building cubic tapered cell region
taper_cells = buildTaperRegion_elliptical(a,a_R,amin,hx,hy,TN,w0,h0,n_f,engine)

# ####################################### cavity without the waveguide region (symmetric) ###############################
# cavity = Cavity1D(
# unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right,
# structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
# engine=engine
# )
####################################### cavity without the waveguide region (asymmetric) ###############################
cavity = Cavity1D(
unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right,
structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
center_cell=centerCell,
center_shift=0,
engine=engine
)

# By setting the save path here, the cavity will save itself after each simulation to this file
cavity.save("cavity.obj")

#define mesh size (use 12nm for accuracy, currently set to 12nm)
# man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 12e-9, dy=None, dz=None)
# man_mesh = MeshRegion(BBox(Vec3(0),Vec3(12e-6,0.7e-6,0.4e-6)), 20e-9, dy=None, dz=None)
man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,3e-6,3e-6)), 15e-9, dy=None, dz=None)

# simulating the resonance and the Q #########################################################
r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
                    analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,4,8))
# r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=60e-15, 
#                     analyze_time=600e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,4,8))
# try simulating without a target frequency 
# r1 = cavity.simulate("resonance", source_pulselength=200e-15, 
#                     analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,8,8))


# Print the reults and plot the electric field profiles
print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
    r1["freq"], r1["vmode"],
    1/(1/r1["qxmin"] + 1/r1["qxmax"]),
    1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
))
r1["xyprofile"].show()
r1["yzprofile"].show()

cavity = Cavity1D(load_path="cavity.obj",engine=engine)
Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
Vmode = r1["vmode"]
F = r1["freq"]

# for debugging the Q 
Qx1 = r1["qxmin"]
Qx2 = r1["qxmax"]
Qy = 1 / (2 / r1["qymax"])
Qz = 1 / (1 / r1["qzmin"] + 1 / r1["qzmax"])
print("Qx1: %f, Qx2: %f, Qy: %f, Qz: %f" % (
    Qx1, Qx2, Qy, Qz
))

Q = 1/((1/Qsc) + (1/Qwvg))
P = (Q*Qsc) / (Vmode*Vmode)
print("Q: %f, P: %f" % ( Q, P))

r1 = cavity.get_results("resonance")[-1]
print(r1['res']["xyprofile"].max_loc())
print(r1['res']["yzprofile"].max_loc())
r1["sess_res"].show()