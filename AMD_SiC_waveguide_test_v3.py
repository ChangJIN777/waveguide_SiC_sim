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
import scipy.constants

# in this script, we will be using the traditional method of putting different number of mirror unit cells on each side of the defect region

# Define geometry paramaters 
#waveguide taper cell number
WN = 6 
#taper cell number (left mirror region)
TN = 8
#mirror cell number (left region) 
MN_L = 24-TN
#mirror cell number (right region)
MN_R = 12-TN
#defect cell number
CN = 0
#lattice constant
a = 290e-9
#hole diameter prefactor 
d = 0.64
#beam width prefactor
w = 1.75
#taper prefactor (for the defect region)
t = 0.7
#beam height (set by epi-layer thickness)
h0 = 220e-9
# cavity beam length
l = 15e-6
# The target resonance frequency, in Hz
# 1280nm = 234.212857812500e12
# 916nm = 327.3e12
target_frequency = 327.3e12
target_wavelength = 9.16e-07

# Define geometry dependencies

#beam width
w0 = w*a
# Radius of the air holes in the cells
r0 = (d*a)/2
#min tapered lattice 
amin = t*a
#min taper hole diameter 
rmin = (t*d*a)/2
#min tapered lattice (waveguide)
amin_wvg = t_wvg*a 
#min taper hole diameter (waveguide)
rmin_wvg = d*amin_wvg/2
#radius taper rate (for the defect taper region)
r_tr = (r0-rmin) / TN
#lattice taper rate
a_tr = (a-amin) / TN
#radius taper rate (for the waveguide region)
r_wvg_tr = (r0-rmin_wvg)/WN
#lattice taper rate (for the waveguide region)
a_wvg_tr = (a-amin_wvg)/WN
#refractive index of the material we are trying to simulate (SiC = 2.6)
n_f = 2.6

# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
engine = LumericalEngine(mesh_accuracy=4, hide=False, lumerical_path=FDTDloc, working_path="./fsps", save_fsp=False)

# the sim material is set to be SiC with refractive index = 2.6 
cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
mirror_cells_left = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * MN_L
mirror_cells_right = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * MN_R
cavity_cells = [UnitCell(structures=[ cell_box ], size=Vec3(amin), engine=engine)] * CN

i = 1
taper_cells_L = []
taper_cells_R = []
while i < TN: 
    taper_box_L = BoxStructure(Vec3(0), Vec3(a-(i*a_tr),w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    taper_hole_L = CylinderStructure(Vec3(0), h0, r0-(i*r_tr), DielectricMaterial(1, order=1, color="blue"))
    taper_cells_L += [UnitCell(structures=[ taper_box_L, taper_hole_L ], size=Vec3(a-(i*a_tr)), engine=engine)]

    taper_box_R = BoxStructure(Vec3(0), Vec3(amin+(i*a_tr),w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    taper_hole_R = CylinderStructure(Vec3(0), h0, rmin+(i*r_tr), DielectricMaterial(1, order=1, color="blue"))
    taper_cells_R += [UnitCell(structures=[ taper_box_R, taper_hole_R ], size=Vec3(amin+(i*a_tr)), engine=engine)]

    i = i+1 

cavity = Cavity1D(
unit_cells=  mirror_cells_left + taper_cells_L + taper_cells_R + mirror_cells_right,
structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
engine=engine
)

# By setting the save path here, the cavity will save itself after each simulation to this file
cavity.save("cavity.obj")

#define mesh size (use 12nm for accuracy, currently set to 50nm)
man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 10e-9, dy=None, dz=None)

# simulating the resonance and the Q =================================================
r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions = [man_mesh], sim_size=Vec3(4,4,4))

# Print the reults and plot the electric field profiles
print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
r1["freq"], r1["vmode"],
1/(1/r1["qxmin"] + 1/r1["qxmax"]),
1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
))
r1["xyprofile"].show()
r1["yzprofile"].show()

Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
Vmode = r1["vmode"]
F = r1["freq"]
resonance_f = float(F) # the resonance frequency 
resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 

Q = 1/((1/Qsc) + (1/Qwvg))
P = (Q*Qsc) / (Vmode*Vmode)
print("Q: %f, P: %f" % ( Q, P))

# debugging 
print("Qsc: %f Qwvg: %f" %(Qsc, Qsc))

Vmode_value = Vmode*((resonance_f/n_f)**3)
fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((target_wavelength-resonance_wavelength)**2)/25))

r1 = cavity.get_results("resonance")[0]
print(r1['res']["xyprofile"].max_loc())
print(r1['res']["yzprofile"].max_loc())
print("Fitness %f"%(fitness))
r1["sess_res"].show()
# ======================================================================================

# # evaluate the quasipotential
# r2 = cavity.simulate("quasipotential", target_freq=target_frequency)
# r2.show()

# file = open("OptimizeList.txt","a") 
# file.write("\n" + str(a) + " " + str(Q) + " " + str(Vmode)+ " " + str(F) + "\n") 
# file.close()
