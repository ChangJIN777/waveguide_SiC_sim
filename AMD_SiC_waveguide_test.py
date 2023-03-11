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

# Define geometry paramaters 
#waveguide taper cell number
WN = 7
#taper cell number (left mirror region)
TN = 8
#mirror cell number (left region) 
MN_L = 20
#mirror cell number (right region)
MN_R = 4
#defect cell number
CN = 0
#lattice constant
a = 2.748043533042073e-07 
#hole diameter prefactor 
d = 0.64
#beam width prefactor
w = 1.75
#taper prefactor (for the defect region)
t = 0.84
#taper prefactor (for the waveguide region)
t_wvg = 0.84
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 15e-6
# The target resonance frequency, in Hz
# 1280nm = 234.212857812500e12
# 916nm = 327.3e12
target_frequency = 327.3e12

###################### define the functions we are using to build the cavity geometry ###################
def cubic_tapering(a,taperPrefac,taperNum):
    """
    a: the lattice constant in the mirror region
    taperNum: the number of taper cells
    taperPrefac: taper prefactor 
    """
    a_taper = np.zeros((taperNum,))
    d = 1-taperPrefac #defined as the depth of the defect (see Jasper Chan's thesis)
    for i in range(taperNum):
        a_taper[i] = a*(1-d*(2*((i/taperNum)**3) - 3*((i/taperNum)**2)+ 1))
    return a_taper

def buildTapering_symmetric(a,taperPrefac,taperNum):
    """
    function for calculating the lattice constants for the tapering region of the cavity 
    Note: this is used to build SYMMETRIC taper cell region
    """
    a_taper_R = cubic_tapering(a,taperPrefac,taperNum)
    a_taper_L = a_taper_R[::-1]
    tapering_region = np.concatenate((a_taper_L, a_taper_R), axis=None)
    return tapering_region

def buildTapering_asymmetric(a,taperPrefac,taperNum_L,taperNum_R):
    """
    function for calculating the lattice constants for the tapering region of the cavity 
    Note: this is used to build ASYMMETRIC taper cell region
    TN_L: the taper cell number on the left cell region 
    TN_R: the taper cell number of the right cell region 
    """
    a_taper_R = cubic_tapering(a,taperPrefac,taperNum=taperNum_R)
    a_taper_L = cubic_tapering(a,taperPrefac,taperNum=taperNum_L)
    a_taper_L = a_taper_L[::-1]
    tapering_region = np.concatenate((a_taper_L, a_taper_R), axis=None)
    return tapering_region

#################################### Define geometry dependencies ######################################
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
############################################################################################################

# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations ###################
FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
engine = LumericalEngine(mesh_accuracy=5, hide=False, lumerical_path=FDTDloc, working_path="./fsps")

# the sim material is set to be SiC with refractive index = 2.6 
cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(2.6, order=2, color="red"))
mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
mirror_cells_left = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * MN_L
mirror_cells_right = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * MN_R
cavity_cells = [UnitCell(structures=[ cell_box ], size=Vec3(amin), engine=engine)] * CN
###########################################################################################################

# ############### building linearly tapered cell region #####################################################
# i = 1
# taper_cells_L = []
# taper_cells_R = []
# while i <= TN: 
#     taper_box_L = BoxStructure(Vec3(0), Vec3(a-(i*a_tr),w0,h0), DielectricMaterial(2.6, order=2, color="red"))
#     taper_hole_L = CylinderStructure(Vec3(0), h0, r0-(i*r_tr), DielectricMaterial(1, order=1, color="blue"))
#     taper_cells_L += [UnitCell(structures=[ taper_box_L, taper_hole_L ], size=Vec3(a-(i*a_tr)), engine=engine)]

#     taper_box_R = BoxStructure(Vec3(0), Vec3(amin+(i*a_tr),w0,h0), DielectricMaterial(2.6, order=2, color="red"))
#     taper_hole_R = CylinderStructure(Vec3(0), h0, rmin+(i*r_tr), DielectricMaterial(1, order=1, color="blue"))
#     taper_cells_R += [UnitCell(structures=[ taper_box_R, taper_hole_R ], size=Vec3(amin+(i*a_tr)), engine=engine)]
#     i = i+1 
# ###################################################################################################################

# ############### building cubic tapered cell region #####################################################
taper_cells = []
aList_taper = buildTapering_symmetric(a,t,TN)
print(aList_taper) # debugging
for i in aList_taper:
    taper_box = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(2.6, order=2, color="red"))
    taper_hole = CylinderStructure(Vec3(0), h0, d*i/2, DielectricMaterial(1, order=1, color="blue"))
    taper_cells += [UnitCell(structures=[ taper_box, taper_hole ], size=Vec3(i), engine=engine)]
############################################################################################################

########################################### set the center of the device ###################################
centerCell = MN_L+TN-1
# centerCell = 0 # testing the code
###########################################################################################################

# ################### adding one sided linearly tapered waveguide region to the cavity ##################################
# waveguide_cells_R = []
# for i in range(WN):
#     waveguide_box_R = BoxStructure(Vec3(0), Vec3(a-((i+1)*a_wvg_tr),w0,h0), DielectricMaterial(2.6, order=2, color="red"))
#     waveguide_hole_R = CylinderStructure(Vec3(0), h0, r0-((i+1)*r_wvg_tr), DielectricMaterial(1, order=1, color="blue"))
#     waveguide_cells_R += [UnitCell(structures=[ waveguide_box_R, waveguide_hole_R ], size=Vec3(a-((i+1)*a_wvg_tr)), engine=engine)]
# #############################################################################################################

################### adding one sided cubic tapered waveguide region to the cavity ##################################
waveguide_cells_R = []
a_wv = cubic_tapering(a,t_wvg,WN)
print(a_wv) # debugging
a_wv = a_wv[::-1]
for i in a_wv:
    waveguide_box_R = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(2.6, order=2, color="red"))
    waveguide_hole_R = CylinderStructure(Vec3(0), h0, d*i/2, DielectricMaterial(1, order=1, color="blue"))
    waveguide_cells_R += [UnitCell(structures=[ waveguide_box_R, waveguide_hole_R ], size=Vec3(i), engine=engine)]
#############################################################################################################


cavity = Cavity1D(
unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right + waveguide_cells_R ,
structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(2.6, order=2, color="red")) ],
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
                    analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,4,10))

# Print the reults and plot the electric field profiles
print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
r1["freq"], r1["vmode"],
1/(1/r1["qxmin"] + 1/r1["qxmax"]),
1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
))
r1["xyprofile"].show()
r1["yzprofile"].show()

cavity = Cavity1D(load_path="cavity_testing.obj",engine=engine)
# Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
# Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
# Vmode = r1["vmode"]
# F = r1["freq"]

# Q = 1/((1/Qsc) + (1/Qwvg))
# P = (Q*Qsc) / (Vmode*Vmode)
# print("Q: %f, P: %f" % ( Q, P))

r1 = cavity.get_results("resonance")[-1]
print(r1['res']["xyprofile"].max_loc())
print(r1['res']["yzprofile"].max_loc())
# r1["sess_res"].show()
# # ======================================================================================

# # # evaluate the quasipotential
# r2 = cavity.simulate("quasipotential", target_freq=target_frequency)
# r2.show()

# file = open("OptimizeList.txt","a") 
# file.write("\n" + str(a) + " " + str(Q) + " " + str(Vmode)+ " " + str(F) + "\n") 
# file.close()

