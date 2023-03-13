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

#define the functions we are using to build the cavity geometry
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

def buildTapering_asymmetric(a_L,a_R,taperPrefac,taperNum):
    """
    function for calculating the lattice constants for the tapering region of the cavity 
    Note: this is used to build ASYMMETRIC taper cell region
    TN_L: the taper cell number on the left cell region 
    TN_R: the taper cell number of the right cell region 
    """
    a_taper_R = cubic_tapering(a_R,taperPrefac,taperNum=taperNum)
    a_taper_L = cubic_tapering(a_L,taperPrefac,taperNum=taperNum)
    a_taper_L = a_taper_L[::-1]
    tapering_region = np.concatenate((a_taper_L, a_taper_R), axis=None)
    return tapering_region

def buildWaveguideRegion_right(a,d,w,t_wvg,h0,WN,n_f,engine):
    """Function used to generate a cubic tapered waveguide region to be added to the right side of the cavity
    
    Args:
        a: the lattice constant used to build the mirror region 
        t_wvg: the tapering prefactor associated with the waveguide region 
        WN: the number of unit cells in the waveguide region 
        n_f: the refractive index associated with the material
        d: hole diameter prefactor 
        w: the beam width prefactor
        h0: the beam height
        engine: the FDTD engine used to simulate the waveguide region
    """
    w0 = w*a
    waveguide_cells_R = []
    a_wv = cubic_tapering(a,t_wvg,WN)
    print(a_wv) # debugging
    a_wv = a_wv[::-1]
    for i in a_wv:
        waveguide_box_R = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        waveguide_hole_R = CylinderStructure(Vec3(0), h0, d*i/2, DielectricMaterial(1, order=1, color="blue"))
        waveguide_cells_R += [UnitCell(structures=[ waveguide_box_R, waveguide_hole_R ], size=Vec3(i), engine=engine)]
    return waveguide_cells_R

def buildMirrorRegion(a,d,w,h0,n_f,MN,engine):
    """the function used to build mirriro region

    Args:
        a (float): the lattice constant used to build the mirror region 
        d (float): hole diameter prefactor 
        w (float): the beam width prefactor
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        MN (int): the number of mirror unit cells
        engine: the FDTD engine used to simulate the waveguide region
    """
    w0 = w*a #beam width
    r0 = d*a/2 #Radius of the air holes in the cells
    cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    mirror_cells = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * MN
    return mirror_cells

def buildTaperRegion(a_L,a_R,d,w,t,h0,n_f,TN,engine):
    """the function used to build mirriro region

    Args:
        a_L (float): the lattice constant used to build the mirror region on the left side 
        a_R (float): the lattice constant used to build the mirror region on the right side 
        d (float): hole diameter prefactor 
        w (float): the beam width prefactor
        t: the tapering prefactor associated with the waveguide region 
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        TN (int): the number of tapering cells 
        engine: the FDTD engine used to simulate the waveguide region
    """
    taper_cells = []
    w0 = w*a #the beam width
    aList_taper = buildTapering_asymmetric(a_L,a_R,t,TN)
    for i in aList_taper:
        taper_box = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        taper_hole = CylinderStructure(Vec3(0), h0, d*i/2, DielectricMaterial(1, order=1, color="blue"))
        taper_cells += [UnitCell(structures=[ taper_box, taper_hole ], size=Vec3(i), engine=engine)]
    return taper_cells

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
t_wvg = 0.84
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 15e-6
# The target resonance frequency, in Hz
# 1280nm = 234.212857812500e12
# 916nm = 327.3e12
target_frequency = 327.3e12
#the prefactor associated with the weaker mirror region
prefactor_mirror_R = 0.9
#the refractive index associated with the material 
n_f = 2.6
# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations ###################
FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
# engine = LumericalEngine(mesh_accuracy=5, hide=False, lumerical_path=FDTDloc, working_path="./fsps")
engine = LumericalEngine(mesh_accuracy=5, hide=True, lumerical_path=FDTDloc, save_fsp=False)

#build the left mirror cell region 
mirror_cells_left = buildMirrorRegion(a,d,w,h0,n_f,MN_L,engine)

#build the right mirror cell region 
a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
mirror_cells_right = buildMirrorRegion(a_R,d,w,h0,n_f,MN_R,engine)

#building cubic tapered cell region
taper_cells = buildTaperRegion(a,a_R,d,w,t,h0,n_f,TN,engine)

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

