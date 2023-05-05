"""
This script is used to generate gds files from the existing designs 
"""
    
from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine
from wvgsolver.parse import DielectricExtrusionFaceGDSParser

import scipy.optimize
import numpy as np
import os
from datetime import datetime

from waveguideSolver_funcs import *
    
#list of lattice constant we want to sweep through 
a_min = (2.843000000000000e-07)*0.9
a_max = (2.843000000000000e-07)*1.1
a_number = 11
a_list = np.linspace(a_min,a_max,a_number)

#hole diameter in the x direction
hx = 72.40e-09
#hole diameter in the y direction
hy = 127.5e-09
#beam width
w0 = 468.6e-09
#taper prefactor (for the defect region)
t = 0.818
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 15e-6
# The target resonance frequency, in Hz
target_frequency = 327.3e12
#the prefactor associated with the weaker mirror region
prefactor_mirror_R = 0.965
#taper prefactor (for the waveguide region)
t_wvg = 0.852
#the refractive index associated with the material 
n_f = 2.6
#the minimum lattice constant in the waveguide region 
amin_wvg = t_wvg*a
# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations 
FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
engine = LumericalEngine(mesh_accuracy=5, hide=False, lumerical_path=FDTDloc, working_path="./fsps")
#the minimum lattice constant in the tapering region
amin = a*t
#the minimum radius prefactor we are tapering to 
d_min = 0.437
hx_min = d_min*hx
hy_min = d_min*hy
#the left mirror cell number 
MN_L = 10
#the right mirror cell number 
MN_R = 3
#the number of taper unit cells 
TN = 5
TN_L = 8
TN_R = 4
#the number of waveguide cells 
WN = 5
# list of minimum tapered hole width 
hy_min_list = [60e-9,120e-9,hy_min,180e-9]

#set the center of the device (for double sided cavities)
centerCell = MN_L+TN_L

def build_cavity_v6(a,hy_min):
    """_summary_

    Args:
        a (float): the lattice constant associated with the cavity

    Returns:
        lumerical object: cavity object built from the parameters given
    """
    # the lattice constant associated with the right mirror region 
    a_R = a*prefactor_mirror_R 
    #build the left mirror cell region 
    mirror_cells_left = buildMirrorRegion_elliptical(a,hx,hy,MN_L,w0,h0,n_f,engine)

    #build the right mirror cell region 
    mirror_cells_right = buildMirrorRegion_elliptical(a_R,hx,hy,MN_R,w0,h0,n_f,engine)

    #building cubic tapered cell region
    taper_cells = buildTaperRegion_elliptical_asymmetric(a,a_R,amin,hx,hy,TN_L,TN_R,w0,h0,n_f,engine)

    #add waveguide region (v1)
    waveguide_cells_R = buildWaveguideRegion_elliptical_right_v2(a,hx,hx_min,hy,hy_min,t_wvg,WN,w0,h0,n_f,engine)

    ####################################### cavity with the waveguide region (asymmetric) ###############################
    cavity = Cavity1D(
    unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right + waveguide_cells_R,
    structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
    center_cell=centerCell,
    center_shift=0,
    engine=engine
    )
    
    return cavity

# # generate GDS from the design
# for i in range(len(a_list)):
#     for j in range(len(hy_min_list)):
#         cavity_temp = build_cavity_v5(a_list[i],hy_min_list[j])
#         parser = DielectricExtrusionFaceGDSParser(cavity_temp)
#         parser.show()
#         file_name = "SiC_cavity_v5_a_%d_hy_%d.gds"%(i,j)
#         parser.save(file_name)
        
# generate GDS from the design
for i in range(len(a_list)):
    cavity_temp = build_cavity_v6(a_list[i],hy_min)
    parser = DielectricExtrusionFaceGDSParser(cavity_temp)
    parser.show()
    file_name = "SiC_cavity_v6_a_%d.gds"%(i)
    parser.save(file_name)