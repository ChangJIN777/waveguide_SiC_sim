"""
This code loop through a number of initial guesses and see which one gives the best result 
"""

from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os
import scipy.constants
import csv 
from datetime import datetime
from tqdm import tqdm
from waveguideSolver_funcs import *

# in this script, we will use a for loop to iterate through different number of unit cells in the weak 
# mirror region and the waveguide region and calculate the properties of the resulting cavities 
def runSim(params):

    # Define geometry paramaters 
    #taper cell number (left mirror region)
    TN_L = 6
    TN_R = 6
    #mirror cell number (left region) 
    MN_L = 10
    #defect cell number
    CN = 0
    #lattice constant
    a = 2.843e-07
    #hole diameter in the x direction 
    hx = 7.24e-08
    #hole diameter in the y direction 
    hy = 1.275e-07
    #beam width prefactor
    w0 = 4.686e-07
    #taper prefactor (for the defect region)
    t = 0.818
    #taper prefactor (for the waveguide region)
    t_wvg = 0.852
    #beam height (set by epi-layer thickness)
    h0 = 250e-9
    # cavity beam length
    l = 15e-6   
    # The target resonance frequency, in Hz
    # 1280nm = 234.212857812500e12
    # 916nm = 327.3e12
    target_frequency = 327.3e12
    target_wavelength = 9.16e-07
    #set the center of the device
    centerCell = MN_L+TN_L-1
    # the number of unit cells in the weaker mirror region
    MN_R = params[0]
    # the taper cell number for the waveguide region 
    WN = params[1]
    #the minimum radius prefactor we are tapering to 
    d_min = 0.437
    hxmin_wvg = d_min*hx
    hymin_wvg = d_min*hy
    #the minimum lattice constant in the tapering region
    amin = a*t
    #the prefactor associated with the weaker mirror region
    prefactor_mirror_R = 0.965

    # Define geometry dependencies
    #refractive index of the material we are trying to simulate (SiC = 2.6)
    n_f = 2.6

    # Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
    FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
    engine = LumericalEngine(mesh_accuracy=5, hide=True, lumerical_path=FDTDloc, save_fsp=False)

    # building the cavity structure =======================================================================
    #build the left mirror cell region 
    mirror_cells_left = buildMirrorRegion_elliptical(a,hx,hy,MN_L,w0,h0,n_f,engine)

    #build the right mirror cell region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    mirror_cells_right = buildMirrorRegion_elliptical(a_R,hx,hy,MN_R,w0,h0,n_f,engine)

    #building cubic tapered cell region
    taper_cells = buildTaperRegion_elliptical_asymmetric(a,a_R,amin,hx,hy,TN_L,TN_R,w0,h0,n_f,engine)

    #add waveguide region 
    waveguide_cells = buildWaveguideRegion_elliptical_right_v2(a,hx,hxmin_wvg,hy,hymin_wvg,t_wvg,WN,w0,h0,n_f,engine)

    cavity = Cavity1D(
    unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right + waveguide_cells,
    structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
    center_cell=centerCell,
    center_shift=0,
    engine=engine
    )
    ##======================================================================================================
    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity.obj")

    #define mesh size (use 12nm for accuracy, currently set to 12nm)
    man_mesh = MeshRegion(BBox(Vec3(0),Vec3(5e-6,2e-6,2e-6)), 15e-9, dy=None, dz=None)

    # simulating the resonance and the Q =================================================
    r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
        analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(1.5,3,8))

    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))

    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Qxmin = r1["qxmin"]
    Qxmax = r1["qxmax"]
    Vmode = r1["vmode"]
    F = r1["freq"]
    resonance_f = float(F) # the resonance frequency 
    resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
    detuning_wavelength = target_wavelength-resonance_wavelength
    detuning_wavelength_nm = detuning_wavelength*1e9
    delta_wavelength = 5e-9 # 5nm tolerance 
    
    Q = 1/((1/Qsc) + (1/Qwvg))
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f" % ( Q, P))

    r1 = cavity.get_results("resonance")[0]
    
    # for debugging purposes  
    print("Qsc: %f Qwvg: %f" %(Qsc, Qwvg))

    fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((target_wavelength-resonance_wavelength)**2)/4))

    # # evaluate the quasipotential
    # r2 = cavity.simulate("quasipotential", target_freq=target_frequency)
    # r2.show()
    
    # writing the data into a csv file instead of a txt file for easier data analysis 
    data = [MN_R,WN,Vmode,Qwvg,Qsc,Qxmin,Qxmax,Q,F,detuning_wavelength,fitness]
    file_name = "OptimizeListFull_with_waveguide_params_loop_v4.csv"
    record_data(data,file_name)
    
    return -1*fitness

# limits on the number of unit cells in the weak mirror region
cellNum_R_min = 3
cellNum_R_max = 10
# limits on the number of unit cells in the waveguide region
waveguide_TN_min = 1
waveguide_TN_max = 8

# # loop through the list of parameters 
# for t_wvg in t_wvg_list:
#     for prefactor_mirror_R in prefactor_mirror_R_list:
#         for cellNum_R in range(cellNum_R_min,cellNum_R_max):
#             for waveguide_TN in tqdm(range(waveguide_TN_min,waveguide_TN_max)):
#                 test_params = [cellNum_R,waveguide_TN,prefactor_mirror_R,t_wvg]
#                 temp = runSim(test_params)
#                 print("the calculated fitness: %f" % (temp))


# # looping over only 2 parameters
# cellNum_R = 3
# waveguide_TN = 6
# for t_wvg in t_wvg_list:    
#     for prefactor_mirror_R in tqdm(prefactor_mirror_R_list):
#         test_params = [cellNum_R,waveguide_TN,prefactor_mirror_R,t_wvg]
#         temp = runSim(test_params)
#         print("the calculated fitness: %f" % (temp))

# looping over only 2 parameters  
for cellNum_R in range(cellNum_R_min,cellNum_R_max):
    for waveguide_TN in tqdm(range(waveguide_TN_min,waveguide_TN_max)):
        test_params = [cellNum_R,waveguide_TN]
        temp = runSim(test_params)
        print("the calculated fitness: %f" % (temp))
