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
import csv 
from datetime import datetime

from waveguideSolver_funcs import *

# in this script, we will use a for loop to iterate through different number of unit cells in the weak 
# mirror region and the waveguide region and calculate the properties of the resulting cavities 
def runSim(params):
    print("Starting sim") # for debugging purpose
    start_time = datetime.now()
    # Define geometry paramaters 
    #taper cell number (left mirror region)
    TN = 8
    #mirror cell number (left region) 
    MN_L = 20
    #defect cell number
    CN = 0
    #beam height (set by epi-layer thickness)
    h0 = 250e-9
    # cavity beam length
    l = 15e-6
    # The target resonance frequency, in Hz
    # 1280nm = 234.212857812500e12
    # 916nm = 327.3e12
    target_frequency = 327.3e12
    target_wavelength = 9.16e-07
    # the number of unit cells in the weaker mirror region
    MN_R = 3
    # the taper cell number for the waveguide region 
    WN = 9
    # uncomment when you acutally want to sweep these parameters #############################################
    #lattice constant
    a = params[2]
    #hole diameter prefactor 
    d = params[3]
    #beam width prefactor
    w = params[4]
    #taper prefactor (for the defect region)
    t = params[5]
    ##########################################################################################################
    #the prefactor characterizing the right mirror region 
    prefactor_mirror_R = params[0]
    # added waveguide region 
    t_wvg = params[1] # set to 7.9075e-01 previously
    
    
    # Define geometry dependencies
    #beam width
    w0 = w*a
    # Radius of the air holes in the cells
    r0 = (d*a)/2
    #min tapered lattice 
    amin = t*a
    #refractive index of the material we are trying to simulate (SiC = 2.6)
    n_f = 2.6
    # added specifically for the waveguide region of the device =============
    #min tapered lattice in the waveguide region 
    amin_wvg = t_wvg*a
    #min taper hole radius 
    rmin_wvg = d*amin_wvg/2

    # Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
    FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
    engine = LumericalEngine(mesh_accuracy=4, hide=True, lumerical_path=FDTDloc, save_fsp=False)

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
    #build the cavity object 
    cavity = Cavity1D(
        unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right + waveguide_cells_R,
        structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
        engine=engine,
        center_cell=centerCell,
        center_shift=0,
    )
    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity.obj")

    #define mesh size (use 12nm for accuracy, currently set to 15)
    man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 15e-9, dy=None, dz=None)

    # simulating the resonance and the Q =================================================
    # r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions = [man_mesh], sim_size=Vec3(4,8,10), source_pulselength=200e-15, analyze_time=1000e-15)
    r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions = [man_mesh], sim_size=Vec3(4,8,10))
    #Note: specify a long pulse to make narrow band and target lossier mode closer to target_frequency
    
    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))
    
    # for debugging purposes 
    Qy =  1/(2/r1["qymax"])
    Qz = 1/(1/r1["qzmin"] + 1/r1["qzmax"])
     
    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Vmode = r1["vmode"]
    F = r1["freq"]
    resonance_f = float(F) # the resonance frequency 
    resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
    detuning_wavelength = target_wavelength-resonance_wavelength
    
    Q = 1/((1/Qsc) + (1/Qwvg))
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f" % ( Q, P))

    r1 = cavity.get_results("resonance")[0]
    
    # for debugging purposes  
    print("Qz: %f Qy: %f Qwvg: %f" %(Qz, Qy, Qwvg))

    fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((detuning_wavelength)**2)/25))

    # # evaluate the quasipotential
    # r2 = cavity.simulate("quasipotential", target_freq=target_frequency)
    # r2.show()
    
    # writing the data into a csv file instead of a txt file for easier data analysis 
    with open("./sim_data/OptimizeListFull_with_waveguide_test_sweep_allParam_v1.csv","a") as file_csv:
        writer = csv.writer(file_csv, delimiter="\t")
        writer.writerow([a,d,w,t,prefactor_mirror_R,t_wvg,Q,Qsc,Qwvg,Vmode,F,detuning_wavelength,fitness])

    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    
    return -1*fitness


p0 = [0.9, 0.84, 3.348909692268754e-07, 0.64, 1.75, 0.84]
# popt = scipy.optimize.minimize(runSim,p0,method='Nelder-Mead')

# # debugging 
temp = runSim(p0)