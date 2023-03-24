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
a = 2.8333e-7
#hole diameter in the x direction 
hx = 7.0833e-8
#hole diameter in the y direction 
hy = 1.5937e-7
#beam width prefactor
w = 1.75
w0 = 5.25e-7
#taper prefactor (for the defect region)
t = 0.8
#taper prefactor (for the waveguide region)
t_wvg = 0.75
#beam height (set by epi-layer thickness)
h0 = 250e-9
# cavity beam length
l = 10e-6
# The target resonance frequency, in Hz
# 916nm = 327.3e12
target_frequency = 327.3e12
target_wavelength = 9.16e-07
#the prefactor associated with the weaker mirror region
prefactor_mirror_R = 0.9
#the refractive index associated with the material 
n_f = 2.6
#the minimum lattice constant in the waveguide region 
amin_wvg = t_wvg*a
# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations 
FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
# engine = LumericalEngine(mesh_accuracy=5, hide=False, lumerical_path=FDTDloc, working_path="./fsps")
engine = LumericalEngine(mesh_accuracy=5, hide=True, lumerical_path=FDTDloc, save_fsp=False)
#the minimum lattice constant in the tapering region
amin = a*t
#the minimum radius prefactor we are tapering to 
d_min = 0.3
#the left mirror cell number 
MN_L = 10 
#the right mirror cell number 
MN_R = 10
#the number of taper unit cells 
TN = 6
#set the center of the device
centerCell = MN_L+TN-1 

def run_Sim(param):
    print("Start sim ==============================")
    start_time = datetime.now()
    a = param[0]
    hx = param[1]
    hy = param[2]
    t = param[3]
    w0 = param[4]
    #build the left mirror cell region 
    mirror_cells_left = buildMirrorRegion_elliptical(a,hx,hy,MN_L,w0,h0,n_f,engine)

    #build the right mirror cell region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    mirror_cells_right = buildMirrorRegion_elliptical(a,hx,hy,MN_R,w0,h0,n_f,engine)

    #building cubic tapered cell region
    taper_cells = buildTaperRegion_elliptical(a,a_R,amin,hx,hy,TN,w0,h0,n_f,engine)

    ####################################### cavity without the waveguide region ###############################
    cavity = Cavity1D(
    unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right,
    structures=[ BoxStructure(Vec3(0), Vec3(l, w*a, h0), DielectricMaterial(n_f, order=2, color="red")) ],
    center_cell=centerCell,
    center_shift=0,
    engine=engine
    )
    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity_elliptical.obj")

    #define mesh size (use 12nm for accuracy, currently set to 12nm)
    # man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 12e-9, dy=None, dz=None)
    man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,2e-6,2e-6)), 20e-9, dy=None, dz=None)

    # simulating the resonance and the Q #########################################################
    # r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
    #                     analyze_time=1000e-15,analyze_fspan=5.0e12,mesh_regions = [man_mesh], sim_size=Vec3(4,8,8))
    r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
                        analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,4,8))

    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))

    cavity = Cavity1D(load_path="cavity_testing.obj",engine=engine)
    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Vmode = r1["vmode"]
    F = r1["freq"]
    resonance_f = float(F) # the resonance frequency 
    resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
    detuning_wavelength = target_wavelength-resonance_wavelength
    detuning_wavelength_nm = detuning_wavelength*1e9
    delta_wavelength = 5e-9 # 5nm tolerance 
    
    #prevent the mode volume from going to unrealistic values 
    if Vmode < 0.48:
        Vmode = 1e6
    
    Q = 1/((1/Qsc) + (1/Qwvg))

    if Q > 500000:
        Q = 500000
    
    if Qwvg > 500000:
        Qwvg = 500000
    
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f, detuning: %f nm" % ( Q, P, detuning_wavelength_nm))

    r1 = cavity.get_results("resonance")[-1]
    
    fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((detuning_wavelength/delta_wavelength)**2)))
    
    # record the data 
    data = [a,hx,hy,t,w0,prefactor_mirror_R,Vmode,Qwvg,Qsc,Q,F,detuning_wavelength,fitness]
    file_name = "OptimizeListFull_elliptical_cavity_sweep_v9.csv"
    record_data(data,file_name)
    
    end_time = datetime.now()
    print('fitness: %f'%(fitness))
    print('Duration: {}'.format(end_time - start_time))
    
    return -1*fitness

# # optimization algorithm
# p0 = [a,hx,hy,t,w0]
# bnd = [(None,None),(None,a),(None,w0),(None,1),(None,None)]
# popt = scipy.optimize.minimize(run_Sim,p0,method='Nelder-Mead')

# debugging 
# run_Sim(p0)

# # sweeping the taper prefactor 
# t_list = np.linspace(0.4,1,20)
# for t in t_list:
#     param = [t]
#     sweep_tapering_elliptical_cavity(param)

# sweeping the beam width  
w_list = np.linspace(4.00e-7,1.00e-8,20)
for w in w_list:
    param = [w]
    w_nm = w*1e9
    print("the beam width: %f nm" %(w_nm))
    sweep_beamWidth_ellipticalCavity_v2(param)

# # optimization algorithm (only the beam width)
# p0 = [w0]
# bnd = [(None,None)]
# popt = scipy.optimize.minimize(sweep_beamWidth_ellipticalCavity_v2,p0,method='Nelder-Mead')

# # sweeping the cell numbers 
# TN_list = [2,3,4,5,6,7,8]
# MN_L_list = [1,2,3,4,5,6,7,8,9,10]
# for TN in TN_list:
#     for MN_L in MN_L_list:
#         param = [MN_L,TN]
#         sweep_cellNum_ellipticalCavity(param)