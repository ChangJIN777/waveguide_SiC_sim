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
from cavity_sim_parameters import *


# import the cavity and simulation parameters 
sim_params = cavity_sim_parameters.rib_sim_params
cavity_params = cavity_sim_parameters.rib_cavity_params

# # define the function used for optimizing the rib unit cells 
# def ribUnitCellOptimization(p):
#     """this function is used to optimize the unit cells of rib cavities
#     Args:
#         params (list): 
#             params[0]: the lattice constant a
#             params[1]: hx
#             params[2]: hy
#     """
#     rib_cavity_params["a"] = p[0]
#     rib_cavity_params["hx"] = p[1]
#     rib_cavity_params["hy"] = p[2]

#     fitness = (rib_cavity_params,rib_sim_params)
#     return -1*fitness


# running the simulation locally 
sim_params["running_cluster"] = False
sim_params["running_local"] = True
sim_params["hide_GUI"] = True
sim_params["save_fsps"] = False


## DEBUGGING ##
# testing the bandgap simulation 
# rib_sim_params["hy"] = 3e-7
# sim_bandGap_rib(rib_cavity_params,rib_sim_params)
# band_structure_rib(rib_cavity_params,rib_sim_params)

# # sweep the dimensions of the rib unit cell 
sim_params["simulationData_fileName"] = "Si_220nm_rib_unitcell_testSweep_TM_t1.txt"
a = 3.e-07
hx = 1.8e-07 # for the TM mode 
hy = 3.84e-07 # for the TM mode 
a_min = a*0.5
a_max = a
hx_min = hx*0.8
hx_max = hx*1.2
hy_min = hy*0.8
hy_max = hy*1.2
a_list = np.linspace(a_min,a_max,15)
hx_list = np.linspace(hx_min,hx_max,10)
hy_list = np.linspace(hy_min,hy_max,10)
sim_data_folder = sim_params["simulationData_loc"]
sim_data_fileName = sim_params["simulationData_fileName"]
for a in a_list:
    cavity_params['a'] = a
    diel_freq, air_freq, mg, bg_mg_rat, delta_k, bg = sim_bandGap_rib(cavity_params,sim_params)
    data = [a, diel_freq, air_freq, mg, bg_mg_rat, delta_k, bg]
    record_data(data,sim_data_fileName,sim_data_folder)


# # optimizing for the mirror unit cells (SWEEPING CODE) ###################
# a0 = rib_cavity_params["a"]
# hx0 = rib_cavity_params["hx"]
# hy0 = rib_cavity_params["hy"]
# p0 = [a0,hx0,hy0]
# bnd = ((2e-07,4e-07),(0,4e-7),(0,5e-7))
# popt = scipy.optimize.minimize(ribUnitCellOptimization,p0,method='Nelder-Mead')



