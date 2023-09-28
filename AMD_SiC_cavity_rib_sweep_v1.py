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
rib_cavity_params = cavity_sim_parameters.rib_cavity_params
rib_sim_params = cavity_sim_parameters.rib_sim_params
rib_sim_params["show_field_profile"] = False
rib_sim_params["save_fsps"] = False
rib_sim_params["hide_GUI"] = True
# testing the side coupling code 
rib_cavity_params["do_sc"] = False
rib_sim_params["running_cluster"] = True 
rib_sim_params["running_local"] = False

# simulation settings 
rib_sim_params["mesh_res"] = 12e-9
rib_sim_params["boundary_condition"] = ['ymin','zmin']
rib_sim_params["mode"] = 'TM'
rib_sim_params["simulationData_fileName"] = "SiC_500nm_rib_testRun_t5.csv"

# geometry settings
rib_cavity_params["hx"] = 1.8e-07 # for the TM mode 
rib_cavity_params["hy"] = 3.84e-07 # for the TM mode 
rib_cavity_params["beam_width"] = 4e-07 # for the TM mode 
rib_cavity_params["C_lattice_tapering_prefactor"] = 0.867
a = 3e-07
hx = 1.8e-07 
hy = 3.84e-07 
w0 = 4e-07
a_min = a*0.95
a_max = a*1.05
a_list = np.linspace(a_min,a_max,10)
for a in a_list:
    rib_cavity_params["a"] = a
    r1 = sim_rib_Cavity_v1(rib_cavity_params,rib_sim_params)
