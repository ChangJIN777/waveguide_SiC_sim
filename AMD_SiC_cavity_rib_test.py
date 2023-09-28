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
rib_sim_params["show_field_profile"] = True
rib_sim_params["save_fsps"] = False
rib_sim_params["hide_GUI"] = False

# testing the side coupling code 
rib_cavity_params["do_sc"] = True
rib_sim_params["running_cluster"] = True  
rib_sim_params["running_local"] = False

# improve the mesh resolution
rib_sim_params["mesh_res"] = 12e-9
# rib_cavity_params["a"] = 2.836499314128942e-07 # for the TE mode 
rib_cavity_params["a"] = 3e-07 # for the TM mode 
rib_cavity_params["hx"] = 1.8e-07 # for the TM mode 
rib_cavity_params["hy"] = 3.84e-07 # for the TM mode 
rib_cavity_params["beam_width"] = 4e-07 # for the TM mode 
rib_sim_params["boundary_condition"] = ['ymin','zmin']
rib_sim_params["mode"] = 'TM'
rib_cavity_params["C_lattice_tapering_prefactor"] = 0.867
rib_sim_params["freq_span"] = 3e14 
rib_sim_params["simulationData_fileName"] = "SiC_500nm_debugging_cavity.csv"

# testing the rib cavity
r1 = sim_rib_Cavity_v1(rib_cavity_params,rib_sim_params)

# sweep the target wavelength 
# rib_sim_params["show_field_profile"] = False
# freq_span = 10e12 
# central_freq = rib_sim_params["target_frequency"]
# freq_start = central_freq - freq_span
# freq_end = central_freq + freq_span
# freq_lists = np.linspace(freq_start,freq_end,5)
# for freq in freq_lists:
#     rib_sim_params["target_frequency"] = freq 
#     r1 = sim_rib_Cavity_v1(rib_cavity_params,rib_sim_params)
