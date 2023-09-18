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

# running the simulation locally 
rib_sim_params["running_cluster"] = False
rib_sim_params["running_local"] = True
rib_sim_params["hide_GUI"] = False

# testing the bandgap simulation 
rib_sim_params["hy"] = 3e-7
sim_bandGap_rib(rib_cavity_params,rib_sim_params)
# band_structure_rib(rib_cavity_params,rib_sim_params)
# # sweep the dimensions of the rib unit cell 
# rib_sim_params["simulationData_fileName"] = "SiC_500nm_rib_unitcell_testSweep_t1.txt"
# hy_min = 3e-7 
# hy_max = 5e-7 
# hy_list = np.linspace(hy_min,hy_max,10)
# for hy in hy_list:
#     rib_cavity_params["hy"] = hy
#     sim_bandGap_rib(rib_cavity_params,rib_sim_params)



