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
rib_sim_params["save_fsps"] = True
rib_sim_params["hide_GUI"] = False
# testing the side coupling code 
rib_cavity_params["do_sc"] = False
rib_sim_params["running_cluster"] = True
rib_sim_params["running_local"] = False

# improve the mesh resolution
rib_sim_params["mesh_res"] = 12e-9
rib_sim_params["a"] = 2.659218106995883e-07
rib_sim_params["boundary_condition"] = ['ymin']
rib_cavity_params["C_lattice_tapering_prefactor"] = 0.85

r1 = sim_rib_Cavity_v1(rib_cavity_params,rib_sim_params)
