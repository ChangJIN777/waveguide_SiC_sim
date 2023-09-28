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
cavity_params = cavity_sim_parameters.cavity_params
sim_params = cavity_sim_parameters.sim_params
sim_params["show_field_profile"] = True
sim_params["hide_GUI"] = False
sim_params["save_fsps"] = True
sim_params["running_cluster"] = True   
sim_params["running_local"] = False
sim_params["boundary_condition"] = ['zmin','xmin']
# testing the side coupling code 
cavity_params["do_sc"] = True
cavity_params["M_lattice_prefactor"] = 1
cavity_params["MN_Right"] = 10
cavity_params["WN"] = 0

# improve the mesh resolution
sim_params["mesh_res"] = 12e-9

r1 = sim_ellipticalCavity_v2(cavity_params,sim_params)
