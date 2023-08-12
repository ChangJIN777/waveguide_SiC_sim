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
sim_params["save_fsps"] = True

# improve the mesh resolution
sim_params["mesh_res"] = 15e-9

r1 = sim_ellipticalCavity_v2(cavity_params,sim_params)
