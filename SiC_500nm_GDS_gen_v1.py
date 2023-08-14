"""
This script is used to generate gds files from the existing designs 
"""
    
from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine
from wvgsolver.parse import DielectricExtrusionFaceGDSParser

import scipy.optimize
import numpy as np
import os
from datetime import datetime

from waveguideSolver_funcs import *
from cavity_sim_parameters import * 

# import the cavity and simulation parameters 
cavity_params = cavity_sim_parameters.cavity_params
sim_params = cavity_sim_parameters.sim_params
#list of lattice constant we want to sweep through 
a_min = (cavity_params["a"])*0.9
a_max = (cavity_params["a"])*1.1
a_number = 5
a_list = np.linspace(a_min,a_max,a_number)

# # generate GDS from the design
# for i in range(len(a_list)):
#     for j in range(len(hy_min_list)):
#         cavity_temp = build_cavity_v5(a_list[i],hy_min_list[j])
#         parser = DielectricExtrusionFaceGDSParser(cavity_temp)
#         parser.show()
#         file_name = "SiC_cavity_v5_a_%d_hy_%d.gds"%(i,j)
#         parser.save(file_name)
        
# generate GDS from the design
for i in range(len(a_list)):
    cavity_temp = build_cavity_500nm_v1(cavity_params,sim_params)
    parser = DielectricExtrusionFaceGDSParser(cavity_temp)
    parser.show()
    file_name = "SiC_500nm_cavity_v3_a_%d.gds"%(i)
    parser.save(file_name)