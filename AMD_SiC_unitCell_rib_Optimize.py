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

# define the function used for optimizing the rib unit cells 
def ribUnitCellOptimization(p):
    """this function is used to optimize the unit cells of rib cavities
    Args:
        params (list): 
            params[0]: the lattice constant a
            params[1]: hx
            params[2]: hy
    """
    rib_cavity_params["a"] = p[0]
    rib_cavity_params["hx"] = p[1]
    rib_cavity_params["hy"] = p[2]

    fitness = ribUnitCellOptimization_SiC(rib_cavity_params,rib_sim_params)
    return -1*fitness


# running the simulation locally 
rib_sim_params["running_cluster"] = True
rib_sim_params["running_local"] = False
rib_sim_params["hide_GUI"] = True

## DEBUGGING ##
# testing the bandgap simulation 
# rib_sim_params["hy"] = 3e-7
# sim_bandGap_rib(rib_cavity_params,rib_sim_params)
# band_structure_rib(rib_cavity_params,rib_sim_params)
# # sweep the dimensions of the rib unit cell 
rib_sim_params["simulationData_fileName"] = "SiC_500nm_rib_unitcell_testSweep_t3.txt"
# a_min = 2.6e-7 
# a_max = 2.7e-7 
# a_list = np.linspace(a_min,a_max,5)
# for a in a_list:
#     rib_cavity_params["a"] = a
#     sim_bandGap_rib(rib_cavity_params,rib_sim_params)

# optimizing for the mirror unit cells (SWEEPING CODE) ###################
a0 = rib_cavity_params["a"]
hx0 = rib_cavity_params["hx"]
hy0 = rib_cavity_params["hy"]
p0 = [a0,hx0,hy0]
bnd = ((2e-07,4e-07),(0,4e-7),(0,5e-7))
popt = scipy.optimize.minimize(ribUnitCellOptimization,p0,method='Nelder-Mead')



