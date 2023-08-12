"""
This code sweeps the tapering prefactors of the defect and the waveguide regions to get an overcoupled cavity 
"""

from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os
from datetime import datetime
from cavity_sim_parameters import *

from waveguideSolver_funcs import *

#filename we are saving the data under 
file_name = "cavity_500nm_optimization_trial1.csv"
# import the cavity and simulation parameters 
cavity_params = cavity_sim_parameters.cavity_params
sim_params = cavity_sim_parameters.sim_params
# setup the simulation 
sim_params["simulationData_fileName"] = file_name
sim_params["hide_GUI"] = True
sim_params["save_fsps"] = False

def optimization(params):
    cavity_params["a"] = params[0] # lattice constant 
    cavity_params["hx"] = params[1] 
    cavity_params["hy"] = params[2]
    r1 = sim_ellipticalCavity_v2(cavity_params,sim_params)
    
    fitness = calculate_fitness(r1,sim_params)
    
    print('fitness: %f'%(fitness))
    
    return -1*fitness

# optimization algorithm
p0 = [cavity_params["a"],cavity_params["hx"],cavity_params["hy"]]
bnd = [(None,300e-9),(50e-9,cavity_params["a"]),(100e-9,cavity_params["a"])]
popt = scipy.optimize.minimize(optimization,p0,method='Nelder-Mead')