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
sim_params = cavity_sim_parameters.sim_params
cavity_params = cavity_sim_parameters.cavity_params

# running the simulation locally 
sim_params["running_cluster"] = False
sim_params["running_local"] = True
sim_params["hide_GUI"] = True
sim_params["save_fsps"] = False

# # sweep the dimensions of the holely unit cell 
sim_params["simulationData_fileName"] = "Si_220nm_hoe_unitcell_testSweep_TE_t2.txt"
a = 346.5e-9 
hx = 96e-9
hy = 96e-9
thickness = 220e-9
beam_width = 284.3e-9
n_f = 3.5
a_min = a
a_max = a*1.3
hx_min = hx*0.8
hx_max = hx*1.2
hy_min = hy*0.8
hy_max = hy*1.2
cavity_params['hx'] = hx
cavity_params['hy'] = hy
cavity_params['thickness'] = thickness
cavity_params['beam_width'] = beam_width
a_list = np.linspace(a_min,a_max,5)
hx_list = np.linspace(hx_min,hx_max,10)
hy_list = np.linspace(hy_min,hy_max,10)
sim_data_folder = sim_params["simulationData_loc"]
sim_data_fileName = sim_params["simulationData_fileName"]
engine, mesh = setup_engine(sim_params)
# for a in a_list:
#     cavity_params['a'] = a
#     diel_freq, air_freq, mg, bg_mg_rat, delta_k, bg = sim_bandGap_elliptical(a,hx,hy,beam_width,thickness,n_f,engine)
#     data = [a, diel_freq, air_freq, mg, bg_mg_rat, delta_k, bg]
#     record_data(data,sim_data_fileName,sim_data_folder)

# optimizing the unit cell
# # debugging 
# p0 = [a,hx,hy]
# ellipseUnitCellOptimization_Si(p0)

# optimizing for the mirror unit cells (SWEEPING CODE) ###################
p0 = [a,hx,hy]
bnd = ((0,None),(None,beam_width/2),(None,beam_width/2))
popt = scipy.optimize.minimize(ellipseUnitCellOptimization_Si,p0,method='Nelder-Mead')
