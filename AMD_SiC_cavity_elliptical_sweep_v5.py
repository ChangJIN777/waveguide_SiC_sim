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
from cavity_sim_parameters import *

from waveguideSolver_funcs import *

#filename we are saving the data under 
file_name = "holeDimSweep_500nm_trial1.csv"

# debugging 
#hx_min = d_min*hx
#hy_min = d_min*hy
#sim_ellipticalCavity_v2(a,hx,hx_min,hy,hy_min,t,t_wvg,target_frequency,file_name,MN_L,MN_R,TN,w0=w0,engine=engine)

# import the cavity and simulation parameters 
cavity_params = cavity_sim_parameters.cavity_params
sim_params = cavity_sim_parameters.sim_params
cavity_params["simulationData_fileName"] = file_name

# sweeping the unit cell parameters
hx_list = np.linspace(50e-09,100e-09,10)
hy_list = np.linspace(100e-09,150e-09,10)
for hx in hx_list:
    for hy in hy_list:
        cavity_params["hx"] = hx
        cavity_params["hy"] = hy
        sim_ellipticalCavity_v2(cavity_params,sim_params)

# # sweeping the taper prefactor 
# t_list = np.linspace(0.4,1,20)
# for t in t_list:
#     param = [t]
#     sweep_tapering_elliptical_cavity(param)

# sweeping the beam width  
# w_list = np.linspace(4.00e-7,8.00e-7,20)
# for w in w_list:
#     param = [w]
#     w_nm = w*1e9
#     print("the beam width: %f nm" %(w_nm))
#     sweep_beamWidth_ellipticalCavity_v2(param)

# optimize the beam width  
# p0 = [w0]
# bnd = [(None,1e-6)]
# popt = scipy.optimize.minimize(sweep_beamWidth_ellipticalCavity_v2,p0,method='Nelder-Mead')
# sweep_beamWidth_ellipticalCavity_v2(param)


# # optimization algorithm (only the beam width)
# p0 = [w0]
# bnd = [(None,None)]
# popt = scipy.optimize.minimize(sweep_beamWidth_ellipticalCavity_v2,p0,method='Nelder-Mead')

# # sweeping the cell numbers 
# TN_list = [2,3,4,5,6,7,8]
# MN_L_list = [1,2,3,4,5,6,7,8,9,10]
# for TN in TN_list:
#     for MN_L in MN_L_list:
#         param = [MN_L,TN]
#         sweep_cellNum_ellipticalCavity(param)