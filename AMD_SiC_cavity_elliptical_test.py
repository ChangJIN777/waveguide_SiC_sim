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

r1 = sim_ellipticalCavity_v2(cavity_params,sim_params)

# Print the reults and plot the electric field profiles
print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
    r1["freq"], r1["vmode"],
    1/(1/r1["qxmin"] + 1/r1["qxmax"]),
    1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
))
r1["xyprofile"].show()
r1["yzprofile"].show()

cavity = Cavity1D(load_path="cavity.obj",engine=engine)
Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
Vmode = r1["vmode"]
F = r1["freq"]

# for debugging the Q 
Qx1 = r1["qxmin"]
Qx2 = r1["qxmax"]
Qy = 1 / (2 / r1["qymax"])
Qz = 1 / (1 / r1["qzmin"] + 1 / r1["qzmax"])
print("Qx1: %f, Qx2: %f, Qy: %f, Qz: %f" % (
    Qx1, Qx2, Qy, Qz
))

Q = 1/((1/Qsc) + (1/Qwvg))
P = (Q*Qsc) / (Vmode*Vmode)
print("Q: %f, P: %f" % ( Q, P))

r1 = cavity.get_results("resonance")[-1]
print(r1['res']["xyprofile"].max_loc())
print(r1['res']["yzprofile"].max_loc())
r1["sess_res"].show()