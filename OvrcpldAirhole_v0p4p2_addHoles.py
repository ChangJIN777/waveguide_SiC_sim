"""
This example creates a cavity, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.utils import BBox
from wvgsolver.geometry import BoxStructure, TriStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.engine import LumericalEngine
from wvgsolver.parse import DielectricExtrusionFaceGDSParser, ObjectParser3D
from waveguideSolver_funcs import *

import numpy as np
import os

#  Initialize Lumerical File Locations
FDTDLoc = '/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/'
FDTDexeLoc = os.path.join(FDTDLoc,'bin/fdtd-solutions')
FDTDmpiLoc = os.path.join(FDTDLoc,'bin/fdtd-engine-ompi-lcl')


hole_radii = np.loadtxt(os.path.join(os.path.curdir, 'v0p4p2/holeStruct.txt'), dtype=float, usecols=(0,1), unpack=False)
hole_radii /= 2
lattice_constants = np.loadtxt(os.path.join(os.path.curdir, 'v0p4p2/periodStruct.txt'), dtype=float, usecols=(0), unpack=False)
center_cell = 12
# Unit cells are triangular prisms
beam_width = 0.482e-6
apex_half_angle = 50*np.pi/180
# lattice_constant = 0.270e-6
beam_height = (beam_width / 2) * np.tan(np.pi/2 - apex_half_angle)
# Radius of the air holes in the cells
# hole_radius = 0.5*0.120e-6
# The length of the cavity beam
beam_length = 15e-6
# The target resonance frequency, in Hz
target_frequency = 406.774e12

# Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
engine = LumericalEngine(mesh_accuracy=5, hide=False, lumerical_path=FDTDLoc, working_path="./fsps")
# engine = LumericalEngine(mesh_accuracy=6, hide=False, working_path="./fsps")

unit_cells = []
for ((radius_x,radius_y), lattice_constant) in zip(hole_radii,lattice_constants):
  temp_unitcell = buildUnitCell_elliptical(lattice_constant,radius_x,radius_y,beam_width,beam_height,engine=engine)
  unit_cells += [temp_unitcell]

cavity = Cavity1D(
  unit_cells=unit_cells,
  structures=[BoxStructure(Vec3(0), Vec3(beam_length,beam_width,beam_height), 
                        DielectricMaterial(2.417, order=2,color="red"))],
  engine=engine,
  center_cell=center_cell,
  center_shift=0
)

# parsed = ObjectParser3D(cavity)
# parsed.show()

# By setting the save path here, the cavity will save itself after each simulation to this file
cavity.save("v0p4p2p2_3-3_120221.obj")

man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,2e-6,2e-6)), 12e-9, dy=None, dz=None)

# specify a long pulse to make narrow band and target lossier mode closer to target_frequency
r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
                    analyze_time=1000e-15, mesh_regions = [man_mesh], sim_size=Vec3(2, 4, 8))
# r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions = [man_mesh], sim_size=Vec3(2, 8, 14))

# Print the reults and plot the electric field profiles
print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
  r1["freq"], r1["vmode"],
  1/(1/r1["qxmin"] + 1/r1["qxmax"]),
  1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
))

r1["xyprofile"].show()
r1["yzprofile"].show()

cavity = Cavity1D(load_path="v0p4p2p2_3-3_120221.obj",engine=engine)
r1 = cavity.get_results("resonance")[-1]
print(r1)
print(r1['res']["xyprofile"].max_loc())
print(r1['res']["yzprofile"].max_loc())
r1["sess_res"].show()

print(r1)

# r2 = cavity.simulate("quasipotential", target_freq=target_frequency)

# # Plot the quasipotential
# r2.show()

# r3 = cavity.simulate("guidedness", target_freq=target_frequency)

# print("Guidedness: %f" % r3)

