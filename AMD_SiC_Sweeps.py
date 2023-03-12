"""
This example creates a cavity, and performs a series of simulations on it, visualizing or printing
out the results of each simulation as it completes.
"""

from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os
import csv
import struct

import matplotlib.pyplot as plt 
import matplotlib.animation as animation 
from matplotlib import style
from datetime import datetime


###################### define the functions we are using to build the cavity geometry ###################
def cubic_tapering(a,taperPrefac,taperNum):
    """
    a: the lattice constant in the mirror region
    taperNum: the number of taper cells
    taperPrefac: taper prefactor 
    """
    a_taper = np.zeros((taperNum,))
    d = 1-taperPrefac #defined as the depth of the defect (see Jasper Chan's thesis)
    for i in range(taperNum):
        a_taper[i] = a*(1-d*(2*((i/taperNum)**3) - 3*((i/taperNum)**2)+ 1))
    return a_taper

def buildTapering_symmetric(a,taperPrefac,taperNum):
    """
    function for calculating the lattice constants for the tapering region of the cavity 
    Note: this is used to build SYMMETRIC taper cell region
    """
    a_taper_R = cubic_tapering(a,taperPrefac,taperNum)
    a_taper_L = a_taper_R[::-1]
    tapering_region = np.concatenate((a_taper_L, a_taper_R), axis=None)
    return tapering_region

def buildTapering_asymmetric(a,taperPrefac,taperNum_L,taperNum_R):
    """
    function for calculating the lattice constants for the tapering region of the cavity 
    Note: this is used to build ASYMMETRIC taper cell region
    TN_L: the taper cell number on the left cell region 
    TN_R: the taper cell number of the right cell region 
    """
    a_taper_R = cubic_tapering(a,taperPrefac,taperNum=taperNum_R)
    a_taper_L = cubic_tapering(a,taperPrefac,taperNum=taperNum_L)
    a_taper_L = a_taper_L[::-1]
    tapering_region = np.concatenate((a_taper_L, a_taper_R), axis=None)
    return tapering_region


def fitness(params):
	print("Starting sim") # for debugging purpose
	start_time = datetime.now()
	a = params[0] # the lattice constant 
	d = params[1] # hole diameter prefactor 
	w = params[2] # beam width prefactor 
	t = params[3] # taper prefactor
	t_wvg = params[4] # the taper prefactor of the waveguide region
	# added parameter for simulating the waveguide 
	# the number of unit cells in the weaker mirror region
	cellNum_R = 3
	# the taper cell number for the waveguide region 
	waveguide_TN = 9
    
    # define geometry parameters
    #taper cell number
	TN = 8
    #mirror cell number (on the left side)
	MN = 20
    #defect cell number
	CN = 0
    #beam height (set by epi-layer thickness)
	h0 = 250e-9
    # cavity beam length
	l = 15e-6
    # The target resonance frequency, in Hz
	target_frequency = 327.3e12 # 327THz
	target_wavelength = 9.16e-07 # 916nm
 
    # Define geometry dependencies
    #beam width
	w0 = w*a
    # Radius of the air holes in the cells
	r0 = (d*a)/2
    #min tapered lattice
	amin = t*a
    #min taper hole radius 
	rmin = (t*d*a)/2
    #radius taper rate
	r_tr = (r0-rmin) / TN
    #lattice taper rate
	a_tr = (a-amin) / TN
	#refractive index of the material we are trying to simulate (SiC = 2.6)
	n_f = 2.6
	# added specifically for the waveguide region of the device =============
	#min tapered lattice in the waveguide region 
	amin_wvg = t_wvg*a
	#min taper hole radius 
	rmin_wvg = d*amin_wvg/2
	#radius taper rate (waveguide region)
	r_tr_wvg = (r0-rmin_wvg) / waveguide_TN
    #lattice taper rate (waveguide region)
	a_tr_wvg = (a-amin_wvg) / waveguide_TN

    # Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
	FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
	engine = LumericalEngine(mesh_accuracy=3, hide=True, lumerical_path=FDTDloc, save_fsp=False) # not saving the file 
	# engine = LumericalEngine(mesh_accuracy=3, hide=True, lumerical_path=FDTDloc, working_path="./fsps", save_fsp=False)
	
	# record the parameter list
	with open("OptimizeListFull_with_waveguide_parameters.csv","a") as file_csv:
		writer = csv.writer(file_csv, delimiter="\t")
		writer.writerow(params)

	cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
	mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
	mirror_cells_left = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * MN
	mirror_cells_right = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * cellNum_R
	cavity_cells = [UnitCell(structures=[ cell_box ], size=Vec3(amin), engine=engine)] * CN
    
	i = 1
	taper_cells_L = []
	taper_cells_R = []
	wvg_cells_R = []
	while i < TN: 
		taper_box_L = BoxStructure(Vec3(0), Vec3(a-(i*a_tr),w0,h0), DielectricMaterial(n_f, order=2, color="red"))
		taper_hole_L = CylinderStructure(Vec3(0), h0, r0-(i*r_tr), DielectricMaterial(1, order=1, color="blue"))
		taper_cells_L += [UnitCell(structures=[ taper_box_L, taper_hole_L ], size=Vec3(a-(i*a_tr)), engine=engine)]

		taper_box_R = BoxStructure(Vec3(0), Vec3(amin+(i*a_tr),w0,h0), DielectricMaterial(n_f, order=2, color="red"))
		taper_hole_R = CylinderStructure(Vec3(0), h0, rmin+(i*r_tr), DielectricMaterial(1, order=1, color="blue"))
		taper_cells_R += [UnitCell(structures=[ taper_box_R, taper_hole_R ], size=Vec3(amin+(i*a_tr)), engine=engine)]

		i = i+1 
	# construct the waveguide region 
	for i in range(waveguide_TN):
		wvg_box_R = BoxStructure(Vec3(0), Vec3(a-((i+1)*a_tr_wvg),w0,h0), DielectricMaterial(n_f, order=2, color="red"))
		wvg_hole_R = CylinderStructure(Vec3(0), h0, r0-((i+1)*r_tr_wvg), DielectricMaterial(1, order=1, color="blue"))
		wvg_cells_R += [UnitCell(structures=[ wvg_box_R, wvg_hole_R ], size=Vec3(a-((i+1)*a_tr_wvg)), engine=engine)]
 
 
	cavity = Cavity1D(
      unit_cells=  mirror_cells_left + taper_cells_L + taper_cells_R + mirror_cells_right  ,
      structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
      engine=engine
    )

    # By setting the save path here, the cavity will save itself after each simulation to this file
	cavity.save("cavity.obj")

    #define mesh size (use 10nm for accuracy, currently set to 12nm)
	man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 12e-9, dy=None, dz=None)

	r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions = [man_mesh], sim_size=Vec3(4,4,4))

    # Print the reults and plot the electric field profiles
	print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
    	r1["freq"], r1["vmode"],
    	1/(1/r1["qxmin"] + 1/r1["qxmax"]),
    	1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))

	Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
	Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
	Vmode = r1["vmode"]
	F = r1["freq"]
 
	resonance_f = float(F) # the resonance frequency 
	resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
	detuning_wavelength = target_wavelength-resonance_wavelength
 
	Q = 1/((1/Qsc) + (1/Qwvg))
	P = (Q*Qsc) / (Vmode*Vmode)
	print("Q: %f, P: %f" % ( Q, P))

	r1 = cavity.get_results("resonance")[0]

	print(a)
	fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((detuning_wavelength)**2)/25))

	# file = open("OptimizeListFull_waveguide_sweep.txt","a") 
	# #file.write("\n" + str(a) + " " + str(Q) + " " + str(Vmode)+ " " + str(F) + "\n")
	# file.write("\n" + str(params) + "\t" + str(Q) + "\t" + str(Vmode)+ "\t" + str(F) + "\n")
	# file.close()
	# writing the data into a csv file instead of a txt file for easier data analysis 
	with open("OptimizeListFull_with_waveguide_char.csv","a") as file_csv:
		writer = csv.writer(file_csv, delimiter="\t")
		writer.writerow([Q,Vmode,F,detuning_wavelength,fitness])


	return -1*fitness

p0 = [3.348909692268754e-7, 0.64, 1.75, 0.84, 0.84]
bnds = ((0,None),(0,None),(0,None),(0,1),(0,1))
popt = scipy.optimize.minimize(fitness,p0,method='Nelder-Mead',bounds=bnds)




