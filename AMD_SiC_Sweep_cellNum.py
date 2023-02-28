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


def fitness(params):

    # added parameter for simulating the waveguide
    cellNum_R = params[0]
    cellNum_R = int(cellNum_R)
    # define geometry parameters
    # taper cell number
    TN = 8
    # mirror cell number (on the left side)
    # MN = 24-TN
    cellNum_L = 24-TN
    # defect cell number
    CN = 0
    # lattice constant
    a = 2.97688965e-07
    # hole diameter prefactor
    d = 6.63014844e-01
    # beam width prefactor
    w = 1.73572998e+00
    # taper prefactor
    t = 7.48911133e-01
    # beam height (set by epi-layer thickness)
    h0 = 220e-9
    # cavity beam length
    l = 15e-6
    # The target resonance frequency, in Hz
    target_frequency = 327.3e12 # 327THz
    target_wavelength = 9.16e-07 # 916nm

    # Define geometry dependencies
    # beam width
    w0 = w*a
    # Radius of the air holes in the cells
    r0 = (d*a)/2
    # min tapered lattice
    amin = t*a
    # min taper hole diameter 
    rmin = (t*d*a)/2
    # radius taper rate
    r_tr = (r0-rmin) / TN
    # lattice taper rate
    a_tr = (a-amin) / TN
    # refractive index of the material we are trying to simulate (SiC = 2.6)
    n_f = 2.6

    # Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
    FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
    engine = LumericalEngine(mesh_accuracy=3, hide=True, lumerical_path=FDTDloc, save_fsp=False) # not saving the file 
    # engine = LumericalEngine(mesh_accuracy=3, hide=True, lumerical_path=FDTDloc, working_path="./fsps", save_fsp=False)
    
    # record the parameter list
    with open("OptimizeListFull_waveguide_Sweep_cellNum.csv","a") as file_csv:
        writer = csv.writer(file_csv, delimiter="\t")
        writer.writerow(params)

    cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    mirror_cells_left = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * cellNum_L
    mirror_cells_right = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * cellNum_R
    cavity_cells = [UnitCell(structures=[ cell_box ], size=Vec3(amin), engine=engine)] * CN
    
    i = 1
    taper_cells_L = []
    taper_cells_R = []
    while i < TN: 
        taper_box_L = BoxStructure(Vec3(0), Vec3(a-(i*a_tr),w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        taper_hole_L = CylinderStructure(Vec3(0), h0, r0-(i*r_tr), DielectricMaterial(1, order=1, color="blue"))
        taper_cells_L += [UnitCell(structures=[ taper_box_L, taper_hole_L ], size=Vec3(a-(i*a_tr)), engine=engine)]

        taper_box_R = BoxStructure(Vec3(0), Vec3(amin+(i*a_tr),w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        taper_hole_R = CylinderStructure(Vec3(0), h0, rmin+(i*r_tr), DielectricMaterial(1, order=1, color="blue"))
        taper_cells_R += [UnitCell(structures=[ taper_box_R, taper_hole_R ], size=Vec3(amin+(i*a_tr)), engine=engine)]

        i = i+1 

    cavity = Cavity1D(
      unit_cells=  mirror_cells_left + taper_cells_L + taper_cells_R + mirror_cells_right ,
      structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
      engine=engine
    )

    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity.obj")

    # define mesh size (use 10nm for accuracy, currently set to 12nm)
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
    with open("OptimizeListFull_waveguide_Sweep_char.csv","a") as file_csv:
        writer = csv.writer(file_csv, delimiter="\t")
        writer.writerow([Q,Vmode,F,detuning_wavelength,fitness])

    return -1*fitness

for i in range(3,16):
    fitness = fitness(i)
        

# p0 = [5,16]
# bnds = ((3,16),(3,30))
# popt = scipy.optimize.minimize(fitness,p0,method='Nelder-Mead',bounds=bnds)

