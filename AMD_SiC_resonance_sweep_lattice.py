from wvgsolver import Cavity1D, UnitCell, Vec3
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine
import datetime as dt 
import matplotlib.pyplot as plt 
import matplotlib.animation as animation 

import scipy.optimize
import numpy as np
import os
import csv
from datetime import datetime


# animation function 
def animate(ax,xi,yi, xs, ys):
    """
    ax: the figure to plot the data on
    xi, yi: the input x and y values to be added to the plot 
    xs, ys: the x and y lists containing all the data
    """
    # Add x and y to lists
    xs.append(xi)
    ys.append(yi)
    
    # Draw x and y lists
    ax.clear()
    ax.plot(xs, ys)

    # Format plot
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('sim run')
    plt.ylabel('Fitness value')

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
    a = params[0] # sweeping the lattice constant 
    # define geometry parameters
    #taper cell number
    TN = 8
    #mirror cell number (the left mirror region)
    MN_L = 20
    #defect cell number
    CN = 0
    #lattice constant
    #a = 320e-9
    #hole diameter prefactor
    d = 0.64
    #beam width prefactor
    w = 1.75
    #taper prefactor
    t = 0.84
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
    #min taper hole diameter 
    rmin = (t*d*a)/2
    #radius taper rate
    r_tr = (r0-rmin) / TN
    #lattice taper rate
    a_tr = (a-amin) / TN
    #the refractive index of SiC 
    n_f = 2.6
    #V mode expected 
    Vmode_exp = 0.6
    #the prefactor characterizing the weaker mirror region 
    prefactor_mirror_R = 1
    #the right (weaker) mirror region 
    cellNum_R = 3
    #the waveguide region 
    waveguide_TN = 9
    # added waveguide region 
    t_wvg = t 
 
    # Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations
    FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
    # engine = LumericalEngine(mesh_accuracy=3, hide=True, lumerical_path=FDTDloc, working_path="./fsps", save_fsp=False)
    engine = LumericalEngine(mesh_accuracy=3, hide=True, lumerical_path=FDTDloc, save_fsp=False)

    ############################# building the mirror structure ###########################################################
    cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    mirror_cells_left = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * MN_L
    # added to modify the quasipotential associated with the right mirror region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    cell_box_R = BoxStructure(Vec3(0), Vec3(a_R,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    mirror_hole_R = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    mirror_cells_right = [UnitCell(structures=[ cell_box_R, mirror_hole_R ], size=Vec3(a), engine=engine)] * cellNum_R
    #########################################################################################################################
    
    ############### building cubic tapered cell region #####################################################
    taper_cells = []
    aList_taper = buildTapering_symmetric(a,t,TN)
    print(aList_taper) # debugging
    for i in aList_taper:
        taper_box = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        taper_hole = CylinderStructure(Vec3(0), h0, d*i/2, DielectricMaterial(1, order=1, color="blue"))
        taper_cells += [UnitCell(structures=[ taper_box, taper_hole ], size=Vec3(i), engine=engine)]
    ############################################################################################################

    ########################################### set the center of the device ###################################
    centerCell = MN_L+TN-1
    
    ################### adding one sided cubic tapered waveguide region to the cavity ##################################
    waveguide_cells_R = []
    a_wv = cubic_tapering(a,t_wvg,waveguide_TN)
    a_wv = a_wv[::-1]
    for i in a_wv:
        waveguide_box_R = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        waveguide_hole_R = CylinderStructure(Vec3(0), h0, d*i/2, DielectricMaterial(1, order=1, color="blue"))
        waveguide_cells_R += [UnitCell(structures=[ waveguide_box_R, waveguide_hole_R ], size=Vec3(i), engine=engine)]
    #############################################################################################################

    cavity = Cavity1D(
        unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right + waveguide_cells_R,
        structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
        engine=engine,
        center_cell=centerCell,
        center_shift=0,
    )

    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity.obj")

    #define mesh size (use 10nm for accuracy, currently set to 15nm)
    man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 15e-9, dy=None, dz=None)

    r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions = [man_mesh], sim_size=Vec3(4,4,12))

    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))

    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Qy = 1/(2/r1["qymax"])
    Qz = 1/(1/r1["qzmin"] + 1/r1["qzmax"])
    Vmode = r1["vmode"]
    F = r1["freq"]

    resonance_f = float(F) # the resonance frequency 
    resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
    detuning_wavelength = target_wavelength-resonance_wavelength
    
    Q = 1/((1/Qsc) + (1/Qwvg))
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f" % ( Q, P))
    print("Qy: %f, Qz: %f" % ( Qy, Qz))
    
    r1 = cavity.get_results("resonance")[0]

    # define the fitness as P with the resonance frequency Gaussian penalty
    fitness = (Q*Qsc)*detuning_wavelength
 
    print(a)
    
    with open("./sim_data/OptimizeListFull_resonance_lattice_sweep_v2.csv","a") as file_csv:
        writer = csv.writer(file_csv, delimiter="\t")
        writer.writerow([a,Q,Qsc,Qwvg,Vmode,detuning_wavelength,fitness])
    
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    
    return -1*fitness

p0 = [2.748043533042073e-07]
popt = scipy.optimize.minimize(fitness,p0,method='Nelder-Mead')

