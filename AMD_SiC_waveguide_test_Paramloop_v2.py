"""
This code loop through a number of initial guesses and see which one gives the best result 
"""

from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os
import scipy.constants
import csv 
from datetime import datetime
from tqdm import tqdm

# in this script, we will use a for loop to iterate through different number of unit cells in the weak 
# mirror region and the waveguide region and calculate the properties of the resulting cavities 
def runSim(params):

    # Define geometry paramaters 
    #taper cell number (left mirror region)
    TN = 8
    #mirror cell number (left region) 
    MN_L = 24-TN
    #defect cell number
    CN = 0
    #lattice constant
    a = 2.921072969062500e-07
    #hole diameter prefactor 
    d = 0.697201546893750
    #beam width prefactor
    w = 1.825228557093750
    #taper prefactor (for the defect region)
    t = 0.684551582507813
    #beam height (set by epi-layer thickness)
    h0 = 220e-9
    # cavity beam length
    l = 15e-6
    # The target resonance frequency, in Hz
    # 1280nm = 234.212857812500e12
    # 916nm = 327.3e12
    target_frequency = 327.3e12
    target_wavelength = 9.16e-07
    # the number of unit cells in the weaker mirror region
    cellNum_R = params[0]
    # the taper cell number for the waveguide region 
    waveguide_TN = params[1]
    #the prefactor characterizing the right mirror region 
    prefactor_mirror_R = params[2]
    #the prefactor characterizing the waveguide region 
    t_wvg = params[3]

    # Define geometry dependencies
    #beam width
    w0 = w*a
    # Radius of the air holes in the cells
    r0 = (d*a)/2
    #min tapered lattice 
    amin = t*a
    #min taper hole diameter 
    rmin = (t*d*a)/2
    #radius taper rate (for the defect taper region)
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
    engine = LumericalEngine(mesh_accuracy=4, hide=True, lumerical_path=FDTDloc, save_fsp=False)

    # building the cavity structure =======================================================================
    # the sim material is set to be SiC with refractive index = 2.6 
    cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    mirror_cells_left = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * MN_L
    # added to modify the quasipotential associated with the right mirror region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    cell_box_R = BoxStructure(Vec3(0), Vec3(a_R,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    mirror_hole_R = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    mirror_cells_right = [UnitCell(structures=[ cell_box_R, mirror_hole_R ], size=Vec3(a), engine=engine)] * cellNum_R
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
    unit_cells=  mirror_cells_left + taper_cells_L + taper_cells_R + mirror_cells_right + wvg_cells_R,
    structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
    engine=engine
    )
    ##======================================================================================================
    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity.obj")

    #define mesh size (use 12nm for accuracy, currently set to 12nm)
    man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 15e-9, dy=None, dz=None)

    # simulating the resonance and the Q =================================================
    r1 = cavity.simulate("resonance", target_freq=target_frequency, mesh_regions = [man_mesh], sim_size=Vec3(4,4,10))

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
    
    # prevent the mode volume from going to an unrealistic value 
    Vmode_exp = 0.6
    
    Q = 1/((1/Qsc) + (1/Qwvg))
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f" % ( Q, P))

    r1 = cavity.get_results("resonance")[0]
    
    # for debugging purposes  
    print("Qsc: %f Qwvg: %f" %(Qsc, Qwvg))

    fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((target_wavelength-resonance_wavelength)**2)/4))

    # # evaluate the quasipotential
    # r2 = cavity.simulate("quasipotential", target_freq=target_frequency)
    # r2.show()
    
    # writing the data into a csv file instead of a txt file for easier data analysis 
    with open("./sim_data/OptimizeListFull_with_waveguide_params_loop_v1.csv","a") as file_csv:
        writer = csv.writer(file_csv, delimiter="\t")
        writer.writerow([cellNum_R,waveguide_TN,prefactor_mirror_R,t_wvg,Q,Qsc,Qwvg,Vmode,F,detuning_wavelength,fitness])
    
    return -1*fitness

# limits on the number of unit cells in the weak mirror region
cellNum_R_min = 3
cellNum_R_max = 10
# limits on the number of unit cells in the waveguide region
waveguide_TN_min = 3
waveguide_TN_max = 10
# taper factor of the weak mirror region
prefactor_mirror_R_min = 0.6
prefactor_mirror_R_max = 1    
prefactor_mirror_R_list = np.linspace(prefactor_mirror_R_min,prefactor_mirror_R_max,5)    
# taper factor of the waveguide region
t_wvg_min = 0.6
t_wvg_max = 1
t_wvg_list = np.linspace(t_wvg_min,t_wvg_max,5)

# loop through the list of parameters 
for t_wvg in tqdm(t_wvg_list):
    for prefactor_mirror_R in prefactor_mirror_R_list:
        for cellNum_R in range(cellNum_R_min,cellNum_R_max):
            for waveguide_TN in range(waveguide_TN_min,waveguide_TN_max):
                test_params = [cellNum_R,waveguide_TN,prefactor_mirror_R,t_wvg]
                temp = runSim(test_params)
                print("the calculated fitness: %f" % (temp))

