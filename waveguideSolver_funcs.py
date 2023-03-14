"""
This file contains all the useful functions for building a overcoupled cavity
"""

from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os
from datetime import datetime
import csv

#define the useful constants 
n_f = 2.6 # for SiC
target_frequency = 327.3e12


#define the functions we are using to build the cavity geometry
def cubic_tapering(a,amin,taperNum):
    """
    a: the lattice constant in the mirror region
    taperNum: the number of taper cells
    amin: the minimum lattice constant in the tapering region
    """
    a_taper = np.zeros((taperNum,))
    d = 1-(amin/a) #defined as the depth of the defect (see Jasper Chan's thesis)
    for i in range(taperNum):
        a_taper[i] = a*(1-d*(2*((i/taperNum)**3) - 3*((i/taperNum)**2)+ 1))
    return a_taper

def buildTapering_symmetric(a,amin,taperNum):
    """
    function for calculating the lattice constants for the tapering region of the cavity 
    Note: this is used to build SYMMETRIC taper cell region
    """
    a_taper_R = cubic_tapering(a,amin,taperNum)
    a_taper_L = a_taper_R[::-1]
    tapering_region = np.concatenate((a_taper_L, a_taper_R), axis=None)
    return tapering_region

def buildTapering_asymmetric(a_L,a_R,amin,taperNum):
    """
    function for calculating the lattice constants for the tapering region of the cavity 
    Note: this is used to build ASYMMETRIC taper cell region
    amin: the minimum lattice constant in the tapering region
    taperNum: the number of tapering cells
    """
    a_taper_R = cubic_tapering(a_R,amin,taperNum=taperNum)
    a_taper_L = cubic_tapering(a_L,amin,taperNum=taperNum)
    a_taper_L = a_taper_L[::-1]
    tapering_region = np.concatenate((a_taper_L, a_taper_R), axis=None)
    return tapering_region

def buildWaveguideRegion_right(a,d,w,t_wvg,h0,WN,n_f,engine):
    """Function used to generate a cubic tapered waveguide region to be added to the right side of the cavity
    
    Args:
        a: the lattice constant used to build the mirror region 
        t_wvg: the tapering prefactor associated with the waveguide region 
        WN: the number of unit cells in the waveguide region 
        n_f: the refractive index associated with the material
        d: hole diameter prefactor 
        w: the beam width prefactor
        h0: the beam height
        engine: the FDTD engine used to simulate the waveguide region
    """
    w0 = w*a
    waveguide_cells_R = []
    amin = t_wvg*a
    a_wv = cubic_tapering(a,amin,WN)
    a_wv = a_wv[::-1]
    print(a_wv) # debugging
    for i in a_wv:
        waveguide_box_R = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        waveguide_hole_R = CylinderStructure(Vec3(0), h0, d*i/2, DielectricMaterial(1, order=1, color="blue"))
        waveguide_cells_R += [UnitCell(structures=[ waveguide_box_R, waveguide_hole_R ], size=Vec3(i), engine=engine)]
    return waveguide_cells_R

def buildMirrorRegion(a,d,w,h0,n_f,MN,engine):
    """the function used to build mirriro region

    Args:
        a (float): the lattice constant used to build the mirror region 
        d (float): hole diameter prefactor 
        w (float): the beam width prefactor
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        MN (int): the number of mirror unit cells
        engine: the FDTD engine used to simulate the waveguide region
    """
    w0 = w*a #beam width
    r0 = d*a/2 #Radius of the air holes in the cells
    cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    mirror_cells = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)] * MN
    return mirror_cells

def buildTaperRegion(a_L,a_R,amin,d,w,h0,n_f,TN,engine):
    """the function used to build mirriro region

    Args:
        a_L (float): the lattice constant used to build the mirror region on the left side 
        a_R (float): the lattice constant used to build the mirror region on the right side 
        d (float): hole diameter prefactor 
        w (float): the beam width prefactor
        amin: the minimum lattice constant in the tapering region
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        TN (int): the number of tapering cells 
        engine: the FDTD engine used to simulate the waveguide region
    """
    taper_cells = []
    w0 = w*a_L #the beam width
    aList_taper = buildTapering_asymmetric(a_L,a_R,amin,TN)
    for i in aList_taper:
        taper_box = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        taper_hole = CylinderStructure(Vec3(0), h0, d*i/2, DielectricMaterial(1, order=1, color="blue"))
        taper_cells += [UnitCell(structures=[ taper_box, taper_hole ], size=Vec3(i), engine=engine)]
    return taper_cells

def buildUnitCell(a,d,w,h0,n_f,engine):
    """the function use the given parameters to build a unit cell 

    Args:
        a (float): the lattice constant 
        d (float): hole diameter prefactor 
        w (float): the beam width prefactor
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        engine : the FDTD engine used to simulate the waveguide region

    Returns:
        cell: the UnitCell object in waveguide solver 
    """
    w0 = w*a
    r0 = a*d/2
    cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    cell = UnitCell(structures=[ cell_box, hole ], size=Vec3(a), engine=engine)
    return cell

def sim_bandGap(a,d,w,h0,n_f,engine):
    """the function generates the bandgap associated with the simulated unit cell

    Args:
        a (float): the lattice constant 
        d (float): hole diameter prefactor 
        w (float): the beam width prefactor
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        engine : the FDTD engine used to simulate the waveguide region
        
    Returns:
        _type_: _description_
    """
    start_time = datetime.now()
    cell = buildUnitCell(a,d,w,h0,n_f,engine)

    r2 = cell.simulate("bandgap", freqs=(0.15e15, 0.4e15, 100000))

    diel_freq = r2[0] # the dielectric band frequency 
    air_freq = r2[1] # the air band frequyency 
    bg = air_freq - diel_freq # the band gap 
    mg = (diel_freq + air_freq) / 2 # the mid band gap 
    bg_mg_rat = bg / mg 

    delta_k = .5-a*mg/2.998e8
    end_time = datetime.now()

    print('--------')
    print('Duration: {}'.format(end_time - start_time))
    print('Lower band edge frequency: %f THz' % (diel_freq / 1e12))
    print("Bandgap ratio: %f percent" % (bg_mg_rat * 100))
    print("Midgap frequency: %f THz" % (mg / 1e12))
    print("Delta k: %f " % delta_k)
    print('\n')

    return diel_freq, air_freq, mg, bg_mg_rat, delta_k

def band_structure(a,d,w,h0,n_f,engine):
    
    start_time = datetime.now()
    cell = buildUnitCell(a,d,w,h0,n_f,engine)
    # r1 = cell.simulate("bandstructure", ks=(0.2, 0.5, 12), freqs=(0.25e15, 0.7e15, 100000),
    #                    dipole_region=Vec3(0.8, 0, 0), window_pos = 0)
    r1 = cell.simulate("bandstructure", ks=(0.2, 0.5, 8), freqs=(0.15e15, 0.4e15, 150000))
    # # # Plot the bandstructure
    r1.show()
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    
def unitCellOptimization_SiC(params):
    print("Starting sim") # for debugging purpose
    a = params[0]
    d = params[1]
    w = params[2]
    h0 = params[3]
    # Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations 
    FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
    engine = LumericalEngine(mesh_accuracy=5, hide=True, lumerical_path=FDTDloc, save_fsp=False)
    # simulate the band gap of the unit cell 
    diel_freq, air_freq, mg, bg_mg_rat, delta_k = sim_bandGap(a,d,w,h0,n_f,engine)
    wavelength_pen = np.exp(-((target_frequency - mg)/target_frequency)**2) # the wavelength detuning penalty
    detuning = target_frequency - mg
    fitness = -1*bg_mg_rat*wavelength_pen
    
    a_nm = a*1e9
    h0_nm = h0*1e9
    print("a: %f (nm), d: %f, w: %f, h0: %f (nm)" % (a_nm, d, w, h0_nm))
    print("Detune: %f" % (detuning))
    print("Fitness: %f" % (fitness))
    
    # writing the data into a csv file instead of a txt file for easier data analysis 
    with open("./sim_data/unitcell_Optimization.csv","a") as file_csv:
        writer = csv.writer(file_csv, delimiter="\t")
        writer.writerow([a,d,w,h0,mg,bg_mg_rat,fitness])
    
    return fitness