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

#define the functions we are using to build the cavity geometry
def linear_tapering(a,amin,taperNum):
    """the function to return an array of linearly tapered parameters

    Args:
        a: the lattice constant in the mirror region
        taperNum: the number of taper cells
        amin: the minimum lattice constant in the tapering region

    Returns:
        a_taper: _description_
    """
    a_taper = np.zeros((taperNum,))
    d = 1-(amin/a) #defined as the depth of the defect (see Jasper Chan's thesis)
    for i in range(taperNum):
        a_taper[i] = a*(1-d*(i/taperNum))
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

def buildTapering_asymmetric_v2(a_L,a_R,amin,TN_L,TN_R):
    """
    function for calculating the lattice constants for the tapering region of the cavity 
    Note: this is used to build ASYMMETRIC taper cell region 
    amin: the minimum lattice constant in the tapering region
    TN_L: the number of tapering cells on the left side of the tapering region 
    TN_R: the number of tapering cells on the right side of the tapering region 
    a_L: the lattice constant on the left side of the tapering region 
    a_R: the lattice constant on the right side of the tapering region
    """
    a_taper_R = cubic_tapering(a_R,amin,taperNum=TN_R)
    a_taper_L = cubic_tapering(a_L,amin,taperNum=TN_L)
    a_taper_L = a_taper_L[::-1]
    tapering_region = np.concatenate((a_taper_L, a_taper_R), axis=None)
    return tapering_region

def buildWaveguideRegion_right(a,hx,hy,w0,t_wvg,WN,h0,n_f,engine):
    """Function used to generate a cubic tapered waveguide region to be added to the right side of the cavity
    
    Args:
        a: the lattice constant used to build the mirror region 
        hx: the radius in the x direction 
        hy: the radius in the y direction 
        t_wvg: the tapering prefactor associated with the waveguide region 
        WN: the number of unit cells in the waveguide region 
        n_f: the refractive index associated with the material
        h0: the beam height
        engine: the FDTD engine used to simulate the waveguide region
    """
    waveguide_cells_R = []
    amin = t_wvg*a
    a_wv = cubic_tapering(a,amin,WN)
    a_wv = a_wv[::-1]
    print(a_wv) # debugging
    for i in a_wv:
        waveguide_box_R = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        waveguide_hole_R = CylinderStructure(Vec3(0), h0, d*i/2, DielectricMaterial(1, order=1, color="blue"))
        waveguide_cells_R += [UnitCell(structures=[ waveguide_box_R, waveguide_hole_R ], size=Vec3(i,w0,h0), engine=engine)]
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
    mirror_cells = [UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a,w0,h0), engine=engine)] * MN
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
        taper_cells += [UnitCell(structures=[ taper_box, taper_hole ], size=Vec3(i,w0,h0), engine=engine)]
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
    cell = UnitCell(structures=[ cell_box, hole ], size=Vec3(a,w0,h0), engine=engine)
    return cell

def sim_bandGap_elliptical(a,hx,hy,w0,h0,n_f,engine):
    """the function generates the bandgap associated with the simulated unit cell

    Args:
        a (float): the lattice constant 
        hx (float): hole diameter in the x direction
        hy (float): hole diameter in the y direction
        w0 (float): the beam width
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        engine : the FDTD engine used to simulate the waveguide region
        
    Returns:
        _type_: _description_
    """
    print("Starting sim ===================================")
    start_time = datetime.now()
    cell = buildUnitCell_elliptical(a,hx,hy,w0,h0,n_f,engine)

    r2 = cell.simulate("bandgap", freqs=(0.15e15, 0.6e15, 100000))

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
    print('Higher band edge frequency: %f THz' % (air_freq / 1e12))
    print("Bandgap ratio: %f percent" % (bg_mg_rat * 100))
    print("Midgap frequency: %f THz" % (mg / 1e12))
    print("Delta k: %f " % delta_k)
    print('\n')

    return diel_freq, air_freq, mg, bg_mg_rat, delta_k, bg

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
    cell = buildUnitCell(a,d,w,h0,n_f,engine=engine)

    r2 = cell.simulate("bandgap", freqs=(0.15e15, 0.5e15, 100000))

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
    r1 = cell.simulate("bandstructure", ks=(0.2, 0.5, 8), freqs=(0.15e15, 0.5e15, 150000))
    # # # Plot the bandstructure
    r1.show()
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    
def unitCellOptimization_SiC(params):
    """this function is used to optimize the unit cells for the mirror regions of the cavity

    Args:
        params (list): 
            params[0]: the lattice constant 
            params[1]: the radius prefactor 
            params[2]: the beam width prefactor

    Returns:
        fitness: the optimization parameter 
    """
    print("Starting sim ===================") # for debugging purpose
    a = params[0]
    d = params[1]
    w = params[2]
    # Use level 4 automeshing accuracy, and show the Lumerical GUI while running simulations 
    FDTDloc="/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/"
    engine = LumericalEngine(mesh_accuracy=5, hide=True, lumerical_path=FDTDloc, save_fsp=False)
    # simulate the band gap of the unit cell 
    diel_freq, air_freq, mg, bg_mg_rat, delta_k = sim_bandGap(a,d,w,h0,n_f,engine)
    wavelength_tolerance = 5e-9
    wavelength_detune = (3e8)/(target_frequency)-(3e8)/(mg)
    wavelength_pen = np.exp(-((wavelength_detune)/wavelength_tolerance)**2) # the wavelength detuning penalty
    detuning = target_frequency - mg
    fitness = -1*bg_mg_rat*wavelength_pen
    
    a_nm = a*1e9
    h0_nm = h0*1e9
    print("a: %f (nm), d: %f, w: %f, h0: %f (nm)" % (a_nm, d, w, h0_nm))
    print("Detune: %f" % (detuning))
    print("Fitness: %f" % (fitness))
    
    # writing the data into a csv file instead of a txt file for easier data analysis 
    with open("./sim_data/unitcell_Optimization_v2.csv","a") as file_csv:
        writer = csv.writer(file_csv, delimiter="\t")
        writer.writerow([a,d,w,h0,mg,bg_mg_rat,fitness])
    
    return fitness

# simulating the unit cell band structures
def bandStructSim(a,d,w,n_f,fmin,fmax,f_grating,kmin,kmax,knum,engine,h0):
    """this function simulate the band structure of a unit cell with parameters given by the inputs
        a: the lattice constant 
        d: hole diameter prefactor 
        w: beam width prefactor
        h0: beam height 
        n_f: the refractive index of the dielectric material 
        engine: the lumerical engine used to simulate the structure
    """
    #beam width
    w0 = w*a
    # Radius of the air holes in the cells
    r0 = (d*a)/2

    # building the unit cell 
    # the sim material is set to be SiC with refractive index = 2.6 
    cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    mirror_hole = CylinderStructure(Vec3(0), h0, r0, DielectricMaterial(1, order=1, color="blue"))
    simulate_unit_cell = UnitCell(structures=[ cell_box, mirror_hole ], size=Vec3(a), engine=engine)
    
    simulate_unit_cell.simulate("bandstructure", ks=(kmin,kmax,knum),freqs=(fmin, fmax, f_grating))

    return 

def unitCellOptimization_SiC(params):
    """this function is used to optimize the unit cells for the waveguide regions of the cavity

    Args:
        params (list): 
            params[0]: the tapered lattice constant (this should be the value the waveguide is tapering to)

    Returns:
        fitness: the optimization parameter (the detuning of the dielectric band and the target frequency)
    """
    print("Starting sim ===================") # for debugging purpose
    a = params[0]
    # simulate the band gap of the unit cell 
    diel_freq, air_freq, mg, bg_mg_rat, delta_k = sim_bandGap(a,d_0,w_0,h0,n_f,engine)
    detuning = np.abs(target_frequency - diel_freq)
    print("Detuning from the dielectric band: %f"%(detuning))
    file_name = "unitcell_Optimization_waveguide_v1.csv"
    data = [a,detuning]
    record_data(data,file_name)
    return detuning

def record_data(data,file_name,file_loc):
    denstination = file_loc + file_name
    with open(denstination,"a") as file_csv:
        writer = csv.writer(file_csv, delimiter="\t")
        writer.writerow(data)
        
def buildWaveguideRegion_right_v2(a,d,d_min,t_wvg,WN,w,h0,n_f,engine):
    """Function used to generate a cubic tapered waveguide region to be added to the right side of the cavity. Note: the new version tapers both the lattice constants and the radius prefactors. 
    
    Args:
        a: the lattice constant used to build the mirror region 
        t_wvg: the tapering prefactor associated with the waveguide region 
        WN: the number of unit cells in the waveguide region 
        d: hole diameter prefactor 
        d_min: the hole diameter prefactor we are tapering to 
        w: the beam width prefactor
        h0: the beam height
        n_f: the refractive index associated with the material
        engine: the FDTD engine used to simulate the waveguide region
    """
    w0 = w*a
    waveguide_cells_R = []
    amin = t_wvg*a
    a_wv = cubic_tapering(a,amin,WN)
    a_wv = a_wv[::-1]
    d_wv = np.linspace(d,d_min,WN)
    for i in range(WN):
        waveguide_box_R = BoxStructure(Vec3(0), Vec3(a_wv[i],w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        waveguide_hole_R = CylinderStructure(Vec3(0), h0, d_wv[i]*a_wv[i]/2, DielectricMaterial(1, order=1, color="blue"))
        waveguide_cells_R += [UnitCell(structures=[ waveguide_box_R, waveguide_hole_R ], size=Vec3(a_wv[i],w0,h0), engine=engine)]
    return waveguide_cells_R

def buildWaveguideRegion_left_v2(a,d,d_min,t_wvg,WN,w,h0,n_f,engine):
    """Function used to generate a cubic tapered waveguide region to be added to the left side of the cavity. Note: the new version tapers both the lattice constants and the radius prefactors. 
    
    Args:
        a: the lattice constant used to build the mirror region 
        t_wvg: the tapering prefactor associated with the waveguide region 
        WN: the number of unit cells in the waveguide region 
        d: hole diameter prefactor 
        d_min: the hole diameter prefactor we are tapering to 
        w: the beam width prefactor
        h0: the beam height
        n_f: the refractive index associated with the material
        engine: the FDTD engine used to simulate the waveguide region
    """
    w0 = w*a
    waveguide_cells_L = []
    amin = t_wvg*a
    a_wv = cubic_tapering(a,amin,WN)
    d_wv = np.linspace(d_min,d,WN)
    for i in range(WN):
        waveguide_box_L = BoxStructure(Vec3(0), Vec3(a_wv[i],w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        waveguide_hole_L = CylinderStructure(Vec3(0), h0, d_wv[i]*a_wv[i]/2, DielectricMaterial(1, order=1, color="blue"))
        waveguide_cells_L += [UnitCell(structures=[ waveguide_box_L, waveguide_hole_L ], size=Vec3(a_wv[i],w0,h0), engine=engine)]
    return waveguide_cells_L

def buildUnitCell_elliptical(a,hx,hy,w0,h0,n_f,engine):
    """the function use the given parameters to build a unit cell with elliptical holes

    Args:
        a (float): the lattice constant 
        hx (float): hole diameter in the x direction
        hy (float): hole diameter in the y direction
        w0 (float): the beam width 
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        engine : the FDTD engine used to simulate the waveguide region

    Returns:
        cell: the UnitCell object in waveguide solver 
    """
    cell_box = BoxStructure(Vec3(0), Vec3(a,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
    hole = CylinderStructure(Vec3(0), h0, hx, DielectricMaterial(1, order=1, color="blue"),radius2=hy)
    cell = UnitCell(structures=[ cell_box, hole ], size=Vec3(a,w0,h0), engine=engine)
    return cell

def buildMirrorRegion_elliptical(a,hx,hy,MN,w0,h0,n_f,engine):
    """the function used to build mirriro region with elliptical holes

    Args:
        a (float): the lattice constant used to build the mirror region 
        hx (float): hole diameter in the x direction
        hy (float): hole diameter in the y direction
        w0 (float): the beam width 
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        MN (int): the number of mirror unit cells
        engine: the FDTD engine used to simulate the waveguide region
    """
    mirror_cell = buildUnitCell_elliptical(a,hx,hy,w0,h0,n_f,engine)
    mirror_cells = [mirror_cell] * MN
    return mirror_cells

def buildTaperRegion_elliptical(a_L,a_R,amin,hx,hy,TN,w0,h0,n_f,engine):
    """the function used to build taper region with elliptical holes

    Args:
        a_L (float): the lattice constant used to build the mirror region on the left side 
        a_R (float): the lattice constant used to build the mirror region on the right side 
        hx (float): hole diameter in the x direction
        hy (float): hole diameter in the y direction
        w0 (float): the beam width 
        amin: the minimum lattice constant in the tapering region
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        TN (int): the number of tapering cells 
        engine: the FDTD engine used to simulate the waveguide region
    """
    taper_cells = []
    aList_taper = buildTapering_asymmetric(a_L,a_R,amin,TN)
    for i in aList_taper:
        temp_cell = buildUnitCell_elliptical(i,hx,hy,w0,h0,n_f,engine)
        taper_cells += [temp_cell]
    return taper_cells

def buildTaperRegion_elliptical_asymmetric(a_L,a_R,amin,hx,hy,TN_L,TN_R,w0,h0,n_f,engine):
    """the function used to build taper region with elliptical holes

    Args:
        a_L (float): the lattice constant used to build the mirror region on the left side 
        a_R (float): the lattice constant used to build the mirror region on the right side 
        hx (float): hole diameter in the x direction
        hy (float): hole diameter in the y direction
        w0 (float): the beam width 
        amin: the minimum lattice constant in the tapering region
        h0 (float): the beam height
        n_f (float): the refractive index associated with the material
        TN_L (int): the number of tapering cells on the left side of the tapering region 
        TN_R (int): the number of tapering cells on the right side of the tapering region
        engine: the FDTD engine used to simulate the waveguide region
    """
    taper_cells = []
    aList_taper = buildTapering_asymmetric_v2(a_L,a_R,amin,TN_L,TN_R)
    for i in aList_taper:
        temp_cell = buildUnitCell_elliptical(i,hx,hy,w0,h0,n_f,engine)
        taper_cells += [temp_cell]
    return taper_cells


def unitCellOptimization_SiC_elliptical(params):
    """this function is used to optimize the unit cells for the waveguide regions of the cavity

    Args:
        params (list): 
            params[0] (float): the tapered lattice constant (this should be the value the waveguide is tapering to)
            params[1] (float): d1 the radius prefactor 1 
            params[2] (float): d2 the radius prefactor 2

    Returns:
        fitness: the optimization parameter (we want large bandgap and small detuning)
    """
    a = params[0]
    hx = params[1]
    hy = params[2]
    w0 = params[3]
    h0 = params[4]
    # simulate the band gap of the unit cell 
    diel_freq, air_freq, mg, bg_mg_rat, delta_k,bg = sim_bandGap_elliptical(a,hx,hy,w0,h0)
    detuning = np.abs((3e8)/target_frequency - (3e8)/mg)
    detuning_nm = detuning*1e9
    print("Detuning from the mid band: %f nm"%(detuning_nm))
    # we want large bandgap and small detuning 
    delta_wv = 1e-9
    fitness = np.exp(-(detuning/delta_wv)**2)*bg_mg_rat
    file_name = "unitcell_Optimization_elliptical_500nm_testRun_1.csv"
    data = [a,hx,hy,w0,h0,detuning,mg,bg_mg_rat,bg,diel_freq, air_freq,fitness]
    record_data(data,file_name)
    return -1*fitness

def unitCellOptimization_SiC_waveguide_elliptical(params):
    """this function is used to optimize the unit cells for the waveguide regions of the cavity

    Args:
        params (list): 
            params[0] (float): the tapered lattice constant (this should be the value the waveguide is tapering to)
            params[1] (float): the radius in the x direction
            params[2] (float): the radius in the y direction

    Returns:
        fitness: the optimization parameter (we want large bandgap and small detuning)
    """
    a = params[0]
    hx = params[1]
    hy = params[2]
    # simulate the band gap of the unit cell 
    diel_freq, air_freq, mg, bg_mg_rat, delta_k, bg = sim_bandGap_elliptical(a,hx,hy)
    detuning = np.abs(target_frequency-diel_freq)
    print("Detuning from the dielectric band: %f"%(detuning))
    file_name = "unitcell_Optimization_elliptical_waveguide_v1.csv"
    data = [a,detuning]
    record_data(data,file_name)
    return detuning

def band_structure_elliptical(a,hx,hy,w0,h0,n_f,engine):
    """This function simulate the band structure of the unit cells with elliptical holes 

    Args:
        a (float): the lattice constant 
        hx (float): the radius in the x direction 
        hy (float): the radius in the y direction
        w0 (float): the beam width 
        h0 (float, optional): the beam height. Defaults to h0.
        n_f (float, optional): the refractive index. Defaults to n_f.
        engine (Lumerical engine object, optional): _description_. Defaults to engine.
    """
    print("Starting simulation =============================")
    start_time = datetime.now()
    cell = buildUnitCell_elliptical(a,hx,hy,w0,h0,n_f,engine)
    # r1 = cell.simulate("bandstructure", ks=(0.2, 0.5, 12), freqs=(0.25e15, 0.7e15, 100000),
    #                    dipole_region=Vec3(0.8, 0, 0), window_pos = 0)
    r1 = cell.simulate("bandstructure", ks=(0.2, 0.5, 8), freqs=(0.15e15, 0.6e15, 100000))
    # # # Plot the bandstructure
    r1.show()
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

def buildWaveguideRegion_elliptical_right(a,hx,hy,t_wvg,WN,w0,h0,n_f,engine):
    """Function used to generate a cubic tapered waveguide region to be added to the right side of the cavity
    
    Args:
        a: the lattice constant used to build the mirror region 
        t_wvg: the tapering prefactor associated with the waveguide region 
        WN: the number of unit cells in the waveguide region 
        n_f: the refractive index associated with the material
        hx: hole diameter in the x direction
        hy: hole diameter in the y direction
        w0: the beam width 
        h0: the beam height
        engine: the FDTD engine used to simulate the waveguide region
    """
    waveguide_cells_R = []
    amin = t_wvg*a
    a_wv = cubic_tapering(a,amin,WN)
    a_wv = a_wv[::-1]
    print(a_wv) # debugging
    for i in a_wv:
        waveguide_box_R = BoxStructure(Vec3(0), Vec3(i,w0,h0), DielectricMaterial(n_f, order=2, color="red"))
        waveguide_hole_R = CylinderStructure(Vec3(0), h0, hx, DielectricMaterial(1, order=1, color="blue"),radius2=hy)
        waveguide_cells_R += [UnitCell(structures=[ waveguide_box_R, waveguide_hole_R ], size=Vec3(i,w0,h0), engine=engine)]
    return waveguide_cells_R

def buildWaveguideRegion_elliptical_right_v2(a,hx,hx_min,hy,hy_min,t_wvg,WN,w0,h0,n_f,engine):
    """Function used to generate a cubic tapered waveguide region to be added to the right side of the cavity. Note: the new version tapers both the lattice constants and the radius prefactors. 
    
    Args:
        a: the lattice constant used to build the mirror region 
        t_wvg: the tapering prefactor associated with the waveguide region 
        WN: the number of unit cells in the waveguide region 
        hx: hole diameter prefactor 1
        hx_min: the hole diameter prefactor 1 we are tapering to 
        hy: hole diameter prefactor 2
        hy_min: the hole diameter prefactor 2 we are tapering to 
        w: the beam width prefactor
        h0: the beam height
        n_f: the refractive index associated with the material
        engine: the FDTD engine used to simulate the waveguide region
    """
    waveguide_cells_R = []
    amin = t_wvg*a
    a_wv = cubic_tapering(a,amin,WN)
    a_wv = a_wv[::-1]
    hx_wv = np.linspace(hx,hx_min,WN)
    hy_wv = np.linspace(hy,hy_min,WN)
    for i in range(WN):
        wvg_unitcell = buildUnitCell_elliptical(a_wv[i],hx_wv[i],hy_wv[i],w0,h0,n_f,engine)
        waveguide_cells_R += [wvg_unitcell]
    return waveguide_cells_R

def buildWaveguideRegion_elliptical_left_v2(a,hx,hx_min,hy,hy_min,t_wvg,WN,w0,h0,n_f,engine):
    """Function used to generate a cubic tapered waveguide region to be added to the left side of the cavity. Note: the new version tapers both the lattice constants and the radius prefactors. 
    
    Args:
        a: the lattice constant used to build the mirror region 
        t_wvg: the tapering prefactor associated with the waveguide region 
        WN: the number of unit cells in the waveguide region 
        hx: hole diameter prefactor 1
        hx_min: the hole diameter prefactor 1 we are tapering to 
        hy: hole diameter prefactor 2
        hy_min: the hole diameter prefactor 2 we are tapering to 
        w: the beam width prefactor
        h0: the beam height
        n_f: the refractive index associated with the material
        engine: the FDTD engine used to simulate the waveguide region
    """
    waveguide_cells_L = []
    amin = t_wvg*a
    a_wv = cubic_tapering(a,amin,WN)
    hx_wv = np.linspace(hx_min,hx,WN)
    hy_wv = np.linspace(hy_min,hy,WN)
    for i in range(WN):
        wvg_unitcell = buildUnitCell_elliptical(a_wv[i],hx_wv[i],hy_wv[i],w0,h0,n_f,engine)
        waveguide_cells_L += [wvg_unitcell]
    return waveguide_cells_L

def sweep_tapering_elliptical_cavity(param):
    """This function sweep through the tapering prefactor and simulate the corresponding resonsance 

    Args:
        param (0): t the tapering prefactor associated with the defect region 
    """
    print("Start sim ==============================")
    start_time = datetime.now()
    t = param[0]
    #build the left mirror cell region 
    mirror_cells_left = buildMirrorRegion_elliptical(a,hx,hy,MN_L,w0,h0,n_f,engine)
    #build the right mirror cell region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    mirror_cells_right = buildMirrorRegion_elliptical(a_R,hx,hy,MN_R,w0,h0,n_f,engine)

    #building cubic tapered cell region
    amin = a*t
    taper_cells = buildTaperRegion_elliptical(a,a_R,amin,hx,hy,TN,w0,h0,n_f,engine)
    
    centerCell = MN_L+TN-1

    #add waveguide region  
    waveguide_cells = buildWaveguideRegion_elliptical_right_v2(a,hx,hxmin_wvg,hy,hymin_wvg,t_wvg,WN,w0,h0,n_f,engine)
    
    ####################################### cavity without the waveguide region ###############################
    cavity = Cavity1D(
    unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right + waveguide_cells,
    structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
    center_cell=centerCell,
    center_shift=0,
    engine=engine
    )
    
    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity_elliptical.obj")

    #define mesh size (use 12nm for accuracy, currently set to 12nm)
    # man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 12e-9, dy=None, dz=None)

    # simulating the resonance and the Q #########################################################
    r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=60e-15, 
                        analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,4,8))

    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))

    cavity = Cavity1D(load_path="cavity_testing.obj",engine=engine)
    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Vmode = r1["vmode"]
    F = r1["freq"]
    
    resonance_f = float(F) # the resonance frequency 
    resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
    detuning_wavelength = target_wavelength-resonance_wavelength
    detuning_wavelength_nm = detuning_wavelength*1e9
    delta_wavelength = 5e-9 # 5nm tolerance 
    
    Q = 1/((1/Qsc) + (1/Qwvg))
    
    #prevent the mode volume from going to unrealistic values 
    if Q > 1000000:
        Q = 1000000
    
    if Qwvg > 1000000:
        Qwvg = 1000000
    
    if Vmode < 0.48:
        Vmode = 1e6
    
    
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f, detuning: %f nm" % ( Q, P, detuning_wavelength_nm))

    r1 = cavity.get_results("resonance")[-1]
    
    fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((detuning_wavelength/delta_wavelength)**2)))
    
    # record the data 
    data = [a,t,hx,hy,Qwvg,Qsc,Q,F,detuning_wavelength,fitness]
    file_name = "OptimizeListFull_elliptical_cavity_sweep_taperingPrefactor_v1.csv"
    record_data(data,file_name)
    
    end_time = datetime.now()
    print('fitness: %f'%(fitness))
    print('Duration: {}'.format(end_time - start_time))
    
    return -1*fitness

def sweep_beam_width_elliptical_cavity(param):
    """This function sweep through the beam width prefactor and simulate the corresponding resonsance 

    Args:
        param (0): w the beam width prefactor associated with the defect region 
    """
    print("Start sim ==============================")
    start_time = datetime.now()
    w = param[0]
    w0=w*a
    w0_nm = w0*1e9
    print("the beam width: %f nm" %(w0_nm))
    #build the left mirror cell region 
    mirror_cells_left = buildMirrorRegion_elliptical(a,hx,hy,MN_L,w0,h0,n_f,engine)
    #build the right mirror cell region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    mirror_cells_right = buildMirrorRegion_elliptical(a,hx,hy,MN_R,w0,h0,n_f,engine)

    #building cubic tapered cell region
    amin = a*t
    taper_cells = buildTaperRegion_elliptical(a,a_R,amin,hx,hy,TN,w0,h0,n_f,engine)
    
    centerCell = MN_L+TN-1

    ####################################### cavity without the waveguide region ###############################
    cavity = Cavity1D(
    unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right,
    structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
    center_cell=centerCell,
    center_shift=0,
    engine=engine
    )
    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity_elliptical.obj")

    #define mesh size (use 12nm for accuracy, currently set to 12nm)
    # man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 12e-9, dy=None, dz=None)
    man_mesh = MeshRegion(BBox(Vec3(0),Vec3(12e-6,0.7e-6,0.4e-6)), 20e-9, dy=None, dz=None)

    # simulating the resonance and the Q #########################################################
    r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=60e-15, 
                        analyze_time=600e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,8,8))

    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))

    cavity = Cavity1D(load_path="cavity_elliptical.obj",engine=engine)
    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Vmode = r1["vmode"]
    F = r1["freq"]
    
    resonance_f = float(F) # the resonance frequency 
    resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
    detuning_wavelength = target_wavelength-resonance_wavelength
    detuning_wavelength_nm = detuning_wavelength*1e9
    delta_wavelength = 5e-9 # 5nm tolerance 
    
    #prevent the mode volume from going to unrealistic values 
    if Vmode < 0.48:
        Vmode = 1e6
    
    Q = 1/((1/Qsc) + (1/Qwvg))

    if Q > 500000:
        Q = 500000
    
    if Qwvg > 500000:
        Qwvg = 500000
    
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f, detuning: %f nm" % ( Q, P, detuning_wavelength_nm))

    r1 = cavity.get_results("resonance")[-1]
    
    fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((detuning_wavelength/delta_wavelength)**2)))
    
    # record the data 
    data = [a,w0,hx,hy,Qwvg,Qsc,Q,F,detuning_wavelength,fitness]
    file_name = "OptimizeListFull_elliptical_cavity_sweep_beamWidth_v2.csv"
    record_data(data,file_name)
    
    end_time = datetime.now()
    print('fitness: %f'%(fitness))
    print('Duration: {}'.format(end_time - start_time))
    
    return -1*fitness

def sweep_beamWidth_ellipticalCavity_v2(param):
    print("Start sim ==============================")
    start_time = datetime.now()
    w0 = param[0]
    #build the left mirror cell region 
    mirror_cells_left = buildMirrorRegion_elliptical(a,hx,hy,MN_L,w0,h0,n_f,engine)

    #build the right mirror cell region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    mirror_cells_right = buildMirrorRegion_elliptical(a_R,hx,hy,MN_R,w0,h0,n_f,engine)

    #building cubic tapered cell region
    taper_cells = buildTaperRegion_elliptical(a,a_R,amin,hx,hy,TN,w0,h0,n_f,engine)
    
    ####################################### cavity without the waveguide region ###############################
    cavity = Cavity1D(
    unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right,
    structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
    center_cell=centerCell,
    center_shift=0,
    engine=engine
    )
    
    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity_elliptical.obj")

    # simulating the resonance and the Q #########################################################
    # r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
    #                     analyze_time=1000e-15,analyze_fspan=5.0e12,mesh_regions = [man_mesh], sim_size=Vec3(4,8,8))
    r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
                        analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,4,8))

    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))

    cavity = Cavity1D(load_path="cavity_testing.obj",engine=engine)
    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Vmode = r1["vmode"]
    F = r1["freq"]
    resonance_f = float(F) # the resonance frequency 
    resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
    detuning_wavelength = target_wavelength-resonance_wavelength
    detuning_wavelength_nm = detuning_wavelength*1e9
    delta_wavelength = 5e-9 # 5nm tolerance 
    
    #prevent the mode volume from going to unrealistic values 
    if Vmode < 0.48:
        Vmode = 1e6
    
    Q = 1/((1/Qsc) + (1/Qwvg))

    if Q > 500000:
        Q = 500000
    
    if Qwvg > 500000:
        Qwvg = 500000
    
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f, detuning: %f nm" % ( Q, P, detuning_wavelength_nm))

    r1 = cavity.get_results("resonance")[-1]
    
    fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((detuning_wavelength/delta_wavelength)**2)))
    
    # record the data 
    data = [a,hx,hy,t,w0,Vmode,Qwvg,Qsc,Q,F,detuning_wavelength,fitness]
    file_name = "OptimizeListFull_elliptical_cavity_sweep_beamWidth_v3.csv"
    record_data(data,file_name)
    
    end_time = datetime.now()
    print('fitness: %f'%(fitness))
    print('Duration: {}'.format(end_time - start_time))
    
    return -1*fitness

def sweep_cellNum_ellipticalCavity(param):
    """this function sweeps through the cell numbers associated with the cavity and calculate the associated Q

    Args:
        param[0] (float): the left mirror cell numbers
        param[1] (float): the tapering cell numbers

    Returns:
        _type_: _description_
    """
    print("Start sim ==============================")
    start_time = datetime.now()
    MN_R = param[0]
    TN = param[1]
    print("MN_R:", MN_R)
    print("TN:", TN)
    #build the left mirror cell region 
    mirror_cells_left = buildMirrorRegion_elliptical(a,hx,hy,MN_L,w0,h0,n_f,engine)

    #build the right mirror cell region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    mirror_cells_right = buildMirrorRegion_elliptical(a_R,hx,hy,MN_R,w0,h0,n_f,engine)

    #building cubic tapered cell region
    taper_cells = buildTaperRegion_elliptical(a,a_R,amin,hx,hy,TN,w0,h0,n_f,engine)
    
    ####################################### cavity without the waveguide region ###############################
    cavity = Cavity1D(
    unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right,
    structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
    center_cell=centerCell,
    center_shift=0,
    engine=engine
    )
    
    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity_elliptical.obj")

    #define mesh size (use 12nm for accuracy, currently set to 12nm)
    # man_mesh = MeshRegion(BBox(Vec3(0),Vec3(4e-6,0.6e-6,0.5e-6)), 12e-9, dy=None, dz=None)
    man_mesh = MeshRegion(BBox(Vec3(0),Vec3(10e-6,2e-6,2e-6)), 12e-9, dy=None, dz=None)

    # simulating the resonance and the Q #########################################################
    # r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, 
    #                     analyze_time=1000e-15,analyze_fspan=5.0e12,mesh_regions = [man_mesh], sim_size=Vec3(4,8,8))
    r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=60e-15, 
                        analyze_time=600e-15,mesh_regions = [man_mesh], sim_size=Vec3(4,4,8))

    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))

    cavity = Cavity1D(load_path="cavity_testing.obj",engine=engine)
    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Qxmin = r1["qxmin"]
    Qxmax = r1["qxmax"]
    Qy = 1/(2/r1["qymax"])
    Qz = 1/(1/r1["qzmin"] + 1/r1["qzmax"])
    
    print("Qx1: %f, Qx2: %f, Qy: %f, Qz: %f" % (
        Qxmin, Qxmax, Qy, Qz
    ))
    
    Vmode = r1["vmode"]
    F = r1["freq"]
    resonance_f = float(F) # the resonance frequency 
    resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
    detuning_wavelength = target_wavelength-resonance_wavelength
    detuning_wavelength_nm = detuning_wavelength*1e9
    delta_wavelength = 5e-9 # 5nm tolerance 
    
    #prevent the mode volume from going to unrealistic values 
    if Vmode < 0.48:
        Vmode = 1e6
    
    Q = 1/((1/Qsc) + (1/Qwvg))

    if Q > 500000:
        Q = 500000
    
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f, detuning: %f nm" % ( Q, P, detuning_wavelength_nm))

    r1 = cavity.get_results("resonance")[-1]
    
    fitness = np.sqrt((Qsc/Qwvg)*P*np.exp(-((detuning_wavelength/delta_wavelength)**2)))
    
    # record the data 
    data = [TN,MN_R,a,hx,hy,t,w0,Vmode,Qwvg,Qsc,Qxmin,Qxmax,Qy,Qz,Q,F,detuning_wavelength,fitness]
    file_name = "OptimizeListFull_elliptical_cavity_sweep_cellNum_v3.csv"
    record_data(data,file_name)
    
    end_time = datetime.now()
    print('fitness: %f'%(fitness))
    print('Duration: {}'.format(end_time - start_time))
    
    return -1*fitness


def report_results(r1,cavity_params,sim_params,file_name,file_loc):
    """
        this function take the simulation result object then report and save all of its properties

        Args:
        r1 (struct): data structure that contains all the simulation results 
    """
    # load the cavity parameters
    target_frequency = sim_params["target_frequency"]
    TN = cavity_params["TN"]
    t = cavity_params["C_lattice_tapering_prefactor"]
    w0 = cavity_params["beam_width"]
    MN_R = cavity_params["MN_Right"]
    a = cavity_params["a"]
    hx = cavity_params["hx"]
    hy = cavity_params["hy"]

    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))

    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Qxmin = r1["qxmin"]
    Qxmax = r1["qxmax"]
    Qy = 1/(2/r1["qymax"])
    Qz = 1/(1/r1["qzmin"] + 1/r1["qzmax"])
    
    print("Qx1: %f, Qx2: %f, Qy: %f, Qz: %f" % (
        Qxmin, Qxmax, Qy, Qz
    ))
    
    Vmode = r1["vmode"]
    F = r1["freq"]
    resonance_f = float(F) # the resonance frequency 
    resonance_wavelength=(3e8)/resonance_f # the resonance wavelength 
    target_wavelength = (3e8)/target_frequency # the resulting wavelength
    detuning_wavelength = target_wavelength-resonance_wavelength
    detuning_wavelength_nm = detuning_wavelength*1e9
    delta_wavelength = 5e-9 # 5nm tolerance 
    
    Q = 1/((1/Qsc) + (1/Qwvg))
    
    P = (Q*Qsc) / (Vmode*Vmode)
    print("Q: %f, P: %f, detuning: %f nm" % ( Q, P, detuning_wavelength_nm))

    # record the data 
    data = [TN,MN_R,a,hx,hy,t,w0,Vmode,Qwvg,Qsc,Qxmin,Qxmax,Qy,Qz,Q,F,detuning_wavelength,resonance_wavelength]
    record_data(data,file_name,file_loc)
    return data
    
def check_detuning(r1,source_frequency,cavity_params,sim_params):
    """
        this function takes in simulation result objects and determine if the detuning is large enough to warrant a rerun 
        
        Args:
        r1 (struct): data structure that contains all the simulation results 
    """
    freq = r1["res"]["freq"]
    wavelen_pen = np.exp(-((source_frequency - freq) / source_frequency) ** 2)
    rerun_thresh = 0.95
    if wavelen_pen < rerun_thresh:
        # shift source frequency to cavity resonance and rerun simulation.
        # (this should help avoid non cavities with artificially low mode volumes)
        print('Source frequency =' + str(source_frequency))
        print('Simmed frequency =' + str(freq))
        print('Wavelen_penalty =' + str(wavelen_pen))
        print('---n')
        sim_params['source frequency'] = freq
        witness_rerun = sim_ellipticalCavity_v2(cavity_params, sim_params)
        return witness_rerun
    else:
        return r1

def setup_engine(sim_params):
    """
        this function is used to setup the engine used to run the simulations
    """
    running_cluster = sim_params["running_cluster"] # check if we are running the simulation on the cluster 
    if running_cluster:
        FDTDloc=sim_params["FDTDloc_cluster"] # for running on the cluster 
    else:
        FDTDloc=sim_params["FDTDloc_local"] # for running on the local desktop
    engine = LumericalEngine(mesh_accuracy=sim_params["mesh_accuracy"], hide=sim_params["hide_GUI"], lumerical_path=FDTDloc, save_fsp=sim_params["save_fsps"])
    man_mesh = MeshRegion(BBox(Vec3(0),sim_params["mesh_box"]), sim_params["mesh_res"], dy=None, dz=None)
    return engine, man_mesh

def sim_ellipticalCavity_v2(cavity_params,sim_params):
    """
    this function sweeps through the hy associated with the cavity and calculate the associated Q

    Args:
        cavity_params: data struct containing all the parameters needed to build the cavity
        sim_params: data struct containing all the parameters needed to run the simulation 
    Returns:
        _type_: none
    """
    
    # read the relevant parameter data 
    a = cavity_params["a"]
    hx = cavity_params["hx"]
    hy = cavity_params["hy"]
    MN_L = cavity_params["MN_Left"]
    MN_R = cavity_params["MN_Right"]
    WN = cavity_params["WN"]
    t_wvg = cavity_params["WG_lattice_tapering_prefactor"]
    w0 = cavity_params["beam_width"]
    h0 = cavity_params["thickness"]
    l = cavity_params["cavity_length"]
    n_f = cavity_params["n_refractive"]
    d_min = cavity_params["WG_hole_tapering_prefactor"]
    t = cavity_params["C_lattice_tapering_prefactor"]
    TN = cavity_params["TN"]
    prefactor_mirror_R = cavity_params["M_lattice_prefactor"]
    engine, man_mesh = setup_engine(sim_params)
    hxmin_wvg = d_min*hx
    hymin_wvg = d_min*hy
    # default location of the data files 
    file_loc = sim_params["simulationData_loc"]
    file_name = sim_params["simulationData_fileName"]
    # The target resonance frequency, in Hz
    # 916nm = 327.3e12 Hz
    target_frequency = sim_params["target_frequency"]

    print("Start sim ==============================")
    centerCell = MN_L+TN-1 
    start_time = datetime.now()
    #build the left mirror cell region 
    mirror_cells_left = buildMirrorRegion_elliptical(a,hx,hy,MN_L,w0,h0,n_f,engine)

    #build the right mirror cell region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    amin = t*a # the defect lattice constant 
    mirror_cells_right = buildMirrorRegion_elliptical(a_R,hx,hy,MN_R,w0,h0,n_f,engine)

    #building cubic tapered cell region
    taper_cells = buildTaperRegion_elliptical(a,a_R,amin,hx,hy,TN,w0,h0,n_f,engine)
    
    #add waveguide region 
    waveguide_cells = buildWaveguideRegion_elliptical_right_v2(a,hx,hxmin_wvg,hy,hymin_wvg,t_wvg,WN,w0,h0,n_f,engine)
    
    ####################################### cavity without the waveguide region ###############################
    cavity = Cavity1D(
    unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right + waveguide_cells,
    structures=[ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ],
    center_cell=centerCell,
    center_shift=0,
    engine=engine,
    boundaries=sim_params["boundary_condition"]
    )

    # By setting the save path here, the cavity will save itself after each simulation to this file
    cavity.save("cavity.obj")

    # simulating the resonance and the Q 
    r1 = cavity.simulate("resonance", target_freq=target_frequency, source_pulselength=200e-15, analyze_time=1000e-15,mesh_regions = [man_mesh], sim_size=Vec3(1.5,3,8))
    
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))

    # Print the reults and plot the electric field profiles
    print("F: %f, Vmode: %f, Qwvg: %f, Qsc: %f" % (
        r1["freq"], r1["vmode"],
        1/(1/r1["qxmin"] + 1/r1["qxmax"]),
        1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    ))
    if sim_params["show_field_profile"]:
        r1["xyprofile"].show()
        r1["yzprofile"].show()

    cavity = Cavity1D(load_path="cavity.obj",engine=engine)
    r1 = cavity.get_results("resonance")[-1]

    # check detuning (prevent unrealistically small mode volumn)
    r1 = check_detuning(r1,target_frequency,cavity_params,sim_params)
    # report and save the results 
    report_results(r1['res'],cavity_params,sim_params,file_name,file_loc)

    return r1

def calculate_fitness(r1,sim_params):
    """this function calculates the fitness associated with the simulation result

    Args:
        r1 (resonance dict): dict containing all the relevant cavity simulation results
    """
    # calculate the fitness 
    target_frequency = sim_params["target_frequency"]
    Vmode = r1["vmode"]
    F = r1["freq"]
    Qwvg = 1/(1/r1["qxmin"] + 1/r1["qxmax"])
    Qsc = 1/(2/r1["qymax"] + 1/r1["qzmin"] + 1/r1["qzmax"])
    Qxmin = r1["qxmin"]
    Qxmax = r1["qxmax"]
    Q = 1/((1/Qsc) + (1/Qwvg))
    delta_wavelength = 1e-9 # 1nm tolerance 
    resonance_wavelength=(3e8)/F
    target_wavelength = (3e8)/target_frequency 
    detuning_wavelength = np.abs(resonance_wavelength-target_wavelength)
    # account for unrealistic mode volumes 
    if Vmode < 0.4:
        Vmode = 1e6
    fitness = np.exp(-((detuning_wavelength/delta_wavelength)**2))*(Q*Qsc) / (Vmode*Vmode)
