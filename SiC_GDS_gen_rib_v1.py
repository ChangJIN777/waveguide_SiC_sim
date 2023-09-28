"""
This script is used to generate gds files from the existing designs 
"""
    
from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine
from wvgsolver.parse import DielectricExtrusionFaceGDSParser

import scipy.optimize
import numpy as np
import os
from datetime import datetime

from waveguideSolver_funcs import *
from cavity_sim_parameters import *

def build_rib_cavity_v1(rib_cavity_params,rib_sim_params):
    # read the relevant parameter data 
    a = rib_cavity_params["a"]
    hx = rib_cavity_params["hx"]
    hy = rib_cavity_params["hy"]
    MN_L = rib_cavity_params["MN_Left"]
    MN_R = rib_cavity_params["MN_Right"]
    WN = rib_cavity_params["WN"]
    t_wvg = rib_cavity_params["WG_lattice_tapering_prefactor"]
    w0 = rib_cavity_params["beam_width"]
    h0 = rib_cavity_params["thickness"]
    l = rib_cavity_params["cavity_length"]
    n_f = rib_cavity_params["n_refractive"]
    d_min = rib_cavity_params["WG_hole_tapering_prefactor"]
    t = rib_cavity_params["C_lattice_tapering_prefactor"]
    TN = rib_cavity_params["TN"]
    prefactor_mirror_R = rib_cavity_params["M_lattice_prefactor"]
    do_sc = rib_cavity_params["do_sc"]
    sc_gap = rib_cavity_params["sc_gap"]
    engine, man_mesh = setup_engine(rib_sim_params)
    hxmin_wvg = d_min*hx
    hymin_wvg = d_min*hy
    # default location of the data files 
    file_loc = rib_sim_params["simulationData_loc"]
    file_name = rib_sim_params["simulationData_fileName"]
    # The target resonance frequency, in Hz
    # 916nm = 327.3e12 Hz
    target_frequency = rib_sim_params["target_frequency"]
    mode_orientation = rib_sim_params["mode"]
    freq_span = rib_sim_params["freq_span"]
    # used for side coupling
    sc_pos = Vec3(0, sc_gap + w0, 0)
    wg_size = Vec3(a, w0, h0)
    sc_cell_box = BoxStructure(sc_pos, wg_size,
                                DielectricMaterial(n_f, order=2, color="blue"))

    print("Start sim ==============================")
    centerCell = MN_L+TN
    start_time = datetime.now()
    #build the left mirror cell region 
    mirror_cells_left = buildMirrorRegion_rib(a,hx,hy,MN_L,w0,h0,n_f,do_sc,sc_gap,sc_cell_box,engine)

    #build the right mirror cell region 
    a_R = a*prefactor_mirror_R # the lattice constant associated with the right mirror region 
    amin = t*a # the defect lattice constant 
    mirror_cells_right = buildMirrorRegion_rib(a_R,hx,hy,MN_R,w0,h0,n_f,do_sc,sc_gap,sc_cell_box,engine)

    #building cubic tapered cell region
    taper_cells = buildTaperRegion_rib(a,a_R,amin,hx,hy,TN,w0,h0,n_f,do_sc,sc_gap,sc_cell_box,engine)
    
    # #add waveguide region 
    # waveguide_cells = buildWaveguideRegion_elliptical_right_v2(a,hx,hxmin_wvg,hy,hymin_wvg,t_wvg,WN,w0,h0,n_f,engine)
    
    # whether we are simulating for the TE or the TM mode 
    if mode_orientation == 'TM':
        E_component = 'Ez'
    else:
        E_component = 'Ey'
    
    # if we are using sidecoupling
    if do_sc:
        structs = [BoxStructure(Vec3(0), Vec3(l, w0, h0),
                                DielectricMaterial(n_f, order=2, color="red")),
                   BoxStructure(sc_pos, Vec3(l, w0, h0),
                                DielectricMaterial(n_f, order=2, color="red"))
                   ]
    else:
        structs = [ BoxStructure(Vec3(0), Vec3(l, w0, h0), DielectricMaterial(n_f, order=2, color="red")) ]
    
    cavity = Cavity1D(
    unit_cells=  mirror_cells_left + taper_cells + mirror_cells_right,
    structures= structs,
    center_cell=centerCell,
    center_shift=0,
    engine=engine,
    boundaries=rib_sim_params["boundary_condition"],
    component=E_component # added for TM mode simulation
    )
    
    return cavity

rib_cavity_params = cavity_sim_parameters.rib_cavity_params
rib_sim_params = cavity_sim_parameters.rib_sim_params
rib_sim_params["running_cluster"] = False  
rib_sim_params["running_local"] = True

# generate GDS from the design (dosage test)
cavity_temp = build_rib_cavity_v1(rib_cavity_params,rib_sim_params)
parser = DielectricExtrusionFaceGDSParser(cavity_temp)
parser.show()
file_name = "SiC_rib_cavity_v1_test.gds"
parser.save(file_name)