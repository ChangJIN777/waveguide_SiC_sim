
from wvgsolver import Cavity1D, UnitCell, Vec3, Waveguide
from wvgsolver.geometry import BoxStructure, CylinderStructure, DielectricMaterial, MeshRegion
from wvgsolver.utils import BBox
from wvgsolver.engine import LumericalEngine

import scipy.optimize
import numpy as np
import os
from datetime import datetime
import csv

class cavity_sim_parameters:

    cavity_params = {
        # unit cell dimensions
        "a": 2.655634610941198e-7, #lattice constant 
        "hx": 9.347181843644768e-08,
        "hy": 1.191433948602531e-07,
        "thickness": 500e-9,
        "beam_width": 5.050493279419303e-07,
        # material parameters
        "n_refractive": 2.6,
        # device parameters 
        "MN_Left": 10, # the number of mirror cells on the left side of the cavity 
        "MN_Right": 6, # the number of mirror cells on the right side of the cavity
        "TN": 6, # symmetric cavity: the number of unit cells in the tapering region 
        "WN": 3, # the number of unit cells in the waveguide region 
        "WG_hole_tapering_prefactor": 0.437, # For the waveguide region: the prefactor that defines the minimum hole dimensions we are tapering to
        "C_lattice_tapering_prefactor": 0.852, # For the cavity region: the prefactor that defines the minimum lattice constant in the cavity region 
        "WG_lattice_tapering_prefactor": 0.852, # For the waveguide region: the prefactor that defines the minimum lattice constant in the waveguide region 
        "cavity_length": 20e-6, # the length of the device (most be larger than the simulation size)
        "M_lattice_prefactor": 0.965, # For the weaker mirror region: the prefactor that specifies the lattice constants of the weaker mirror region 
        # the parameters used for doing sidecouplings
        "do_sc": False, # Boolean used to determine if we are doing side coupling or not 
        "sc_gap": 100e-9, # the gap between the side coupling waveguide and the caivty
    }
        
    sim_params = {
        # simulation parameters 
        "target_frequency": 327.3e12, # the target resonance frequency in Hz (specify the dipole source frequency)
        "FDTDloc_local": 'C:/Program Files/Lumerical/v221/', # the location of the FDTD engine on local PC
        "FDTDloc_cluster": "/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/", # the location of the FDTD engine on the cluster
        "hide_GUI": False, # specify if we are going to hide the guide 
        "save_fsps": True, # specify if we are going to save the fsps file 
        "simulationData_loc": "./sim_data/", # the folder we are going to store the simulation data files in 
        "simulationData_fileName": "SiC_500nm_testRun_t1.txt",
        "mesh_res": 15e-9, # specify the resolution of the meshing we are using 
        "mesh_box": Vec3(4e-6,2e-6,2e-6), # specify the size of the fine meshing box (12nm for accuracy)
        "running_cluster": False, # specify if we are running the code on the cluster 
        "running_local": True, # specify if we are running on local PC
        "mesh_accuracy": 4, # specify the mesh accuracy of the engine 
        "boundary_condition": ['ymin','zmin'], # specify the symmetry to save time on simulation
        "show_field_profile": False
    }
    
    rib_cavity_params = {
        # unit cell dimensions
        "a": 2.655634610941198e-7, #lattice constant 
        "hx": 9.347181843644768e-08,
        "hy": 1.191433948602531e-07,
        "thickness": 500e-9,
        "beam_width": 5.050493279419303e-07,
        "nPoints": 20, # specify the number of points used to make up the gaussian curve
        # material parameters
        "n_refractive": 2.6,
        # device parameters 
        "MN_Left": 10, # the number of mirror cells on the left side of the cavity 
        "MN_Right": 6, # the number of mirror cells on the right side of the cavity
        "TN": 6, # symmetric cavity: the number of unit cells in the tapering region 
        "WN": 3, # the number of unit cells in the waveguide region 
        "WG_hole_tapering_prefactor": 0.437, # For the waveguide region: the prefactor that defines the minimum hole dimensions we are tapering to
        "C_lattice_tapering_prefactor": 0.852, # For the cavity region: the prefactor that defines the minimum lattice constant in the cavity region 
        "WG_lattice_tapering_prefactor": 0.852, # For the waveguide region: the prefactor that defines the minimum lattice constant in the waveguide region 
        "cavity_length": 20e-6, # the length of the device (most be larger than the simulation size)
        "M_lattice_prefactor": 0.965, # For the weaker mirror region: the prefactor that specifies the lattice constants of the weaker mirror region 
        "shift": -1 # DONT KNOW WHAT THIS IS YET <-- ask Michael 
    }
        
    rib_sim_params = {
        # simulation parameters 
        "target_frequency": 327.3e12, # the target resonance frequency in Hz (specify the dipole source frequency)
        "FDTDloc_local": 'C:/Program Files/Lumerical/v221/', # the location of the FDTD engine on local PC
        "FDTDloc_cluster": "/n/sw/lumerical-2021-R2-2717-7bf43e7149_seas/", # the location of the FDTD engine on the cluster
        "hide_GUI": True, # specify if we are going to hide the guide 
        "save_fsps": False, # specify if we are going to save the fsps file 
        "simulationData_loc": "./sim_data/", # the folder we are going to store the simulation data files in 
        "simulationData_fileName": "SiC_500nm_testRun_t1.txt",
        "mesh_res": 15e-9, # specify the resolution of the meshing we are using 
        "mesh_box": Vec3(4e-6,2e-6,2e-6), # specify the size of the fine meshing box (12nm for accuracy)
        "running_cluster": False, # specify if we are running the code on the cluster 
        "running_local": True, # specify if we are running on local PC
        "mesh_accuracy": 4, # specify the mesh accuracy of the engine 
        "boundary_condition": ['ymin','zmin'], # specify the symmetry to save time on simulation
        "show_field_profile": False
    }
