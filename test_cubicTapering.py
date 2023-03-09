import numpy as np 
import matplotlib.pyplot as plt 
import sys
import os 


def cubic_defect(a,taperPrefac,taperNum):
    """
    a: the lattice constant in the mirror region
    taperNum: the number of taper cells
    taperPrefac: taper prefactor 
    """
    a_taper = np.zeros((taperNum,))
    for i in range(taperNum-1):
        a_taper[i] = a*(1-(1-taperPrefac) * (2 * ((i) / taperNum) ** 3 - 3 * ((i) / taperNum) ** 2+1))
    a_taper[taperNum-1] = a
    return a_taper

def buildTapering(a,taperPrefac,taperNum):
    """
    function for calculating the lattice constants for the tapering region of the cavity 
    Note: this is used to build SYMMETRIC taper cell region
    """
    a_taper_R = cubic_defect(a,taperPrefac,taperNum)
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
    a_taper_R = cubic_defect(a,taperPrefac,taperNum=taperNum_R)
    a_taper_L = cubic_defect(a,taperPrefac,taperNum=taperNum_L)
    a_taper_L = a_taper_L[::-1]
    tapering_region = np.concatenate((a_taper_L, a_taper_R), axis=None)
    return tapering_region

# # testing the building the tapering region of the cell #######################################
# a = 290 # nm 
# taperNum = 8
# taperPrefac = 0.5 
# taperingNum = np.linspace(0,taperNum*2,taperNum*2)
# a_taper = buildTapering(a,taperPrefac,taperNum)
# plt.plot(taperingNum,a_taper)
# plt.show()
# # debugging
# print(a_taper) 

# testing the building the tapering region of the cell #######################################
a = 290 # nm 
taperNum_L = 8
taperNum_R = 3
taperNum = taperNum_L + taperNum_R
taperPrefac = 0.84
taperingNum = np.linspace(0,taperNum,taperNum)
a_taper = buildTapering_asymmetric(a,taperPrefac,taperNum_L,taperNum_R)
plt.plot(taperingNum,a_taper)
plt.show()
# debugging
print(a_taper) 
