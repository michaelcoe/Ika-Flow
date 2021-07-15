import os
import sys
import re
import numpy as np
import pandas as pd

import processingIO as pio

from pathlib import Path

data_location = r'/home/mco143/Documents/Current_OpenFOAM_Simulations/gridInd_Ebrahim_snappyHexMesh/Ebrahim_validation3'
save_location = r'/home/mco143/Documents/Current_OpenFOAM_Simulations/gridInd_Ebrahim_snappyHexMesh/Ebrahim_validation3'

'''
Example path to the postProcessing folder
/path/to/simulation/files/backGround/postProcessing/---/0/data

--- is:

forceCoeffs_object --> forceCoeff

forces_object --> forces

solverInfo --> solverInfo

yplus --> yPlus
'''

force_paths = pio.get_files(data_location, 'force.dat')
forceCoeff_paths = pio.get_files(data_location, 'coefficient.dat')
yPlus_paths = pio.get_files(data_location, 'yPlus.dat')
solverInfo_paths = pio.get_files(data_location, 'solverInfo.dat')

force_paths.sort()
force_sheets = [x.parts[-6] for x in force_paths]

forceCoeff_paths.sort()
forceCoeff_sheets = [x.parts[-6] for x in forceCoeff_paths]

yPlus_paths.sort()
yPlus_sheets = [x.parts[-6] for x in yPlus_paths]

solverInfo_paths.sort()
solverInfo_sheets = [x.parts[-6] for x in solverInfo_paths]

startTime = 0.0
endTime = 20.0

forceCoeff_save = pio.write_excel_files(save_location, forceCoeff_paths, 
                                        forceCoeff_sheets, 'forceCoeff', startTime, endTime, None)
force_save = pio.write_excel_files(save_location, force_paths, 
                                        force_sheets, 'forces', startTime, endTime, None)
yplus_save = pio.write_excel_files(save_location, yPlus_paths, 
                                        yPlus_sheets, 'yPlus', startTime, endTime, 'wing')
solverInfo_save = pio.write_excel_files(save_location, solverInfo_paths, 
                                        solverInfo_sheets, 'solverInfo', startTime, endTime, None)
