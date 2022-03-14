import os
import sys
import re
import numpy as np
import pandas as pd

from pathlib import Path

def get_files(file_path, extension):
    
    paths = []
    for path in Path(file_path).rglob(extension):
        paths.append(path)
    
    return paths

'''
Retrieve the data with an array of arrays for each parameter that's important.

Format of the data:

# Force Data
Time (total_x total_y total_z)	(pressure_x pressure_y pressure_z)	(viscous_x viscous_y viscous_z)

# Force Coefficient Data
Time    Cd  Cs  Cl  CmRoll  CmPitch CmYaw   Cd(f)   Cd(r)   Cs(f)   Cs(r)   Cl(f)   Cl(r)

# Solver Info (Residuals)
Time    U_solver    Ux_initial  Ux_final    Ux_iters    Uy_initial  Uy_final    Uy_iters    U_converged k_solver    k_initial   k_final        	k_iters           	k_converged       	p_solver          	p_initial         	p_final           	p_iters           	p_converged       	omega_solver      	omega_initial     	omega_final       	omega_iters       	omega_converged

# Y Plus
Time    patch   min max average
'''
# function to process force.dat files
def process_force_file(file_name, startTime = 0.0, endTime = 20.0):
    raw = []

    with open(file_name, 'r') as f:
        for line in f:
            tmp = [x.strip('(').strip(')') for x in line.split()]
            if len(tmp) == 0:
                continue
            elif tmp[0] == '#':
                continue
            else:
                if float(tmp[0]) < startTime:
                    continue
                elif float(tmp[0]) > endTime:
                    break
                else:
                    try:
                        raw.append([ float(i) for i in tmp ])
                    except:
                        print("could not convert string to float in line:")
                        print("\t" + line)
                        print("in file:")
                        print("\t" + file_name)

    raw = np.array(raw)

    return np.array([raw[:,0], raw[:,1], raw[:,2], raw[:,3], raw[:,4], raw[:,5], raw[:,6], raw[:,7], raw[:,8], raw[:,9]])

# function to process coefficient.dat files
def process_forceCoeff_file(file_name, startTime = 0.0, endTime = 20.0):
    data = np.loadtxt(file_name, comments='#', skiprows=13)
    mask = np.logical_and(data[:,0] > startTime, data[:,0] < endTime)
    raw = data[mask]
    return np.array([raw[:,0], raw[:,1], raw[:,2], raw[:,3], raw[:,4], raw[:,5], raw[:,6], raw[:,7], raw[:,8], raw[:,9], 
                    raw[:,10], raw[:,11], raw[:,12]])

# function to process solverInfo.dat files
def process_solverInfo_file(file_name, startTime = 0.0, endTime = 20.0):
    data = np.loadtxt(file_name, comments='#', skiprows=2, usecols=(0, 2, 3, 5, 6, 10, 11, 15, 16, 20, 21))
    mask = np.logical_and(data[:,0] > startTime, data[:,0] < endTime)
    raw = data[mask]
    return np.array([raw[:,0], raw[:,1], raw[:,2], raw[:,3], raw[:,4], raw[:,5], raw[:,6], raw[:,7], raw[:,8], raw[:,9], 
                    raw[:,10]])

# function to process yplus.dat files
def process_yPlus_file(file_name, patch, startTime = 0.0, endTime = 20.0):
    raw = []

    with open(file_name, 'r') as f:
        for line in f:
            tmp = line.split()
            if len(tmp) == 0:
                continue
            elif tmp[0] == '#':
                continue
            else:
                if float(tmp[0]) < startTime:
                    continue
                elif float(tmp[0]) > endTime:
                    break
                elif tmp[1] != patch:
                    continue
                else:
                    try:
                        matches = re.findall("[+-]?\d+\.\d+|\d+", line)
                        raw.append([ float(i) for i in matches ])
                    except:
                        print("could not convert string to float in line:")
                        print("\t" + line)
                        print("in file:")
                        print("\t" + file_name)

    raw = np.array(raw)

    return np.array([raw[:,0], raw[:,1], raw[:,2], raw[:,3]])


# write excel files based on object type
def write_excel_files(save_location, files, sheet_names, object_type, startTime=0.0, endTime=20.0, patch_name=None):
    
    save_path = Path(save_location).joinpath(object_type)
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    writer = pd.ExcelWriter(Path(save_path).joinpath(object_type + '.xlsx'))

    # Save Force Data
    if object_type == 'forces':
        for i, path in enumerate(files):
            sheet = sheet_names[i]
            
            [time, total_x, total_y, total_z, pressure_x, pressure_y, pressure_z, viscous_x, viscous_y, viscous_z] = process_force_file(path, 
            startTime, endTime)
            
            data_template = {'Time':time, 'total_x':total_x, 'total_y':total_y, 'total_z':total_z, 'pressure_x':pressure_x, 
            'pressure_y':pressure_y, 'pressure_z':pressure_z, 'viscous_x':viscous_x, 'viscous_y':viscous_y, 'viscous_z':viscous_z}
            
            mesh_df = pd.DataFrame(data_template)
            mesh_df.to_excel(writer, sheet_name=sheet, index=False)

    # Save Force Coefficient Data
    elif object_type == 'forceCoeff':
        for i, path in enumerate(files):
            sheet = sheet_names[i]
            
            [time, Cd, Cs, Cl, CmRoll, CmPitch, CmYaw, Cdf, Cdr, Csf, Csr, Clf, Clr] = process_forceCoeff_file(path, startTime, endTime)
            
            data_template = {'Time':time, 'Cd':Cd, 'Cs':Cs, 'Cl':Cl, 'CmRoll':CmRoll, 'CmPitch':CmPitch, 'CmYaw':CmYaw, 
            'Cd(f)':Cdf, 'Cd(r)':Cdr, 'Cs(f)':Csf, 'Cs(r)':Csr, 'Cl(f)':Clf, 'Cl(r)':Clr}
            
            mesh_df = pd.DataFrame(data_template)
            mesh_df.to_excel(writer, sheet_name=sheet, index=False)

    # Save Solver Info Data
    elif object_type == 'solverInfo':
        for i, path in enumerate(files):
            sheet = sheet_names[i]
            
            [time, Ux_initial, Ux_final, Uy_initial, Uy_final, k_initial, k_final, p_initial, p_final, omega_initial, 
            omega_final] = process_solverInfo_file(path, startTime, endTime)
            
            data_template = {'Time':time, 'Ux_initial':Ux_initial, 'Ux_final':Ux_final, 'Uy_initial':Uy_initial, 'Uy_final':Uy_final, 
            'k_initial':k_initial, 'k_final':k_final, 'p_initial':p_initial, 'p_final':p_final, 'omega_initial':omega_initial, 
            'omega_final':omega_final}
            
            mesh_df = pd.DataFrame(data_template)
            mesh_df.to_excel(writer, sheet_name=sheet, index=False)

    # Save yPlus Data
    elif object_type == 'yPlus':
        for i, path in enumerate(files):
            sheet = sheet_names[i]
            
            [time, min_yplus, max_yplus, average_yplus] = process_yPlus_file(path, patch_name, startTime, endTime)
            
            data_template = {'Time':time, 'min_yplus':min_yplus, 'max_yplus':max_yplus, 'average_yplus':average_yplus}
            
            mesh_df = pd.DataFrame(data_template)
            mesh_df.to_excel(writer, sheet_name=sheet, index=False)
    
    else:
        print("No valid data object file specified")

    writer.save()

    return save_path

# write csv files based on object type
def write_csv_files(save_location, files, file_names, object_type, startTime=0.0, endTime=20.0, path_name=None):

    save_path = Path(save_location).joinpath(object_type)
    if not os.path.exists(save_path):
        os.makedirs(save_path)

    # Save Force Data
    if object_type == 'forces':
        for i, path in enumerate(files):
            csv_save_location = Path(save_path).joinpath(file_names[i] + '-' + object_type + '.csv')
            raw = process_force_file(path, startTime, endTime)
            np.savetxt(csv_save_location, raw, delimiter=',', 
            header='time, total_x, total_y, total_z, pressure_x, pressure_y, pressure_z, viscous_x, viscous_y, viscous_z')

    # Save Force Coefficient Data
    elif object_type == 'forceCoeff':
        for i, path in enumerate(files):
            csv_save_location = Path(save_path).joinpath(file_names[i] + '-' + object_type + '.csv')
            raw = process_forceCoeff_file(path, startTime, endTime)
            np.savetxt(csv_save_location, raw, delimiter=',', 
            header='time, Cd, Cs, Cl, CmRoll, CmPitch, CmYaw, Cdf, Cdr, Csf, Csr, Clf, Clr')

    # Save Solver Info Data
    elif object_type == 'solverInfo':
        for i, path in enumerate(files):
            csv_save_location = Path(save_path).joinpath(file_names[i] + '-' + object_type + '.csv')
            raw = process_solverInfo_file(path, startTime, endTime)
            np.savetxt(csv_save_location, raw, delimiter=',', 
            header='time, Ux_initial, Ux_final, Uy_initial, Uy_final, k_initial, k_final, p_initial, p_final, omega_initial, omega_final')

    # Save yPlus Data
    elif object_type == 'yPlus':
        for i, path in enumerate(files):
            csv_save_location = Path(save_path).joinpath(file_names[i] + '-' + object_type + '.csv')
            raw = process_yPlus_file(path, patch_name, startTime, endTime)
            np.savetxt(csv_save_location, raw, delimiter=',', header='time, min_yplus, max_yplus, average_yplus')
    else:
        print("No valid data object file specified")
    
    return save_path