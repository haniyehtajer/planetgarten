import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u


#This function reads .out files, and extracts collision times, collision types, and number of bodies versus time
def extract_data_outfile(file_path):
    """
    Extracts collision time, types, and t versus number of bodies in simulation from the output file.
    """
    coll_time = []
    coll_types = []
    n_bodies = []                     
    n_bodies_list = []
    
    # Open and read the file
    with open(file_path, 'r') as file:
        for line in file:
            # Check if the line contains the desired phrase
            if "TIME OF COLLISION:" in line:
                # Extract the value after "TIME OF COLLISION:"
                time_value = float(line.split(":")[1].strip())
                coll_time.append(time_value)
            
            if "COLLISION TYPE:" in line:
                # Extract the value after "COLLISION TYPE:"
                type_value = line.split(":")[1].strip()
                coll_types.append(type_value)
                
 
            if "N_tot=" in line and "t=" in line:
                # Extract N_tot and time values
                parts = line.split()
                n_bodies = int(parts[1])  # Extract the value after "N_tot="
                t = int(float(parts[3]))  # Extract the time value and convert to integer
                # Append the (t, n_bodies) as a tuple to the list
                n_bodies_list.append((t, n_bodies))
                
    n_bodies = np.array(n_bodies_list, dtype=[('t', 'i4'), ('n_bodies', 'i4')])
                
    
    coll_time = np.array(coll_time)
    coll_types = np.array(coll_types)
    
    
    return coll_time, coll_types, n_bodies


#extract vi/vesc and b/b_crit

#This function reads .out files, and extracts collision times, collision types, and number of bodies versus time
def extract_data_impact(file_path):
    """
    Extracts collision time, types, and t versus number of bodies in simulation from the output file.
    """
    vi_vesc = []
    b_bcrit = []
    
    # Open and read the file
    with open(file_path, 'r') as file:
        for line in file:
            # Check if the line contains the desired phrase
            if "Vimp/Vesc:" in line:
                # Extract the value after "TIME OF COLLISION:"
                vi_vesc_value = float(line.split(":")[1].strip())
                vi_vesc.append(vi_vesc_value)
            
            if "b/Rtarg:" in line:
                # Extract the value after "COLLISION TYPE:"
                b_bcrit_value = float(line.split(":")[1].strip())
                b_bcrit.append(b_bcrit_value)
                

    vi_vesc = np.array(vi_vesc)
    b_bcrit = np.array(b_bcrit)
    
    
    return vi_vesc, b_bcrit

#This function can be used to extract data from collision_report files, with the
#format used as input files for Differentiated Body Composition Tracker
def extract_data_report(file_path):
    time = []
    type = []
    
    with open(file_path, "r") as file:
        for line in file:
            # Split the line into elements
            values = line.split()
            time.append(float(values[0]))
            type.append(int(values[1]))
    
    time = np.array(time)
    type = np.array(type)
    return time, type


#Function to get data from DBCT output file
def extract_dbct_report(file_path):
    mass = []
    cmf = []
    hash = []
    
    with open(file_path, "r") as file:
        for line in file:
            # Split the line into elements
            values = line.split()
            hash.append(int(values[0]))
            mass.append(float(values[1]))
            cmf.append(float(values[2]))
            
    hash = np.array(hash)
    mass = np.array(mass)
    cmf = np.array(cmf)
    return hash, mass, cmf

