import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy import units as u


#This function reads .out files, and extracts collision times, collision types, and number of bodies versus time
def extract_data_outfile(file_path):
    """
    Extracts collision time, types, and t versus number of bodies in simulation from the output file.
    Returns two pandas dataframes:
    1. t_type: with columns ['coll_times', 'coll_types']
    2. t_n: with columns ['t', 'n_bodies']
    """
    coll_time = []
    coll_types = []
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
                
    # Convert lists to numpy arrays for structured data
    n_bodies = np.array(n_bodies_list, dtype=[('t', 'i4'), ('n_bodies', 'i4')])
    
    # Create pandas DataFrame for collision time and types
    t_type = pd.DataFrame({
        'coll_times': coll_time,
        'coll_types': coll_types
    })
    
    # Create pandas DataFrame for n_bodies (t and n_bodies)
    t_n = pd.DataFrame(n_bodies)
    
    return t_type, t_n



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

#function to plot v versus v
def plot_b_v(b, v, type):
    shape = np.zeros(len(type), dtype=object)
    color = np.zeros(len(type), dtype=object)
    
    for i in range(len(type)):
        if type[i] == 'EFFECTIVELY MERGED':
            shape[i] = 'D'
            color[i] = 'blue'
        
        elif type[i] == 'SIMPLY MERGED':
            shape[i] = 'D'
            color[i] = 'red'
        
        elif type[i] == 'PARTIAL ACCRETION':
            shape[i] = 'o'
            color[i] = 'black'
        
        elif type[i] == 'PARTIAL EROSION':
            shape[i] = '^'
            color[i] = 'grey'
    
        elif type[i] == 'SUPER-CATASTROPHIC':
            shape[i] = '^'
            color[i] = 'magenta'
            
        elif type[i] == 'GRAZE AND MERGE':
            shape[i] = 'o'
            color[i] = 'green'
        
        elif type[i] == 'ELASTIC BOUNCE':
            shape[i] = '*'
            color[i] = 'brown'
            
        elif type[i] == 'HIT AND RUN':
            shape[i] = '^'
            color[i] = 'orange'
            
        else:
            shape[i] = '*'
            color[i] = 'cyan'
            
    plt.figure(figsize=(10,6))

    # Track the types already used in the legend
    used_labels = set()

    for i in range(len(type)):
        label = type[i] if type[i] not in used_labels else None  # Only add label if it hasn't been used
        if label:
            used_labels.add(type[i])  # Add the label to the set

        plt.scatter(
        b[i], v[i], 
        label=label,
        marker=shape[i], 
        color=color[i]
        )

    plt.yscale('log')
    plt.xlabel('b/R_target')
    plt.ylabel('v/v_esc')
    plt.legend(loc='upper left')  # Show legend if needed

    plt.show()
    
    
#Function to make CDF
def make_cdf(data):
    data_sorted = np.sort(data)
    cdf = np.arange(1, len(data_sorted) + 1) / len(data_sorted)
    return data_sorted, cdf
    

def plot_v_cdf(v,color, label):
    v_plot, cdf = make_cdf(v)
    
    plt.plot(np.log(v_plot), cdf, color=color, label=label)

    plt.xlabel('log(vi_vesc)')
    plt.ylabel('CDF')
    plt.title('CDF of vi/vesc')
    plt.legend()
    plt.grid(True)
    plt.show()