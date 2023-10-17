################################################################
#           Plot NMR observables from DFT on 2D Spectra        #
#                                                              #
#            Ettore Bartalucci, v 0.1, 12.10.23 Aachen         #
################################################################

# =============================================================================
# Import the necessary modules
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re # Regex
import nmrglue as ng # NMRglue module for handling NMR spectra 
# =============================================================================

#------------------------------------------------------------------------------#
# SECTION 1: READ INPUT FILE AND PARSE BASED BY EXPERIMENT TYPE

def read_input(input_file):
    """
    Read the input file and do a keyword match
    """

    # Get the absolute path of the input file
    input_path = os.path.abspath(input_file)
    # Get the directory path of the input file
    input_dir = os.path.dirname(input_path)

    # Open the input and search for keywords:
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # one dimensional nmr experiments
    nmr_1d = []
    # two dimensional nmr experiments
    nmr_2d = []

    current_experiment = None

    for line in lines:
        line = line.strip()

        if line.startswith('$1d_nmr_experiments'):
            current_experiment = nmr_1d
        elif line.startswith('$2d_nmr_experiments'):
            current_experiment = nmr_2d
        elif current_experiment is not None:
            current_experiment.append(line)
        
        print(current_experiment)
        
    return nmr_1d, nmr_1d

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# SECTION 2: GENERATE NMR PLOTS

def make_nmr_plot(input_file, nucleus, nmr_shifts):
    """
    Generate either a 1D or 2D NMR plot based on input data, nucleus type, and NMR shifts.

    Parameters:
    - input_file (str): Input file containing keywords for plot settings.
    - nucleus (str): Nucleus type (e.g., "1H" or "13C").
    - nmr_shifts (list): List of NMR shifts from DFT calculations.

    Returns:
    - None (displays the plot).
    """

    # Define static shift ranges for various nuclei
    shiftranges = {
        "1H": list(range(12, -2, -1)),
        "7Li": list(range(12, -20, -1)),
        "11B": list(range(90, -60, -1)),
        "13C": list(range(200, 0, -1)),
        "15N": list(range(220, 0, -1)),
        "19F": list(range(150, -300, -1)),
        "29Si": list(range(80, -50, -1)),
        "31P": list(range(200, -200, -1)),
    }

    # Open the input file and search for keywords:
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Check if 1D or 2D NMR experiment
    plot_1d_nmr = None
    plot_2d_nmr = None

    for line in lines:
        line = line.strip()

        if "$1d_nmr_experiments" in line:
            plot_1d_nmr = line.split()[-1]

        elif "$2d_nmr_experiments" in line:
            plot_2d_nmr = line.split()[-1]

    if plot_1d_nmr:
        if plot_1d_nmr in shiftranges:
            title = f"{nucleus} 1D NMR Plot"
            x_axis_length = shiftranges[plot_1d_nmr]
            create_1d_nmr_plot(nmr_shifts, title, x_axis_length)
        else:
            print(f"Unsupported 1D plot option: {plot_1d_nmr}")
    
    elif plot_2d_nmr:
        if plot_2d_nmr == "1H-1H":
            title = "1H-1H 2D NMR Plot"
            x_axis_length = shiftranges["1H"]
            create_2d_nmr_plot(nmr_shifts, title, x_axis_length)
        elif plot_2d_nmr == "13C-13C":
            title = "13C-13C 2D NMR Plot"
            x_axis_length = shiftranges["13C"]
            create_2d_nmr_plot(nmr_shifts, title, x_axis_length)
        # Add more 2D plot options as needed
        else:
            print(f"Unsupported 2D plot option: {plot_2d_nmr}")
    
    else:
        print("No valid plot keyword found in the input.")

# Function for 1D NMR plots
def create_1d_nmr_plot(data, title, x_axis_length):
    """
    Create and display a 1D NMR plot.

    Parameters:
    - data (list): List of NMR shifts.
    - title (str): Plot title.
    - x_axis_length (list): X-axis length.

    Returns:
    - None (displays the plot).
    """
    plt.figure(figsize=(x_axis_length, 4))
    plt.plot(data)
    plt.title(title)
    plt.xlabel("Chemical Shift (ppm)")
    plt.ylabel("Intensity")
    plt.grid(True)
    plt.show()

# Function for 2D NMR plots
def create_2d_nmr_plot(data, title, x_axis_length):
    """
    Create and display a 2D NMR plot.

    Parameters:
    - data (list): List of NMR shifts.
    - title (str): Plot title.
    - x_axis_length (list): X-axis length.

    Returns:
    - None (displays the plot).
    """
    plt.figure(figsize=(x_axis_length, x_axis_length))
    plt.imshow(data, cmap='viridis', aspect='auto', origin='lower')
    plt.title(title)
    plt.xlabel("X-axis / ppm")
    plt.ylabel("Y-axis / ppm")
    plt.colorbar(label="Intensity")
    plt.grid(False)
    plt.show()




#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# SECTION : IMPORT NMR SPECTRA WITH NMRGLUE

def find_nmr_files(working_dir):

    # Get the current working directory
    current_directory = os.getcwd()

    # Function to recursively search for files in a directory
    def search_nmr_files(directory):
        found_1r = False
        found_2rr = False

        for root, _, files in os.walk(directory):
            for filename in files:
                if "1r" in filename:
                    found_1r = True
                elif "2rr" in filename:
                    found_2rr = True

                # If both types of files are found, return immediately
                if found_1r and found_2rr:
                    return "both"

        if found_1r:
            return "1r"
        elif found_2rr:
            return "2rr"
        else:
            return None

    # Search for files in the current directory and its subfolders
    nmr_files = search_nmr_files(current_directory)

    # Handle based on the result
    if nmr_files == "both":
        
        # We have only 1D NMR experiments 
        with open(log_file, 'w') as f:
            f.write(f"Found 1r and 2rr files in folder {working_dir}. We will work with both 1D and 2D NMR data.\n")
        print("Found 1r and 2rr files. We will work with 1D and 2D NMR data.")

        # define data folders nmr_spectra\test_1d_1h\pdata\1\1r
        nmr_1d = 'D:/PhD/Data/CODES/PYTHON/(NC)2Ipy/DFT_output_on_ssNMR_spectra/nmr_spectra/test_1d_1h/pdata/1'

        nmr_2d = 'D:/PhD/Data/CODES/PYTHON/(NC)2Ipy/DFT_output_on_ssNMR_spectra/nmr_spectra/test_2d_1h_1h/pdata/1'
        
        # read data from Bruker
        dic, data_1d = ng.bruker.read_pdata(nmr_1d) 

        u = ng.bruker.guess_udic(dic, data_1d) # go look in the u variable for acq info

        plt.plot(data_1d)
        plt.show() 




    elif nmr_files == "1r": 
        # We have only 1D NMR experiments
        with open(log_file, 'w') as f:
            f.write(f"Found 1r files in folder {current_dir}. We will work with 1D NMR data.\n")
        print("Found 1r files. We will work with 1D NMR data.")
        

    elif nmr_files == "2rr":
        # We have only 2D NMR experiments
        with open(log_file, 'w') as f:
            f.write(f"Found 2rr files in folder {current_dir}. We will work with 2D NMR data.\n")
        print("Found 2rr files. We will work with 2D NMR data.")

    else:
        # No matching files found
        with open(error_file, 'w') as f:
            f.write(f"No 1r or 2rr files found in folder {current_dir}. Please check your folder. \n")
        print("No 1r or 2rr files found. Please check your folder")












#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# SECTION 1:












#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# SECTION : INITIALIZING ERROR AND LOG FILES

# Errors and logs
error_file = 'error_file.txt' # here you can find all the errors from the run
log_file = 'log_file.txt' # here saving everything necessary for publications

# Header log file
with open(log_file, 'w') as f:
    f.write(f"LOG FILE FOR PLOTTING DFT PARAMETERS ON NMR SPECTRA \n")
    f.write(f"ETTORE BARTALUCCI, AACHEN (DE) \n")
print("Generated log file")

# Paths
current_dir = os.getcwd()
with open(log_file, 'w') as f:
    f.write(f"We work in folder: {current_dir}")
print(f"Current working directory is: {current_dir}")

#------------------------------------------------------------------------------#



#------------------------------------------------------------------------------#
# SECTION MAIN

def main():
    """
    Main function
    """

    # Input file
    input_file = os.path.join(current_dir, 'DFT_output_on_ssNMR_spectra/input.txt')

    # Shielding file


    # Read input file
    read_input(input_file)

    # Read shielding file

    # find nmr files
    find_nmr_files(current_dir)
    


if __name__ == '__main__':
    main()
#------------------------------------------------------------------------------#






















