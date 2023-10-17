###################################################
# PLOT ORCA OUTPUT FOR ITERATIVE NMR CALCULATIONS #
#                                                 #
#               Ettore Bartalucci                 #
#               First: 16.09.2023                 #
#               Last:  03.10.2023                 #
#                                                 #
###################################################

# Import modules 
import os
import sys 
import numpy as np
import matplotlib.pyplot as plt
import re # RegEx

# -------------------------- DONE --------------------------------#
# ----------------------------------------------------------------#
# SECTION 1: READ ORCA OUTPUT FILES AND COUNT NUMBER OF JOBS
def count_jobs_number(output_file):
    """
    Read the output file from an ORCA calculation and count how many jobs have been run.
    """
    with open(output_file, 'r') as f:
        # Count number of jobs in output file
        lines = f.readlines()
        count = 0
        for line in lines:
            if line.strip().startswith("$$$$$$$$$$$$$$$$  JOB NUMBER"):
                count += 1
        return count
# ----------------------------------------------------------------#


# -------------------------- DONE --------------------------------#
# ----------------------------------------------------------------#
# SECTION 2: READ ORCA OUTPUT FILES AND EXTRACT LEVEL OF THEORY
def extract_level_of_theory(output_file):
    try:
        with open(output_file, 'r') as f:
            lines = f.readlines()
            
            # Search for the line containing "# Level of theory"
            for i, line in enumerate(lines):
                if "# Level of theory" in line:
                    # Extract the line immediately after it, this wont work if ppl dont use my syntax
                    level_of_theory_line = lines[i + 1].strip()

                    # Remove the line number from the line - do i want to keep this?
                    level_of_theory = level_of_theory_line.replace("| 10>", "").strip()

                    return level_of_theory

            return "Level of theory not found in the file."
    
    except FileNotFoundError:
        return f"File '{output_file}' not found."
    except Exception as e:
        return f"An error occurred: {str(e)}"

    # here write info to file?            
# ----------------------------------------------------------------#


# -------------------------- DONE --------------------------------#
# ----------------------------------------------------------------#
# SECTION 3: READ ORCA OUTPUT FILES AND SPLIT FOR NUMBER OF JOBS
def split_orca_output(output_file):
    """
    This function splits the huge ORCA multi-job file into individual job files.
    Then sew it back with the initial output lines from ORCA so that the files can
    be opened in Avogadro, otherwise it doesnt work.
    """
    if not os.path.isfile(output_file):
        print(f"Error in SECTION 3: ORCA output file '{output_file}' not found, please define.")
        return

    # Make use of RegEx for matching JOB lines
    job_matching = re.compile(r'\$+\s*JOB\s+NUMBER\s+(\d+) \$+')

    # Extract initial ORCA text before job files to append to splitted files
    initial_content = extract_initial_output_content(output_file, job_matching)

    with open(output_file, 'r') as f:
        current_job_n = None # current job number
        current_job_content = [] # job specific info

        for line in f:
            match = job_matching.search(line) # regex match search
            if match:
                # if match is found, write to file
                if current_job_n is not None:
                    output_file_path = f'split_output/splitted_orca_job{current_job_n}.out'
                    with open(output_file_path, 'w') as out:
                        out.writelines(initial_content + current_job_content) # initial orca info + job specific info
                    
                    print(f'Wrote job {current_job_n} to {output_file_path}')

                current_job_n = match.group(1)
                current_job_content = [] # reset

            current_job_content.append(line)

        # write last job to file
        if current_job_n is not None:
            output_file_path = f'split_output/splitted_orca_job{current_job_n}.out'
            with open(output_file_path, 'w') as out:
                out.writelines(initial_content + current_job_content)
            
            print(f'Wrote job {current_job_n} to {output_file_path}')

    print(f'ORCA output has been split into {current_job_n} sub files for further analysis')

# This adds the initial content necessary for avogadro visualization to each splitted file
def extract_initial_output_content(output_file, job_matching):
    """
    Add the initial ORCA file until JOB NUMBER to each file to be read by Avogadro
    """
    initial_content = []
    with open(output_file, 'r') as f:
        for line in f:
            match = job_matching.search(line) # if u match stop there
            if match:
                break # break as soon as you see the first job line 
            initial_content.append(line) # and append to file
    return initial_content
# ----------------------------------------------------------------#


# -------------------------- DONE --------------------------------#
# ----------------------------------------------------------------#
# SECTION 4: READ ORCA PROPERTY FILES AND EXTRACT SHIELDINGS
def read_property_file(property_file, job_number):
    """
    Read the property file from an ORCA calculation. Extract CSA shieldings (shielding tensor (ppm))
    and shifts (P(iso)) for each nucleus.

    Input:
    property_file:
        is the file that comes as outcome from ORCA calculations containing the most important informations
        on the simulations. We want the NMR parameter, but in principle it contains a summary of the ORCA
        .mpi8.out file.
    job_number:
        number of ORCA jobs ran
    Output:
    shielding_{job_number}.txt:
        file with condensed NMR data ready for plotting.
        It contains csa_tensor_data for the shielding tensor, shifts_data for the isotropic shifts
    """

    # Dictionary to store NMR data for each nucleus
    csa_tensor_data = {}
    shifts_data = {}

    nucleus_info = None
    reading_shielding = False

    with open(property_file, 'r') as f:
        for line in f:
            if "Nucleus:" in line:
                nucleus_info = line.strip().split()
                nucleus_index = int(nucleus_info[1])
                nucleus_name = nucleus_info[2]
                csa_tensor_data[(nucleus_index, nucleus_name)] = []
                reading_shielding = True
            elif reading_shielding and "Shielding tensor (ppm):" in line:
                # Skip the header line
                next(f)
                shielding_values = []
                for _ in range(3):
                    tensor_line = next(f).split()
                    shielding_values.append([float(val) for val in tensor_line])
                csa_tensor_data[(nucleus_index, nucleus_name)] = shielding_values
            elif "P(iso)" in line:
                shifts = float(line.split()[-1])
                shifts_data[(nucleus_index, nucleus_name)] = shifts
                reading_shielding = False

    # Write the extracted data to shieldings.txt
    output_shielding_path = f'nmr_data/shieldings_{job_number}.txt'
    for job_n in range(1, job_number +1):
        with open(output_shielding_path, "w") as output_f:
            for (nucleus_index, nucleus_name), shielding_values in csa_tensor_data.items():
                output_f.write(f"Nucleus {nucleus_index} {nucleus_name}\nShielding tensor (ppm):\n")
                for values in shielding_values:
                    values = values[1:]
                    output_f.write("\t".join(map(str, values)) + "\n")
                shifts = shifts_data[(nucleus_index, nucleus_name)]
                output_f.write(f"Isotropic Shift: {shifts}\n")
                output_f.write("\n")
        
        print(f"Shieldings extracted and saved to 'nmr_data/shieldings_job{job_number}.txt'.")


# ----------------------------------------------------------------#


# -------------------------- Start --------------------------------#
# ----------------------------------------------------------------#
# SECTION 5: READ ORCA OUTPUT FILE FOR EXTRACTING SCALAR COUPLINGS
def read_couplings(output_file):
    """
    Read the output file from an ORCA calculation. Extract scalar couplings for each nucleus

    Input:
    output_file:
        is the file that comes as outcome from ORCA calculations containing the most important informations
        on the simulations. We want the NMR parameter, which are towards the end of the ORCA .mpi8.out file.
    Output:
    j_couplings.txt:
        file with condensed NMR data ready for plotting.
        It contains scalar J couplings in pseudo-table format
    """

    # Dictionary to store NMR data for each nucleus
    j_coupling_data = {}

    # List to store the order of nuclei
    nuclei_order = []

    reading_couplings = False

    # Define regular expressions for start and end markers
    start_marker = re.compile(r'^\s*SUMMARY\s+OF\s+ISOTROPIC\s+COUPLING\s+CONSTANTS\s+\(Hz\)')
    end_marker = re.compile(r'Maximum memory used throughout the entire EPRNMR-calculation:')

    with open(output_file, 'r') as f:
        for line in f:
            line = line.strip()

            # Check for start marker
            if start_marker.search(line):
                reading_couplings = True
                print("Entering the reading couplings block.")
                continue  # Start reading couplings
            
            # Check for end marker
            if end_marker.search(line):
                reading_couplings = False
                print("Exiting the reading couplings block.")
                break  # Stop reading couplings when this line is encountered

            # If we're reading couplings, process the line
            if reading_couplings and line:
                data_values = re.split(r'\s+', line)
                nucleus = data_values[0]
                j_couplings = data_values[1:]
                nuclei_order.append(nucleus)  # Add the nucleus to the order list
                j_coupling_data[nucleus] = j_couplings
                
    # Write the formatted J coupling data to j_couplings.txt
    with open('nmr_data/j_couplings.txt', 'w') as output_file:
        # Write the header row with nuclei information
        output_file.write('\t'.join(nuclei_order) + '\n')
        
        # Write the data rows
        for nucleus, j_couplings in j_coupling_data.items():
            output_file.write(f"{nucleus}\t{' '.join(j_couplings)}\n")

    print("J couplings extracted and saved to 'nmr_data/j_couplings.txt'.")

# ----------------------------------------------------------------#


# ----------------------------------------------------------------#
# SECTION X: PRINT SIMULATION INFO FROM FILES


# ----------------------------------------------------------------#


# ----------------------------------------------------------------#
# SECTION X: PLOTTING NMR DATA
def plot_csa(shielding_tensor):
    """
    Load the previously extracted shielding values and plot them as a function of distance
    """


def plot_isotropic_shifts(shielding_tensor):

    nucleus_data = {}  # Dictionary to store data for each nucleus type
    
    with open(shielding_tensor, 'r') as f:
        lines = f.readlines()

    current_nucleus = None
    current_data = []

    for line in lines:
        line = line.strip()
        if line.startswith('Nucleus'):
            if current_nucleus is not None:
                nucleus_data[current_nucleus] = current_data
            current_nucleus = line.split()[-1]
            current_data = []
        elif line.startswith('Isotropic Shift:'):
            isotropic_shift = float(line.split()[-1])
            current_data.append(isotropic_shift)

    if current_nucleus is not None:
        nucleus_data[current_nucleus] = current_data

    # Plot one plot per nucleus type
    for nucleus, shifts in nucleus_data.items():
        plt.figure(figsize=(8, 6))
        plt.plot(shifts, marker='o', linestyle='-', label=nucleus)
        plt.title(f'Isotropic Shifts for {nucleus}')
        plt.xlabel('Data Point')
        plt.ylabel('Isotropic Shift (ppm)')
        plt.legend()
        plt.grid(True)
        plt.show()
# ----------------------------------------------------------------#

# ----------------------------------------------------------------#
# LOGS AND ERRORS
error_log_file = 'error_log_file.txt' 
log_file = 'log_file.txt'  # this file will be the summary for the publications and for database for machine learning

# RELATIVE PATH
current_dir = os.getcwd()
print(f'Current working directory is: {current_dir}')

# SECTION MAIN: Call all modules
def main():
    """
    Main function.
    """

    # --- Wroking with the orca MPI8 output file --- #
    # MPI8 output file
    orca_output = os.path.join(current_dir, 'myoutput.out')

    # Count number of simulation jobs that were ran
    n_jobs = count_jobs_number(orca_output)
    print(f'Number of ORCA DFT calculations in file: {n_jobs}')

    # Extract level of theory
    lot_out = extract_level_of_theory(orca_output)
    print(lot_out)

    # Read J couplings
    # read_couplings(orca_output)

    # Split orca output in several subfiles
    #split_orca_output(orca_output)
    # ---------------------------------------------- #

    # --- Working with input file (future feature) --- #
    # Input file
    input_file = os.path.join(current_dir, 'input_file.txt')
    # ------------------------------------------------ #

    # --- Working with property files from ORCA --- #
    # Property files = number of jobs
    for job_number in range (1, n_jobs + 1):

        # Property files from all jobs (except first)
        orca_properties = os.path.join(current_dir, f'properties/run_all_displaced_distances_job{job_number}_property.txt')

        # Read each of them
        read_property_file(orca_properties, job_number)

        # Plot NMR data 
        shielding_data = os.path.join(current_dir, f'nmr_data/shieldings_{job_number}.txt')
        plot_isotropic_shifts(shielding_data)


if __name__ == '__main__':
    main()