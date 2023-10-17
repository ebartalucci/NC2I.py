import nmrglue as ng
import os
import matplotlib.pyplot as plt



#D:\PhD\Data\DFT\NONCOV\Scripts\DFT_output_on_ssNMR_spectra\nmr_spectra\test_1d_1h\fid

# Read Bruker nmr FID data, you need the path written with '/' instead of '\'
#dic,data = ng.bruker.read_pdata('D:/PhD/Data/DFT/NONCOV/Scripts/DFT_output_on_ssNMR_spectra/nmr_spectra/test_1d_1h/pdata/1')

# Set the spectral parameters
# Guess universal dictionaries for spectra info based on Bruker specifications
#u = ng.bruker.guess_udic(dic, data) # go look in the u variable for acq info


#plt.plot(data)
#plt.show()



def find_nmr_files():
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
        # Do something for both "1r" and "2rr" files
        print("Found both 1r and 2rr files. Handling both...")
        print(f'{current_directory}')
        # Replace this with your specific action for both cases

    elif nmr_files == "1r":
        # Do something for "1r" files
        print("Found 1r files. Handling 1r...")
        # Replace this with your specific action for "1r" files

    elif nmr_files == "2rr":
        # Do something for "2rr" files
        print("Found 2rr files. Handling 2rr...")
        # Replace this with your specific action for "2rr" files

    else:
        # No matching files found
        print("No 1r or 2rr files found.")

# Example usage:
find_nmr_files()