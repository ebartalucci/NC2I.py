##########################################################################
#      DISTANCE SCANNER FOR MOLECULAR FRAGMENTS IN DFT CALCULATIONS      #
#                          ------------------                            #
#                          v.1.3 / 10.09.23                              #
#                          ETTORE BARTALUCCI                             #
##########################################################################

# Importing necessary modules
import os
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Required for 3D plot
from sklearn.cluster import KMeans
from timeit import default_timer as timer

#-------------------------------------------------------------------#
#------------------------------FUNCTIONS----------------------------#
#-------------------------------------------------------------------#

# SECTION 1: READING XYZ ATOMIC COORDINATES FROM FILE AND SPLIT THE FRAGMENTS
def read_atomic_coord(file_path):
    """
    This function reads the geometry optimized atomic coordinates of the two fragments
    you want to displace, the input is the classical .xyz coordinate file that can
    be then feeded to ORCA after displacement.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()[2:]  # skip the first two lines (atom count and comment)
        coordinates = []
        atom_identities = [] # this is the one storing atom identity info, need to append it at a later step
        for line in lines:
            atom_data = line.split()
            atom_identity = atom_data[0] 
            coord = [float(atom_data[1]), float(atom_data[2]), float(atom_data[3])]
            coordinates.append(coord)
            atom_identities.append(atom_identity)
        coordinates = np.array(coordinates)
        atom_identities = np.array(atom_identities)
        return coordinates, atom_identities

#-------------------------------------------------------------------#

# SECTION 2: CALCULATION OF CENTROIDS AND STORING VALUES IN FILE
def calculate_centroids(coordinates):
    """
    This function computes centroids for the defined fragments.
    """
    num_atoms = len(coordinates)
    centroids = []
    for fragment in coordinates:
        centroid = np.sum(fragment, axis=0) / len(fragment)
        centroids.append(centroid)
    return np.array(centroids)

def write_centroids(file_path, centroids):
    """
    This function writes centroids to file for furhter manipulation.
    """
    with open(file_path, 'w') as f:
        f.write(f'Centroid coordinates:\n')
        for centroid in centroids:
            f.write(f'{centroid[0]:.6f} {centroid[1]:.6f} {centroid[2]:.6f}\n')

#-------------------------------------------------------------------#

# SECTION 3: CHECKPOINT COMPUTE TOPOLOGY AND K-NEAREST CLUSTERING FOR MOLECULAR CENTROIDS
def plot_starting_molecular_fragments(coords1, coords2, centroids):
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    ax1.scatter(coords1[:, 0], coords1[:, 1], coords1[:, 2], color='blue', label='Molecule 1')
    ax1.scatter(coords2[:, 0], coords2[:, 1], coords2[:, 2], color='red', label='Molecule 2')
    ax1.scatter(centroids[:, 0], centroids[:, 1], centroids[:, 2], color='green', marker='x', s=100, label='Centroids')
    ax1.legend()

#-------------------------------------------------------------------#

# SECTION 4: READ USER PROVIDED INPUT FILE, SPLIT AND ASSIGN MOLECULAR FRAGMENTS TO INDIVIDUAL MOLECULES
def assign_molecule_fragments(coordinates, input_file):
    """
    This function does
    """
    with open(input_file, 'r') as f:
        lines = f.readlines()
        fragment1_indices = []
        fragment2_indices = []
        current_fragment = None
        for line in lines:
            if line.strip() == "$fragment1":
                current_fragment = fragment1_indices
            elif line.strip() == "$fragment2":
                current_fragment = fragment2_indices
            elif line.strip() == "$displacement":
                break
            else:
                index = int(line.strip()) - 1  # Convert 1-based index to 0-based index
                current_fragment.append(index)    
    coords1 = coordinates[fragment1_indices]
    coords2 = coordinates[fragment2_indices]
    return coords1, coords2

#-------------------------------------------------------------------#

# SECTION 5: DISPLACE THE TWO FRAGMENT ALONG THE CENTROID LINE MOVING ONE AND FIXING THE OTHER
def displace_fragment(coords1, displacement_direction, displacement_step, i):
    """
    This function displaces the fragment along the displacement_direction vector.
    """
    displacement_direction /= np.linalg.norm(displacement_direction)  # Normalize the displacement direction vector
    displacement_vector = displacement_direction * displacement_step * i # Displace along the normalized direction
    print(f'Displacement vector: {displacement_vector}')
    print(type(displacement_vector))
    return coords1 - displacement_vector  # Apply displacement by adding the vector
    
    # From olivia
    #displacement_vector = coords1 - displacement_direction * i  # displace along one axis 
    #print(f'Displacement vector:{displacement_vector}') # write this to log file
    #print(type(displacement_vector)) # write this to log file
    #return displacement_vector

#-------------------------------------------------------------------#

# SECTION 6: WRITE NEW DISPLACED COORDINATES TO FILES
def write_displaced_xyz_file(file_path, coords_fixed, coords_displaced, atom_identities):
    """
    This function writes both the fixed and displaced coordinates to an XYZ file.
    """
    with open(file_path, 'w') as f:
        num_atoms = len(coords_fixed) + len(coords_displaced)
        f.write(f'{num_atoms}\n')
        f.write(f'Step {file_path.stem}\n')

        # Write fixed fragment coordinates
        for i, atom in enumerate(coords_fixed):
            f.write(f'{atom_identities[i+len(coords_displaced)]} {atom[0]:.6f} {atom[1]:.6f} {atom[2]:.6f}\n') # i+len(coord_fixed) to skip to the fixed fragment indices
        
        # Write displaced fragment coordinates
        for i, atom in enumerate(coords_displaced):
            f.write(f'{atom_identities[i]} {atom[0]:.6f} {atom[1]:.6f} {atom[2]:.6f}\n') 
        

#-------------------------------------------------------------------#

# SECTION 7: count fragments for K-means clustering
def count_fragments(input_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()
        count = 0
        for line in lines:
            if line.strip().startswith("$fragment"):
                count += 1
        return count

#-------------------------------------------------------------------#

# SECTION 8: COMPUTE DISTANCE OF DISPLACED ATOMS FROM FIXED CENTROID
def compute_distance_from_centroid(coord_displaced, centroids):

    # Convert the lists to NumPy arrays for element-wise operations
    coord_displaced = np.array(coord_displaced)

    # Assuming fragment_centroids contains two sets of centroids
    if len(centroids) == 2:
        fixed_centroid = centroids[1]
        fixed_centroid = np.array(fixed_centroid)

        # Compute distance 
        distance_to_centroid = np.linalg.norm(coord_displaced - fixed_centroid, axis=1)

        return distance_to_centroid
    else:
        print('Error: fragment_centroids should contain two sets of centroids, if more please modify Section 8')

#-------------------------------------------------------------------#

# SECTION 9: WRITE DISTANCES TO FILES
def write_distances_file(file_path, coords_displaced, distance_to_centroid, atom_identities, displacement_step):
    """
    This function writes the distances between fixed centroid and displaced coordinates to a file.
    """
    with open(file_path, 'w') as f:
        num_atoms = len(coords_displaced)
        f.write(f'Number of atoms: {num_atoms}\n')
        f.write(f'Step {file_path.stem}\n')
        f.write(f'Displacement step: {displacement_step} A\n')
         
        # Write displaced fragment coordinates
        for i, atom in enumerate(coords_displaced):
            f.write(f'{atom_identities[i]} {distance_to_centroid[i]}\n') 


#-------------------------------------------------------------------#

#-------------------------------------------------------------------#
#----------------------------END FUNCTIONS--------------------------#
#-------------------------------------------------------------------#

# WORKING WITH RELATIVE PATHS 
current_dir = os.getcwd()
print(current_dir)

# LOGS AND ERRORS
error_log_file = 'error_log_file.txt' # to finish
log_file = 'log_file.txt' # to finish


# START TIMER: COMPUTE EFFECTIVE WALL TIME
start = timer() # this is not in [sec] i think

# SECTON: MAIN
def main():

    # Relative paths
    start_structure = os.path.join(current_dir, 'Distances/input_structures/final_pi-pi-T_shape_full_opt_r2scan3c_def2mtzvpp.xyz')
    centroid_out = os.path.join(current_dir, 'Distances/centroid_output/centroid_file.xyz')
    input_file = os.path.join(current_dir, 'Distances/input_file/input_file.txt')

    # Read xyz file: this should be either a fully optimized geometry or one with relaxed H
    coordinates, atom_identities = read_atomic_coord(start_structure)
    print(f'Starting coordinates: {coordinates}')
    print(f'Atom identities: {atom_identities}')

    # Assign coordinates to molecular fragments, check nomenclature of your atoms in avogadro or any other molecular graphics soft
    coords1, coords2 = assign_molecule_fragments(coordinates, input_file)

    # Concatenate coordinates for k-means clustering
    all_coords = np.concatenate((coords1, coords2), axis=0)
    # print(f'All coords: {all_coords}')

    # Count how many fragments you have defined in the input file, important for accurate K-means clustering
    n_fragments = count_fragments(input_file)
    print(f"Number of '$fragment' occurrences: {n_fragments}")

    # Perform k-means clustering to compute centroids
    kmeans = KMeans(n_clusters=n_fragments) # K-means clusters = number of centroids = number of fragments
    kmeans.fit(all_coords)
    centroids = kmeans.cluster_centers_

    # Compute centroids for each fragment
    fragment_centroids = calculate_centroids([coords1, coords2])

    # Write centroid coordinates to file
    write_centroids(centroid_out, fragment_centroids)
    print(f'Centroid coordinates: {fragment_centroids}')

    # Calculate displacement direction (line connecting centroids)
    displacement_direction = centroids[1] - centroids[0]
    displacement_direction /= np.linalg.norm(displacement_direction)
    print(f'Displacement direction:{displacement_direction}')

    # Read displacement step size from input file
    displacement_step = None
    with open(input_file, 'r') as f:
        lines = f.readlines()
        read_displacement = False
        for line in lines:
            if read_displacement:
                displacement_values = line.strip().split()
                if displacement_values:
                    displacement_step = float(displacement_values[0])
                    break
            elif line.strip() == "$displacement":
                read_displacement = True

    if displacement_step is None:
        print('ERROR: displacement step size not found in input file, please specify it! Syntax => $displacement + number')
        return
    print(f'Displacement step is: {displacement_step}') # please doublecheck that it is the same value you defined in the input

    # Displace the first fragment iteratively and save each structure
    displaced_fragment_coords = coords1.copy()  # Make a copy of the original coordinates of the fragment that is displaced
    print(f'Original coordinates displaced fragment:', displaced_fragment_coords)

    # Initialize the coordinates for the fixed fragment (e.g., coords2)
    coords_fixed = coords2.copy() # make a copy of the fixed fragment coordinates to append to the displaced ones
    print(f'Original coordinates fixed fragment:', coords_fixed)

    all_displaced_fragment_coords = [displaced_fragment_coords]  # List to store all displaced structures

    # Combine displaced coordinates with original ones
    all_combined_coords = [np.concatenate((coords_fixed, displaced_fragment_coords), axis=0)]  # List to store all combined structures

    fragment_centroids = [fragment_centroids[0]]  # List to store all centroids

    # Dissociation limit
    diss_lim = 20 # change with the output value in agnstrom from func(dissociation_limit)

    for i in range(1, diss_lim):  # Iterate 10 times (adjust the number as needed) put this as to be the dissociation limit of each DFT run
        
        displacement_vector = [] 

        # Compute new set of coordinates for displaced fragments, change $displacement value in input file to tune the displacement
        displaced_fragment_coords = displace_fragment(coords1, displacement_direction, displacement_step, i)
        #print(f'Displaced fragment coord is: {displaced_fragment_coords}')

        combined_coords = np.concatenate((coords_fixed, displaced_fragment_coords), axis=0)
        all_combined_coords.append(combined_coords)

        # Update centroids for the displaced structure
        fragment_centroid = calculate_centroids([displaced_fragment_coords])
        fragment_centroids.append(fragment_centroid[0])
        print(f'Updated centroids:', fragment_centroid)

        # Write displaced structure to file
        output_file = Path(os.path.join(current_dir, f'Distances/displaced_structures/displaced_structure_{i}.xyz'))
        write_displaced_xyz_file(output_file, coords_fixed, displaced_fragment_coords, atom_identities)

        all_displaced_fragment_coords.append(displaced_fragment_coords)

        # Compute distance between the fixed fragment centroid and all the atoms from the displaced fragment
        centroid_to_displaced_distance = compute_distance_from_centroid(displaced_fragment_coords, centroids)
        print(f'Distance between displaced coordinates and centroid is: {centroid_to_displaced_distance}')

        # Write distances to file - needed for DFT calculations outputs
        distance_output_file = Path(os.path.join(current_dir, f'Distances/distance_files/distances_structures_{i}.xyz'))
        write_distances_file(distance_output_file, displaced_fragment_coords, centroid_to_displaced_distance, atom_identities, displacement_step)


    # Plot initial topology for molecular fragments and centroids
    fig = plot_starting_molecular_fragments(coords1, coords2, centroids)

    # Generate colors for the plots based on displacement iteration
    num_iterations = len(all_displaced_fragment_coords)
    colors = plt.cm.viridis(np.linspace(0.2, 1.0, num_iterations))

    # Plot displaced molecular fragments and centroids
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot original fragments and centroids
    ax.scatter(coords1[:, 0], coords1[:, 1], coords1[:, 2], color=colors[0], label='Molecule 1 (Original)')
    ax.scatter(coords2[:, 0], coords2[:, 1], coords2[:, 2], color=colors[0], label='Molecule 2 (Original)')
    ax.scatter(centroids[:, 0], centroids[:, 1], centroids[:, 2], color=colors[0], marker='x', s=100, label='Centroids (Original)')

    # Plot displaced fragments and centroids
    for i, displaced_coords in enumerate(all_displaced_fragment_coords[1:], start=1):
        color = colors[i]
        label = f'Iteration {i}'
        ax.scatter(displaced_coords[:, 0], displaced_coords[:, 1], displaced_coords[:, 2], color=color, label=label)
        ax.scatter(fragment_centroids[i][0], fragment_centroids[i][1], fragment_centroids[i][2], color=color, marker='x', s=100, label=f'Centroids ({label})')

    ax.legend()
    plt.show()


    # END TIMER: STOP TIMER AND PRINT
    elapsed_time = timer() - start  # in seconds
    print(f'Elapsed time for the code to run is: {elapsed_time}')


if __name__ == '__main__':
    main()
