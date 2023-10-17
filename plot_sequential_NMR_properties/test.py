def main():
    """
    Main function.
    """

    # MPI8 output file
    orca_output = os.path.join(current_dir, 'myoutput.out')

    # Input file
    input_file = os.path.join(current_dir, 'input_file.txt')

    # Count number of simulation jobs that were ran
    n_jobs = count_jobs_number(orca_output)
    print(f'Number of ORCA DFT calculations in file: {n_jobs}')

    for job_number in range(1, n_jobs + 1):
        # Property file for the current job
        orca_properties = os.path.join(current_dir, f'myoutput_property_{job_number}.txt')

        # Read the property file for each job to extract shieldings
        read_property_file(orca_properties)

        # Read the output file to extract couplings
        read_couplings(orca_output)

        # Split orca output in several subfiles (if needed)
        # split_orca_output(orca_output)

        # Extract level of theory (if needed)
        lot_out = extract_level_of_theory(orca_output)
        print(lot_out)

        # Plot NMR data (if needed)
        shielding_data = os.path.join(current_dir, 'nmr_data/shieldings.txt')
        plot_isotropic_shifts(shielding_data)

if __name__ == '__main__':
    main()
