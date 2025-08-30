import argparse
import os

def separate_pdb_chains(pdb_file_path):
    """
    Reads a PDB file and separates ATOM, HETATM, and TER records for
    chains A and B into two new, separate PDB files.

    Args:
        pdb_file_path (str): The full path to the input PDB file.
    """
    if not os.path.isfile(pdb_file_path):
        print(f"Error: The file '{pdb_file_path}' was not found.")
        return

    directory = os.path.dirname(pdb_file_path)
    base_name = os.path.splitext(os.path.basename(pdb_file_path))[0]

    # Define the output file names based on the input file name
    output_path_a = os.path.join(directory, f"{base_name}_chain_A.pdb")
    output_path_b = os.path.join(directory, f"{base_name}_chain_B.pdb")

    print(f"Reading from: {pdb_file_path}")

    try:
        with open(pdb_file_path, 'r') as infile, \
             open(output_path_a, 'w') as outfile_a, \
             open(output_path_b, 'w') as outfile_b:

            last_written_chain = None

            for line in infile:
                # look at lines that describe atomic coordinates
                if line.startswith(('ATOM', 'HETATM')):
                    # chain identifier is in column 22
                    chain_id = line[21]

                    if chain_id == 'A':
                        outfile_a.write(line)
                        last_written_chain = 'A'
                    elif chain_id == 'B':
                        outfile_b.write(line)
                        last_written_chain = 'B'

                # look for 'ter' record, shows termination of a chain
                elif line.startswith('TER'):
                    if last_written_chain == 'A':
                        outfile_a.write(line)
                    elif last_written_chain == 'B':
                        outfile_b.write(line)
                    last_written_chain = None

        print(f"Successfully created Chain A file: {output_path_a}")
        print(f"Successfully created Chain B file: {output_path_b}")

    except IOError as e:
        print(f"An error occurred during file processing: {e}")


# run from command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Separate a PDB file with chains A and B into two distinct files. "
                    "The output files will be named based on the input file."
    )
    # arg for input file path
    parser.add_argument(
        "pdb_file",
        type=str,
        help="Path to the input PDB file."
    )

    # parse user input args
    args = parser.parse_args()

    # call main function
    separate_pdb_chains(args.pdb_file)
