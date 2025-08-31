#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import sys
from collections import defaultdict

def separate_pdb_chains(pdb_file_path):
    """
    Reads a PDB file and separates it into multiple new PDB files, one for
    each chain found in the original file.

    Args:
        pdb_file_path (str): The full path to the input PDB file.
    """
    if not os.path.isfile(pdb_file_path):
        print(f"Error: The file '{pdb_file_path}' was not found.", file=sys.stderr)
        sys.exit(1)

    directory = os.path.dirname(pdb_file_path)
    if not directory:
        directory = '.' # Handle case where file is in the current directory
    base_name = os.path.splitext(os.path.basename(pdb_file_path))[0]

    # Use a dictionary to hold the lines for each chain
    chains = defaultdict(list)
    found_chains = set()

    print(f"Reading from: {pdb_file_path}")
    
    try:
        with open(pdb_file_path, 'r') as infile:
            for line in infile:
                # Process ATOM, HETATM, and TER records which contain chain IDs
                if line.startswith(('ATOM', 'HETATM', 'TER')):
                    # The chain identifier is at a fixed position (index 21)
                    chain_id = line[21]
                    if chain_id.strip(): # Ensure there is a chain ID
                        chains[chain_id].append(line)
                        found_chains.add(chain_id)
    except IOError as e:
        print(f"An error occurred while reading the file: {e}", file=sys.stderr)
        sys.exit(1)

    if not found_chains:
        print("Warning: No chains with valid ATOM/HETATM records found in the file.")
        return

    print(f"Found chains: {', '.join(sorted(list(found_chains)))}")

    # Write the collected lines to their respective output files
    for chain_id, lines in chains.items():
        output_path = os.path.join(directory, f"{base_name}_chain_{chain_id}.pdb")
        try:
            with open(output_path, 'w') as outfile:
                for line in lines:
                    outfile.write(line)
            print(f"Successfully created Chain {chain_id} file: {output_path}")
        except IOError as e:
            print(f"An error occurred while writing file for chain {chain_id}: {e}", file=sys.stderr)

# This block allows the script to be run from the command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Separate a PDB file into distinct files for each chain. "
                    "The script automatically detects all chains present."
    )
    parser.add_argument(
        "pdb_file",
        type=str,
        help="Path to the input PDB file."
    )

    args = parser.parse_args()
    separate_pdb_chains(args.pdb_file)

