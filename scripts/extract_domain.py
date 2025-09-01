#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import sys

def extract_domain(pdb_file_path, start_res, end_res, output_file=None):
    """
    Extracts a specific range of residues from a PDB file and saves it
    to a new PDB file.

    Args:
        pdb_file_path (str): The full path to the input PDB file.
        start_res (int): The starting residue number of the domain.
        end_res (int): The ending residue number of the domain.
        output_file (str, optional): The path for the output PDB file.
                                     If None, a name is generated automatically.
    """
    if not os.path.isfile(pdb_file_path):
        print(f"Error: The file '{pdb_file_path}' was not found.", file=sys.stderr)
        sys.exit(1)

    if start_res > end_res:
        print(f"Error: Start residue ({start_res}) cannot be greater than end residue ({end_res}).", file=sys.stderr)
        sys.exit(1)

    directory = os.path.dirname(pdb_file_path) or '.'
    base_name = os.path.splitext(os.path.basename(pdb_file_path))[0]

    # Generate an output filename if one isn't provided
    if not output_file:
        output_file = os.path.join(directory, f"{base_name}_domain_{start_res}_{end_res}.pdb")

    print(f"Reading from: {pdb_file_path}")
    print(f"Extracting residues from {start_res} to {end_res}.")

    atom_count = 0
    try:
        with open(pdb_file_path, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        # Residue sequence number is in columns 23-26
                        res_num = int(line[22:26].strip())
                        if start_res <= res_num <= end_res:
                            outfile.write(line)
                            atom_count += 1
                    except (ValueError, IndexError):
                        # Ignore lines with malformed residue numbers
                        continue
        
        if atom_count > 0:
            print(f"Successfully wrote {atom_count} ATOM/HETATM records to: {output_file}")
        else:
            print(f"Warning: No residues found in the range {start_res}-{end_res}. An empty file was created.")

    except IOError as e:
        print(f"An error occurred during file processing: {e}", file=sys.stderr)
        sys.exit(1)

# This block allows the script to be run from the command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract a specific domain (residue range) from a single-chain PDB file."
    )
    parser.add_argument(
        "pdb_file",
        type=str,
        help="Path to the input PDB file."
    )
    parser.add_argument(
        "start_residue",
        type=int,
        help="The starting residue number of the domain to extract."
    )
    parser.add_argument(
        "end_residue",
        type=int,
        help="The ending residue number of the domain to extract."
    )
    parser.add_argument(
        "--output_file",
        type=str,
        help="Optional: Specify a name for the output PDB file."
    )

    args = parser.parse_args()
    extract_domain(args.pdb_file, args.start_residue, args.end_residue, args.output_file)
