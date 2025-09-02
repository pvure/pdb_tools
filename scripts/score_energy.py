#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import sys
import time
from pyrosetta import *

def score_pdb(pdb_path, relax):
    """
    Calculates the total Rosetta score for a PDB file.

    Args:
        pdb_path (str): The path to the PDB file.
        relax (bool): If True, run FastRelax before scoring.
    """
    if not os.path.isfile(pdb_path):
        print(f"Error: The file '{pdb_path}' was not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading PDB: {pdb_path}")
    
    # Initialize PyRosetta in silent mode and add common options
    init_options = "-ignore_unrecognized_res -ex1 -ex2 -out:level 100"
    init(extra_options=init_options)
    
    # Load the pose
    pose = pose_from_pdb(pdb_path)
    
    # Create the standard ref2015 score function
    scorefxn = create_score_function('ref2015')

    if relax:
        print("\n--- Starting Structural Relaxation (FastRelax) ---")
        relax_start_time = time.time()
        
        fast_relax = rosetta.protocols.relax.FastRelax(scorefxn, 5)
        
        print(f"Score before relaxation: {scorefxn(pose):.2f} REU")
        fast_relax.apply(pose)
        print(f"Score after relaxation:  {scorefxn(pose):.2f} REU")
        
        relax_end_time = time.time()
        print(f"Relaxation finished in {relax_end_time - relax_start_time:.2f} seconds.")
    else:
        # Just score the pose if not relaxing
        scorefxn(pose)

    total_score = pose.energies().total_energy()

    print("\n" + "="*50)
    print(f"Final Calculated Score: {total_score:.2f} REU")
    print("="*50)

# This block allows the script to be run from the command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate the total Rosetta score for a given PDB file."
    )
    parser.add_argument(
        "pdb_file",
        type=str,
        help="Path to the input PDB file to score."
    )
    parser.add_argument(
        "--relax",
        action="store_true",
        help="Perform a FastRelax on the structure before scoring."
    )

    args = parser.parse_args()
    
    start_time = time.time()
    score_pdb(args.pdb_file, args.relax)
    end_time = time.time()
    
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds.")
