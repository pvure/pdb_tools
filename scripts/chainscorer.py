#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import sys
import time
from pyrosetta import *

def score_complex_and_chains(pdb_path, chain_a_id, chain_b_id, relax):
    """
    Calculates the Rosetta score for a two-chain complex and for each
    chain individually.

    Args:
        pdb_path (str): The path to the PDB file.
        chain_a_id (str): The chain identifier for the first partner (e.g., 'A').
        chain_b_id (str): The chain identifier for the second partner (e.g., 'B').
        relax (bool): If True, run FastRelax on the complex before scoring.
    """
    if not os.path.isfile(pdb_path):
        print(f"Error: The file '{pdb_path}' was not found.", file=sys.stderr)
        sys.exit(1)

    print(f"Loading PDB: {pdb_path}")
    
    # Initialize PyRosetta in silent mode
    init_options = "-ignore_unrecognized_res -ex1 -ex2 -out:level 100"
    init(extra_options=init_options)
    
    # Load the full complex pose
    pose = pose_from_pdb(pdb_path)
    
    # Create the standard ref2015 score function
    scorefxn = create_score_function('ref2015')

    if relax:
        print("\n--- Starting Structural Relaxation on the Complex ---")
        relax_start_time = time.time()
        fast_relax = rosetta.protocols.relax.FastRelax(scorefxn, 5)
        
        print(f"Score before relaxation: {scorefxn(pose):.2f} REU")
        fast_relax.apply(pose)
        print(f"Score after relaxation:  {scorefxn(pose):.2f} REU")
        
        relax_end_time = time.time()
        print(f"Relaxation finished in {relax_end_time - relax_start_time:.2f} seconds.")

    # 1. Score the entire complex
    scorefxn(pose)
    complex_score = pose.energies().total_energy()
    
    # 2. Split the pose into individual chains and score them
    chain_poses = pose.split_by_chain()
    pose_a = None
    pose_b = None

    # Find the correct pose for each chain ID
    for chain_pose in chain_poses:
        # Get the chain ID from the first residue of the split pose
        current_chain_id = chain_pose.pdb_info().chain(1)
        if current_chain_id == chain_a_id:
            pose_a = chain_pose
        elif current_chain_id == chain_b_id:
            pose_b = chain_pose
    
    if pose_a is None:
        print(f"Error: Could not find chain '{chain_a_id}' in the PDB file.", file=sys.stderr)
        sys.exit(1)
    if pose_b is None:
        print(f"Error: Could not find chain '{chain_b_id}' in the PDB file.", file=sys.stderr)
        sys.exit(1)

    # Score chain A
    scorefxn(pose_a)
    chain_a_score = pose_a.energies().total_energy()

    # Score chain B
    scorefxn(pose_b)
    chain_b_score = pose_b.energies().total_energy()

    # 3. Print all results
    print("\n" + "="*50)
    print("      Rosetta Energy Score Results")
    print("="*50)
    print(f"Total Complex (Chains {chain_a_id}+{chain_b_id}): {complex_score:>15.2f} REU")
    print(f"Chain {chain_a_id} alone:                     {chain_a_score:>15.2f} REU")
    print(f"Chain {chain_b_id} alone:                     {chain_b_score:>15.2f} REU")
    print("="*50)

# This block allows the script to be run from the command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate Rosetta scores for a complex and its individual chains."
    )
    parser.add_argument("pdb_file", type=str, help="Path to the input PDB file.")
    parser.add_argument("--chain_a", required=True, type=str, help="Identifier for the first chain (e.g., 'A').")
    parser.add_argument("--chain_b", required=True, type=str, help="Identifier for the second chain (e.g., 'B').")
    parser.add_argument("--relax", action="store_true", help="Perform FastRelax on the complex before scoring.")

    args = parser.parse_args()
    
    start_time = time.time()
    score_complex_and_chains(args.pdb_file, args.chain_a, args.chain_b, args.relax)
    end_time = time.time()
    
    print(f"\nTotal execution time: {end_time - start_time:.2f} seconds.")
