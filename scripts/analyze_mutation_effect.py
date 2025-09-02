#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import time
from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax

def get_user_input():
    """
    Prompts the user for all necessary inputs to run the analysis.
    """
    inputs = {}
    print("--- Please provide the required information ---")
    
    # Get PDB file path
    while True:
        path = input("Enter the full path to your PDB file: ").strip()
        if os.path.isfile(path):
            inputs['pdb_path'] = path
            break
        else:
            print(f"Error: File not found at '{path}'. Please try again.")

    # Get chain IDs
    inputs['target_chain'] = input("Enter the target chain ID (e.g., A): ").strip().upper()
    inputs['binder_chain'] = input("Enter the binder chain ID (e.g., B): ").strip().upper()

    # Get mutation details
    print("\n--- Mutation Details ---")
    inputs['mutation_chain'] = input(f"Which chain to mutate? ({inputs['target_chain']} or {inputs['binder_chain']}): ").strip().upper()
    inputs['residue_id'] = int(input("Enter the residue number to mutate: ").strip())
    inputs['original_res'] = input("Enter the original 3-letter amino acid code (e.g., SER): ").strip().upper()
    inputs['new_res'] = input("Enter the new 3-letter amino acid code (e.g., TRP): ").strip().upper()

    # Ask about relaxation
    while True:
        relax_choice = input("\nPerform structural relaxation? (yes/no): ").strip().lower()
        if relax_choice in ['yes', 'y']:
            inputs['relax'] = True
            break
        elif relax_choice in ['no', 'n']:
            inputs['relax'] = False
            break
        else:
            print("Invalid input. Please enter 'yes' or 'no'.")
            
    return inputs

def perform_mutation(pdb_path, chain, res_id, original, new):
    """
    Creates a new temporary PDB file with the specified mutation.
    This uses a simple text-replacement method.
    """
    mutated_pdb_path = "mutated_temp.pdb"
    print(f"\nCreating mutated PDB file: {mutated_pdb_path}")

    with open(pdb_path, 'r') as infile, open(mutated_pdb_path, 'w') as outfile:
        for line in infile:
            if line.startswith('ATOM'):
                line_chain = line[21]
                try:
                    line_res_id = int(line[22:26].strip())
                    line_res_name = line[17:20].strip()

                    if line_chain == chain and line_res_id == res_id and line_res_name == original:
                        modified_line = line[:17] + new.ljust(3) + line[20:]
                        outfile.write(modified_line)
                    else:
                        outfile.write(line)
                except ValueError:
                    outfile.write(line)
            else:
                outfile.write(line)
    return mutated_pdb_path

def analyze_pose(pose, scorefxn, target_chain, binder_chain, relax_protocol):
    """
    Analyzes a pose to get complex and individual chain scores.
    Optionally relaxes the pose first.
    """
    if relax_protocol:
        print("\n--- Starting Structural Relaxation ---")
        relax_start_time = time.time()
        print(f"Score before relaxation: {scorefxn(pose):.2f} REU")
        relax_protocol.apply(pose)
        print(f"Score after relaxation:  {scorefxn(pose):.2f} REU")
        relax_end_time = time.time()
        print(f"Relaxation finished in {relax_end_time - relax_start_time:.2f} seconds.")
    
    # Score the complex
    scorefxn(pose)
    complex_score = pose.energies().total_energy()

    # Split and score individual chains
    chain_poses = pose.split_by_chain()
    pose_target, pose_binder = None, None
    for p in chain_poses:
        if p.pdb_info().chain(1) == target_chain:
            pose_target = p
        elif p.pdb_info().chain(1) == binder_chain:
            pose_binder = p
    
    scorefxn(pose_target)
    target_score = pose_target.energies().total_energy()
    
    scorefxn(pose_binder)
    binder_score = pose_binder.energies().total_energy()

    return {
        'complex': complex_score,
        'target': target_score,
        'binder': binder_score
    }

def main():
    """
    Main function to run the interactive mutation and scoring workflow.
    """
    total_start_time = time.time()
    
    # Get all parameters from the user
    params = get_user_input()

    # Initialize PyRosetta
    init_options = "-ignore_unrecognized_res -ex1 -ex2 -out:level 100"
    init(extra_options=init_options)
    
    scorefxn = create_score_function('ref2015')
    relax_protocol = FastRelax(scorefxn, 5) if params['relax'] else None

    # --- 1. Analyze the Wild-Type (Original) Structure ---
    print("\n" + "="*50)
    print("    ANALYZING WILD-TYPE STRUCTURE")
    print("="*50)
    wt_pose = pose_from_pdb(params['pdb_path'])
    wt_scores = analyze_pose(wt_pose, scorefxn, params['target_chain'], params['binder_chain'], relax_protocol)

    # --- 2. Create and Analyze the Mutated Structure ---
    mutated_pdb_file = perform_mutation(
        params['pdb_path'], params['mutation_chain'], params['residue_id'], 
        params['original_res'], params['new_res']
    )
    
    print("\n" + "="*50)
    print("      ANALYZING MUTATED STRUCTURE")
    print("="*50)
    mut_pose = pose_from_pdb(mutated_pdb_file)
    mut_scores = analyze_pose(mut_pose, scorefxn, params['target_chain'], params['binder_chain'], relax_protocol)

    # --- 3. Display Final Results ---
    print("\n" + "#"*60)
    print("              FINAL ENERGY COMPARISON")
    print("#"*60 + "\n")
    print(f"{'Structure':<25} {'Wild-Type (REU)':>20} {'Mutant (REU)':>15}")
    print("-"*60)
    print(f"{'Complex':<25} {wt_scores['complex']:>20.2f} {mut_scores['complex']:>15.2f}")
    print(f"{'Target Chain ('+params['target_chain']+')':<25} {wt_scores['target']:>20.2f} {mut_scores['target']:>15.2f}")
    print(f"{'Binder Chain ('+params['binder_chain']+')':<25} {wt_scores['binder']:>20.2f} {mut_scores['binder']:>15.2f}")
    
    # Clean up the temporary file
    os.remove(mutated_pdb_file)
    print(f"\nCleaned up temporary file: {mutated_pdb_file}")
    
    total_end_time = time.time()
    print(f"\nTotal execution time: {total_end_time - total_start_time:.2f} seconds.")

if __name__ == "__main__":
    main()
