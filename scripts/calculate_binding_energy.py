#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import time
from pyrosetta import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.relax import FastRelax

def relax_pose(pose, scorefxn):
    """
    Performs a structural relaxation on the pose to relieve clashes.
    """
    print("\n" + "="*60)
    print("--- Starting Structural Relaxation (FastRelax) ---")
    
    fast_relax = FastRelax(scorefxn, 5) # 5 rounds is standard
    relaxed_pose = pose.clone()
    
    score_before = scorefxn(relaxed_pose)
    fast_relax.apply(relaxed_pose)
    score_after = scorefxn(relaxed_pose)
    
    print("Relaxation complete.")
    print(f"Score before relaxation: {score_before:.2f} REU")
    print(f"Score after relaxation:  {score_after:.2f} REU")
    print("="*60)
    
    return relaxed_pose

def calculate_binding_metrics(pose, chain_a, chain_b, scorefxn):
    """
    Uses the InterfaceAnalyzerMover to calculate detailed binding metrics.
    """
    print("\n" + "="*60)
    print("--- Calculating Binding Energy with InterfaceAnalyzerMover ---")

    # Define the interface between the two chains
    interface_definition = f"{chain_a}_{chain_b}"
    print(f"Defining interface as: {interface_definition}")

    # Setup the InterfaceAnalyzerMover
    iam = InterfaceAnalyzerMover(interface_definition)
    iam.set_compute_packstat(True)
    iam.set_pack_separated(True)
    iam.set_scorefunction(scorefxn)

    print("Running InterfaceAnalyzerMover...")
    iam.apply(pose)

    print("\n--- Interface Analysis Results ---")
    # This prints a detailed summary of all metrics
    iam.show()
    
    # The key metric for binding energy is dG
    binding_energy = iam.get_interface_dG()
    
    print("\n" + "="*60)
    print(f"Final Calculated Binding Energy (dG_separated): {binding_energy:.2f} REU")
    print("A more negative value indicates stronger binding.")
    print("="*60)
    
def get_energy_breakdown(pose, binder_chain, scorefxn):
    """
    Calculates the energy breakdown by score term for the binder chain.
    """
    print("\n" + "="*60)
    print(f"--- Calculating Energy Breakdown for Binder Chain: {binder_chain} ---")

    # Score the pose to populate the energy graph
    scorefxn(pose)
    
    pdb_info = pose.pdb_info()
    binder_residues = [i for i in range(1, pose.total_residue() + 1) if pdb_info.chain(i) == binder_chain]

    if not binder_residues:
        print(f"Error: Could not find any residues for binder chain '{binder_chain}'")
        return

    energies = pose.energies()
    energy_breakdown = {}

    # Sum energies for all residues in the binder chain
    for res_idx in binder_residues:
        residue_energies = energies.residue_total_energies(res_idx)
        for score_term_index in range(1, int(rosetta.core.scoring.ScoreType.n_score_types) + 1):
            score_term = rosetta.core.scoring.ScoreType(score_term_index)
            energy_value = residue_energies[score_term] * scorefxn.weights()[score_term]
            
            if energy_value != 0:
                term_name = score_term.name
                energy_breakdown[term_name] = energy_breakdown.get(term_name, 0.0) + energy_value
    
    print(f"Energy Breakdown for Binder Chain '{binder_chain}':")
    print("-" * 45)
    print(f"{'Score Term':<20} {'Weighted Energy (REU)':>20}")
    print("-" * 45)
    # Sort by energy value for clarity
    sorted_breakdown = sorted(energy_breakdown.items(), key=lambda item: item[1])
    for term, energy in sorted_breakdown:
        print(f"{term:<20} {energy:>20.4f}")
    print("="*60)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate binding energy and interface metrics for a protein complex using PyRosetta's InterfaceAnalyzerMover."
    )
    parser.add_argument("pdb_file", type=str, help="Path to the input PDB file.")
    parser.add_argument("--chain_a", required=True, type=str, help="Identifier for the first chain (e.g., target).")
    parser.add_argument("--chain_b", required=True, type=str, help="Identifier for the second chain (e.g., binder).")
    parser.add_argument("--relax", action="store_true", help="Perform a structural relaxation (FastRelax) before analysis.")
    parser.add_argument("--breakdown", action="store_true", help="Show a detailed energy breakdown for the binder chain (chain_b).")

    args = parser.parse_args()

    start_time = time.time()

    if not os.path.isfile(args.pdb_file):
        print(f"Error: PDB file not found at '{args.pdb_file}'")
    else:
        # Initialize PyRosetta with options to silence verbose output
        init(extra_options="-ignore_unrecognized_res -ex1 -ex2 -out:level 100")
        
        # Create a standard score function (ref2015 is default)
        scorefxn = get_score_function(True)
        
        print(f"Loading PDB: {args.pdb_file}")
        pose = pose_from_pdb(args.pdb_file)
        
        analysis_pose = pose
        if args.relax:
            analysis_pose = relax_pose(pose, scorefxn)
            
        calculate_binding_metrics(analysis_pose, args.chain_a, args.chain_b, scorefxn)
        
        if args.breakdown:
            get_energy_breakdown(analysis_pose, args.chain_b, scorefxn)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("\n" + "="*60)
    print(f"Total execution time: {elapsed_time:.2f} seconds")
    print("="*60)

