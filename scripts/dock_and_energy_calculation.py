#!/usr/bin/env python

from __future__ import print_function
import argparse
import os
import time
import tempfile
from pyrosetta import *
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.protocols.docking import DockingProtocol

def combine_pdbs(target_path, binder_path, output_path):
    """
    Combines two PDB files into a single file, re-chaining them to A and B.
    """
    print(f"Combining {os.path.basename(target_path)} (Chain A) and {os.path.basename(binder_path)} (Chain B)...")
    with open(output_path, 'w') as outfile:
        with open(target_path, 'r') as infile:
            for line in infile:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    outfile.write(line[:21] + 'A' + line[22:])
        outfile.write("TER\n")
        with open(binder_path, 'r') as infile:
            for line in infile:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    outfile.write(line[:21] + 'B' + line[22:])
    print(f"Successfully created combined PDB at {output_path}")

def run_docking_protocol(pose, scorefxn_low, scorefxn_high):
    """
    Performs a full Rosetta docking simulation on the input pose.
    """
    docking_protocol = DockingProtocol()
    docking_protocol.set_lowres_scorefxn(scorefxn_low)
    docking_protocol.set_highres_scorefxn(scorefxn_high)
    
    # This sets up the FoldTree for docking between chains A and B
    protocols.docking.setup_foldtree(pose, "A_B", Vector1([1]))
    
    docked_pose = pose.clone()
    docking_protocol.apply(docked_pose)
    
    return docked_pose

def calculate_binding_metrics(pose, chain_a, chain_b, scorefxn):
    """
    Uses the InterfaceAnalyzerMover to calculate detailed binding metrics.
    """
    print("\n" + "="*60)
    print(f"--- Analyzing Best Docked Pose ---")

    interface_definition = f"{chain_a}_{chain_b}"
    iam = InterfaceAnalyzerMover(interface_definition)
    iam.set_compute_packstat(True)
    iam.set_pack_separated(True)
    iam.set_scorefunction(scorefxn)

    iam.apply(pose)

    print("\n--- Interface Analysis Results ---")
    iam.show()
    
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
    scorefxn(pose)
    
    pdb_info = pose.pdb_info()
    binder_residues = [i for i in range(1, pose.total_residue() + 1) if pdb_info.chain(i) == binder_chain]

    if not binder_residues:
        print(f"Error: Could not find any residues for binder chain '{binder_chain}'")
        return

    energies = pose.energies()
    energy_breakdown = {}

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
    sorted_breakdown = sorted(energy_breakdown.items(), key=lambda item: item[1])
    for term, energy in sorted_breakdown:
        print(f"{term:<20} {energy:>20.4f}")
    print("="*60)

def main():
    parser = argparse.ArgumentParser(
        description="Combine, dock, and analyze a protein complex from separate PDBs using the Rosetta DockingProtocol."
    )
    parser.add_argument("target_pdb", type_str, help="Path to the target PDB file (will become chain A).")
    parser.add_argument("binder_pdb", type_str, help="Path to the binder PDB file (will become chain B).")
    parser.add_argument("--n_decoys", type=int, default=5, help="Number of docking simulations (decoys) to generate.")
    parser.add_argument("--output_pdb", type=str, help="Optional: Path to save the best-scoring docked complex PDB file.")
    parser.add_argument("--breakdown", action="store_true", help="Show a detailed energy breakdown for the binder chain (B) of the best decoy.")
    
    args = parser.parse_args()
    start_time = time.time()

    temp_fd, temp_pdb_path = tempfile.mkstemp(suffix=".pdb")
    os.close(temp_fd)

    try:
        combine_pdbs(args.target_pdb, args.binder_pdb, temp_pdb_path)

        init(extra_options="-ignore_unrecognized_res -ex1 -ex2 -out:level 100")
        
        # Score functions for docking
        scorefxn_high = create_score_function('ref2015') # Standard full-atom
        scorefxn_low = create_score_function('interchain_cen') # Low-res centroid
        
        print(f"\nLoading combined PDB into PyRosetta...")
        initial_pose = pose_from_pdb(temp_pdb_path)
        
        best_score = float('inf')
        best_pose = None

        print("\n" + "="*60)
        print(f"--- Starting Docking Simulation ({args.n_decoys} decoys) ---")
        docking_start_time = time.time()
        for i in range(args.n_decoys):
            print(f"Generating decoy {i+1}/{args.n_decoys}...")
            docked_pose = run_docking_protocol(initial_pose, scorefxn_low, scorefxn_high)
            current_score = scorefxn_high(docked_pose)
            print(f"  Decoy {i+1} score: {current_score:.2f} REU")
            
            if current_score < best_score:
                best_score = current_score
                best_pose = docked_pose
                print(f"  *** New best score found! ***")
        
        docking_end_time = time.time()
        docking_elapsed = docking_end_time - docking_start_time
        print(f"Docking simulation finished in {docking_elapsed:.2f} seconds.")
        print("="*60)

        if best_pose:
            print(f"\nDocking complete. Best score: {best_score:.2f} REU")
            
            analysis_start_time = time.time()
            calculate_binding_metrics(best_pose, 'A', 'B', scorefxn_high)
            
            if args.breakdown:
                get_energy_breakdown(best_pose, 'B', scorefxn_high)
            
            analysis_end_time = time.time()
            analysis_elapsed = analysis_end_time - analysis_start_time
            print(f"Final analysis finished in {analysis_elapsed:.2f} seconds.")
                
            if args.output_pdb:
                best_pose.dump_pdb(args.output_pdb)
                print(f"\nSuccessfully saved best decoy to {args.output_pdb}")
        else:
            print("\nError: Docking simulation failed to produce a valid pose.")

    finally:
        os.remove(temp_pdb_path)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("\n" + "="*60)
    print(f"Total execution time: {elapsed_time:.2f} seconds")
    print("="*60)

if __name__ == "__main__":
    main()

