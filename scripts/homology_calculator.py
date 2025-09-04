#!/usr/bin/env python3
"""
PDB Chain B Homology Analysis using MMseqs2

This script extracts protein sequences from chain B of PDB files,
then uses the MMseqs2 command-line tool to perform a rapid all-vs-all 
sequence search and calculate a pairwise identity matrix.

This version is compatible with older Biopython versions.

Requirements:
- MMseqs2: conda install -c bioconda mmseqs2
- Python libs: pip install biopython matplotlib seaborn numpy pandas

Usage:
python homology_analyzer_mmseqs2.py /path/to/pdb/folder
"""
import os
import sys
import argparse
import subprocess
import shutil
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import PDB, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# This specific import is used for backward compatibility with older Biopython versions
from Bio.PDB.Polypeptide import protein_letters_3to1

class PDBSequenceAnalyzer:
    def __init__(self, pdb_folder):
        self.pdb_folder = Path(pdb_folder)
        self.sequences = {}
        self.pairwise_identities = None

    def extract_chain_b_sequences(self):
        print("--- Step 1: Extracting Chain B sequences ---")
        parser = PDB.PDBParser(QUIET=True)
        
        for pdb_file in self.pdb_folder.glob("*.pdb"):
            try:
                structure = parser.get_structure(pdb_file.stem, pdb_file)
                chain_b = next((model['B'] for model in structure if 'B' in model), None)
                
                if chain_b is None:
                    print(f"  Warning: Chain B not found in {pdb_file.name}")
                    continue
                
                # This .join() with a generator expression is an efficient way to build the sequence string.
                # It uses the imported 'protein_letters_3to1' dictionary to perform the conversion,
                # which avoids the 'seq1' error.
                sequence = "".join(
                    protein_letters_3to1.get(res.get_resname(), 'X')
                    for res in chain_b.get_residues()
                    if PDB.is_aa(res, standard=True)
                )
                
                if sequence:
                    # Sanitize sequence ID for MMseqs2 (no special characters)
                    clean_name = pdb_file.stem.replace(" ", "_")
                    self.sequences[clean_name] = sequence
                    print(f"  Extracted sequence from {pdb_file.name} ({len(sequence)} residues)")
                else:
                    print(f"  Warning: No standard amino acids in chain B of {pdb_file.name}")
                    
            except Exception as e:
                print(f"  Error processing {pdb_file.name}: {e}")
        
        print(f"\nSuccessfully extracted {len(self.sequences)} sequences.")
        return self.sequences

    def calculate_identities_with_mmseqs2(self):
        """
        Uses MMseqs2 to perform an all-vs-all search and creates an identity matrix.
        """
        print("\n--- Step 2: Calculating Pairwise Identities with MMseqs2 ---")
        if len(self.sequences) < 2:
            print("Comparison requires at least 2 sequences. Skipping.")
            return None

        # Create a temporary directory for MMseqs2 files
        tmp_dir = Path("./mmseqs_tmp")
        tmp_dir.mkdir(exist_ok=True)
        
        input_fasta = tmp_dir / "input_sequences.fasta"
        
        # Write extracted sequences to the FASTA file
        records = [SeqRecord(Seq(seq), id=name, description="") for name, seq in self.sequences.items()]
        SeqIO.write(records, input_fasta, "fasta")
        
        try:
            print("  Running MMseqs2 for all-vs-all comparison...")
            db_path = tmp_dir / "sequenceDB"
            result_path = tmp_dir / "search_results"
            
            # MMseqs2 command pipeline
            subprocess.run(['mmseqs', 'createdb', str(input_fasta), str(db_path)], check=True, capture_output=True, text=True)
            subprocess.run(['mmseqs', 'search', str(db_path), str(db_path), str(result_path), str(tmp_dir)], check=True, capture_output=True, text=True)
            output_tsv = tmp_dir / "results.tsv"
            subprocess.run(['mmseqs', 'convertalis', str(db_path), str(db_path), str(result_path), str(output_tsv), '--format-output', 'query,target,pident'], check=True, capture_output=True, text=True)

            # Read the results into a pandas DataFrame
            results_df = pd.read_csv(output_tsv, sep='\t', header=None, names=['query', 'target', 'pident'])
            
            # Build the identity matrix
            names = list(self.sequences.keys())
            identity_matrix = pd.DataFrame(np.identity(len(names)) * 100, index=names, columns=names)
            
            for _, row in results_df.iterrows():
                identity_matrix.loc[row['query'], row['target']] = row['pident'] 
                identity_matrix.loc[row['target'], row['query']] = row['pident'] 
            
            self.pairwise_identities = identity_matrix
            print("  Pairwise identity matrix calculated successfully.")

        except FileNotFoundError:
            print("\nError: 'mmseqs' command not found.", file=sys.stderr)
            print("Please make sure MMseqs2 is installed and in your system's PATH.", file=sys.stderr)
            sys.exit(1)
        except subprocess.CalledProcessError as e:
            print(f"\nError running MMseqs2: {e}", file=sys.stderr)
            print(f"STDERR:\n{e.stderr}", file=sys.stderr)
            sys.exit(1)
        finally:
            # Clean up the temporary directory
            shutil.rmtree(tmp_dir)

        return self.pairwise_identities
    
    def generate_plots(self):
        """Generates and saves all plots."""
        print("\n--- Step 3: Generating Plots ---")
        if self.pairwise_identities is not None:
            self.plot_identity_heatmap()
        else:
            print("Skipping plots as no identity matrix is available.")

    def plot_identity_heatmap(self):
        """Plots the pairwise identity heatmap."""
        plt.figure(figsize=(12, 10))
        sns.heatmap(self.pairwise_identities, annot=True, fmt='.1f', cmap='viridis', cbar_kws={'label': 'Sequence Identity (%)'})
        plt.title('Pairwise Sequence Identity Matrix (via MMseqs2)')
        plt.tight_layout()
        plt.savefig('pairwise_identity_heatmap_mmseqs2.png', dpi=300)
        print("  Saved heatmap to 'pairwise_identity_heatmap_mmseqs2.png'")
        plt.close()
    
    def save_results(self):
        """Saves identity matrix to a file."""
        print("\n--- Step 4: Saving Results ---")
        if self.pairwise_identities is not None:
            self.pairwise_identities.to_csv("pairwise_identities_mmseqs2.csv")
            print("  Saved identity matrix to 'pairwise_identities_mmseqs2.csv'")

def main():
    parser = argparse.ArgumentParser(description='Analyze Chain B sequences from PDB files using MMseqs2.')
    parser.add_argument('pdb_folder', help='Path to the folder containing PDB files.')
    args = parser.parse_args()
    
    pdb_path = Path(args.pdb_folder)
    if not pdb_path.is_dir():
        print(f"Error: Folder '{args.pdb_folder}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    analyzer = PDBSequenceAnalyzer(pdb_path)
    analyzer.extract_chain_b_sequences()
    analyzer.calculate_identities_with_mmseqs2()
    analyzer.generate_plots()
    analyzer.save_results()
    
    print("\n Analysis complete!")

if __name__ == "__main__":
    main()