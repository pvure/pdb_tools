# PDB Tools

This repository contains a collection of Python scripts for common PDB file manipulations.

## `separate_structures.py`

Separates a multi-chain PDB file into individual PDB files for each specified chain.

### Usage
python separate_structures.py /path/to/your/protein.pdb

**Example Output:**
If `protein.pdb` contains chains A and B, this script will generate `protein_chain_A.pdb` and `protein_chain_B.pdb`.

## `mutate_single_residue.py`

Performs an in-silico mutation by renaming a specified residue in a PDB file and saves the result as a new file.

### Usage

The script requires the input file, chain, target residue ID, original residue name, and new residue name.

python mutate_single_residue.py <pdb_file> --chain <ID> --residue_id <num> --original_resname <AAA> --new_resname <BBB>

**Example:**
To mutate Tryptophan (TRP) at position 191 on chain B to Serine (SER) in a file named `my_protein.pdb`:
python mutate_single_residue.py my_protein.pdb --chain B --residue_id 191 --original_resname TRP --new_resname SER

**Example Output:**
This will create a new file named `my_protein_chainB_TRP191SER.pdb`.

The script includes a `tolerance` feature to find residues even if the numbering is slightly off.