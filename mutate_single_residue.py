import argparse
import os

def find_target_residue(pdb_file_path, chain_id, target_residue_id, original_resname, tolerance):
    """
    Scans a PDB file to find the best matching residue ID within a tolerance.

    Args:
        pdb_file_path (str): Path to the PDB file.
        chain_id (str): The chain to search in.
        target_residue_id (int): The desired residue number.
        original_resname (str): The 3-letter code of the residue we expect to find.
        tolerance (int): The +/- range to search within.

    Returns:
        int or None: The actual residue ID of the best match, or None if not found.
    """
    best_match = {'id': None, 'distance': float('inf')}
    
    seen_residues = set()

    with open(pdb_file_path, 'r') as infile:
        for line in infile:
            if line.startswith('ATOM'):
                line_chain_id = line[21]
                if line_chain_id != chain_id:
                    continue

                try:
                    line_residue_id = int(line[22:26].strip())
                    # unique identifier for a residue is number and name
                    residue_key = (line_residue_id, line[17:20].strip())
                except ValueError:
                    continue
                
                # skip over already seen residues
                if residue_key in seen_residues:
                    continue
                seen_residues.add(residue_key)

                line_resname = residue_key[1]
                
                if line_resname == original_resname:
                    distance = abs(line_residue_id - target_residue_id)
                    if distance <= tolerance and distance < best_match['distance']:
                        best_match['id'] = line_residue_id
                        best_match['distance'] = distance

    return best_match['id']


def mutate_residue(pdb_file_path, chain_id, residue_id, new_resname, original_resname, tolerance, output_file_path=None):
    """
    Reads a PDB file, mutates a single specified residue by changing its name,
    and writes the result to a new PDB file. It allows for a small tolerance
    in the residue ID.
    """
    if not os.path.isfile(pdb_file_path):
        print(f"Error: Input file not found at '{pdb_file_path}'")
        return

    # standardize residue names to uppercase
    new_resname = new_resname.upper()
    original_resname = original_resname.upper()

    # find residue to mutate using 1 aa tolerance
    actual_residue_to_mutate = find_target_residue(
        pdb_file_path, chain_id, residue_id, original_resname, tolerance
    )

    if actual_residue_to_mutate is None:
        print(f"Error: Could not find residue '{original_resname}' on chain {chain_id} "
              f"near position {residue_id} (tolerance: {tolerance}).")
        return
        
    if actual_residue_to_mutate != residue_id:
        print(f"Note: Target residue {residue_id} not found. "
              f"Found matching residue '{original_resname}' at position {actual_residue_to_mutate}. Proceeding with mutation.")
    
    # generate output path
    if output_file_path is None:
        directory = os.path.dirname(pdb_file_path)
        base_name = os.path.splitext(os.path.basename(pdb_file_path))[0]
        output_file_path = os.path.join(
            directory,
            f"{base_name}_chain{chain_id}_{original_resname}{actual_residue_to_mutate}{new_resname}.pdb"
        )

    print(f"\nReading from: {pdb_file_path}")
    print(f"Mutating Chain {chain_id}, Residue {actual_residue_to_mutate} ({original_resname}) to {new_resname}")

    # write new PDB file with mutation
    try:
        with open(pdb_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
            for line in infile:
                if line.startswith('ATOM'):
                    line_chain_id = line[21]
                    try:
                        line_residue_id = int(line[22:26].strip())
                    except ValueError:
                        outfile.write(line)
                        continue

                    if line_chain_id == chain_id and line_residue_id == actual_residue_to_mutate:
                        modified_line = line[:17] + new_resname.ljust(3) + line[20:]
                        outfile.write(modified_line)
                    else:
                        outfile.write(line)
                else:
                    outfile.write(line)
        
        print(f"Successfully created mutated file: {output_file_path}")

    except IOError as e:
        print(f"An error occurred during file processing: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform a single-point mutation on a residue in a PDB file with tolerance."
    )
    parser.add_argument("pdb_file", type=str, help="Path to the input PDB file.")
    parser.add_argument("--chain", required=True, type=str, help="Chain ID (e.g., A).")
    parser.add_argument("--residue_id", required=True, type=int, help="Residue number to target (e.g., 142).")
    parser.add_argument("--original_resname", required=True, type=str, help="3-letter code of the original residue (e.g., SER).")
    parser.add_argument("--new_resname", required=True, type=str, help="3-letter code for the new residue (e.g., TRP).")
    parser.add_argument("--tolerance", type=int, default=1, help="Search tolerance (+/-) for the residue ID. Default is 1.")
    parser.add_argument("--output_file", type=str, help="Optional: Path for the output PDB file.")

    args = parser.parse_args()
    
    mutate_residue(
        args.pdb_file,
        args.chain,
        args.residue_id,
        args.new_resname,
        args.original_resname,
        args.tolerance,
        args.output_file
    )

