#!/usr/bin/env python
import subprocess
import argparse
import os

def parse_batch_file(batch_file):
    """Parse drug IDs and SMILES from a tab/whitespace separated file."""
    entries = []
    seen = set()
    with open(batch_file, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            if "\t" in raw_line:
                drug_id, smiles = (part.strip() for part in raw_line.split("\t", 1))
            else:
                parts = raw_line.split()
                if len(parts) < 2:
                    continue
                smiles = parts[-1].strip()
                drug_id = " ".join(parts[:-1]).strip()

            if not drug_id or not smiles:
                continue

            key = (drug_id, smiles)
            if key in seen:
                continue
            seen.add(key)
            entries.append(key)
    return entries

def run_docking(drug_id, smiles, pdb_option, pdb_value, chain):
    """Run dockpipe.py with the given arguments."""
    cmd = [
        "python", "molecular_docking.py",
        f"--{pdb_option}", pdb_value,
        "--smiles", smiles,
        "--prefix", drug_id,
        "--chain", chain
    ]
    
    print(f"Running docking for {drug_id}...")
    try:
        subprocess.run(cmd, check=True)
        print(f"Completed docking for {drug_id}")
    except subprocess.CalledProcessError as e:
        print(f"Error running docking for {drug_id}: {e}")

def main():
    parser = argparse.ArgumentParser(description="Batch molecular docking using dockpipe.py")
    
    # Only one of these arguments should be provided
    pdb_group = parser.add_mutually_exclusive_group(required=True)
    pdb_group.add_argument("--pdb", help="PDB identifier")
    pdb_group.add_argument("--local_pdb", help="Path to local PDB file")
    
    parser.add_argument("--batch_file", required=True, help="File containing drug IDs and SMILES strings")
    parser.add_argument("--chain", required=True, help="Chain identifier")
    
    args = parser.parse_args()
    
    # Determine which PDB option to use
    if args.pdb:
        pdb_option = "pdb"
        pdb_value = args.pdb
    else:
        pdb_option = "local_pdb"
        pdb_value = args.local_pdb
    
    # Parse the batch file
    entries = parse_batch_file(args.batch_file)
    
    if not entries:
        print("No valid entries found in batch file.")
        return
    
    print(f"Found {len(entries)} compounds to dock.")
    
    # Run docking for each entry
    for i, (drug_id, smiles) in enumerate(entries, 1):
        print(f"\nProcessing compound {i}/{len(entries)}")
        run_docking(drug_id, smiles, pdb_option, pdb_value, args.chain)

if __name__ == "__main__":
    main()

