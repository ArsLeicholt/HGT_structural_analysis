#!/usr/bin/env python3
"""
Extract residues within 5 Å of ligands from AlphaFold3 predictions.
Generates a single CSV file with structure names and their contact residues.

Handles:
- Multiple ligands (including Mg2+)
- Different ligand types across predictions
- Cases with and without Mg2+
"""

import os
import sys
import csv
import argparse
from pathlib import Path
import numpy as np
from collections import defaultdict


def parse_cif_file(cif_file):
    """Extract protein residues and ligands from CIF file"""
    residues = {}  # (chain, resid) -> {atoms: [...], comp_id: str}
    ligands = defaultdict(list)  # ligand_type -> [atoms]

    with open(cif_file, 'r') as f:
        in_atom_section = False
        for line in f:
            # Check if we're entering the atom section
            if '_atom_site' in line:
                in_atom_section = True
                continue

            if in_atom_section and (line.startswith('ATOM') or line.startswith('HETATM')):
                parts = line.split()
                if len(parts) < 16:
                    continue

                record_type = parts[0]
                atom_id = int(parts[1])
                atom_type = parts[2]
                atom_name = parts[3]
                comp_id = parts[5]  # residue/ligand name
                chain_id = parts[6]
                entity_id = parts[7]
                seq_id = parts[8] if parts[8] != '.' else None
                x = float(parts[10])
                y = float(parts[11])
                z = float(parts[12])
                b_factor = float(parts[15])  # pLDDT

                atom_info = {
                    'record_type': record_type,
                    'atom_id': atom_id,
                    'atom_type': atom_type,
                    'atom_name': atom_name,
                    'comp_id': comp_id,
                    'chain_id': chain_id,
                    'entity_id': entity_id,
                    'seq_id': seq_id,
                    'coords': np.array([x, y, z]),
                    'plddt': b_factor
                }

                # Store protein residues
                if record_type == 'ATOM' and seq_id is not None:
                    res_key = (chain_id, seq_id)
                    if res_key not in residues:
                        residues[res_key] = {
                            'atoms': [],
                            'comp_id': comp_id,
                            'chain_id': chain_id,
                            'seq_id': seq_id
                        }
                    residues[res_key]['atoms'].append(atom_info)

                # Store ligand atoms (including MG and excluding HOH)
                if record_type == 'HETATM' and comp_id != 'HOH':
                    ligands[comp_id].append(atom_info)

            elif in_atom_section and line.startswith('#'):
                # End of atom section
                break

    return {
        'residues': residues,
        'ligands': ligands
    }


def get_contact_residues(residues, ligands, distance_threshold=5.0):
    """Get residues within distance_threshold of any ligand"""
    contact_residues = set()

    # Get all ligand atoms
    all_ligand_atoms = []
    for ligand_type, atoms in ligands.items():
        all_ligand_atoms.extend(atoms)

    if not all_ligand_atoms:
        return []

    # For each residue, check if any atom is within threshold
    for res_key, res_data in residues.items():
        chain_id, seq_id = res_key
        res_comp_id = res_data['comp_id']

        is_contact = False
        for res_atom in res_data['atoms']:
            for lig_atom in all_ligand_atoms:
                distance = np.linalg.norm(res_atom['coords'] - lig_atom['coords'])
                if distance <= distance_threshold:
                    is_contact = True
                    break
            if is_contact:
                break

        if is_contact:
            # Format: chain:number:type (e.g., "A:37:ARG")
            residue_id = f"{chain_id}:{seq_id}:{res_comp_id}"
            contact_residues.add(residue_id)

    # Sort residues by chain and number
    sorted_residues = sorted(contact_residues, key=lambda x: (x.split(':')[0], int(x.split(':')[1])))

    return sorted_residues


def process_prediction(result_dir):
    """Process one AF3 prediction and extract contact residues"""
    result_path = Path(result_dir)
    original_name = result_path.name

    # Find best model from ranking_scores.csv
    ranking_file = result_path / "ranking_scores.csv"
    best_sample_dir = None

    if ranking_file.exists():
        with open(ranking_file, 'r') as f:
            reader = csv.DictReader(f)
            best_score = -1
            best_seed = None
            best_sample = None
            for row in reader:
                score = float(row['ranking_score'])
                if score > best_score:
                    best_score = score
                    best_seed = row['seed']
                    best_sample = row['sample']

        if best_seed and best_sample:
            best_sample_dir = result_path / f"seed-{best_seed}_sample-{best_sample}"

    # Fallback to first sample if ranking not found
    if not best_sample_dir or not best_sample_dir.exists():
        sample_dirs = sorted([d for d in result_path.iterdir()
                             if d.is_dir() and d.name.startswith('seed-')])
        if not sample_dirs:
            return None, None
        best_sample_dir = sample_dirs[0]

    # Look for CIF file
    cif_file = best_sample_dir / "model.cif"
    if not cif_file.exists():
        return None, None

    # Parse CIF file
    data = parse_cif_file(cif_file)

    if not data['residues'] or not data['ligands']:
        return None, None

    # Get contact residues
    contact_residues = get_contact_residues(data['residues'], data['ligands'])

    if not contact_residues:
        return None, None

    return original_name, contact_residues


def process_all_predictions(base_dir, output_file, distance_threshold=5.0):
    """Process all predictions and create single CSV"""
    base_path = Path(base_dir)

    if not base_path.exists():
        print(f"Error: Directory {base_dir} does not exist")
        return

    # Find all subdirectories
    result_dirs = [d for d in base_path.iterdir() if d.is_dir()]

    if not result_dirs:
        print(f"Warning: No subdirectories found in {base_dir}")
        return

    print(f"Processing {len(result_dirs)} predictions from {base_dir}")
    print(f"Output file: {output_file}")
    print(f"Distance threshold: {distance_threshold} Å")
    print("-" * 80)

    all_results = []
    successful = 0
    failed = 0

    for result_dir in sorted(result_dirs):
        # Find nested directory if it exists
        nested_dirs = [d for d in result_dir.iterdir()
                      if d.is_dir() and not d.name.startswith('.')]

        if nested_dirs:
            # Process nested directory
            for nested_dir in nested_dirs[:1]:  # Just use first one
                structure_name, contact_residues = process_prediction(nested_dir)
                if structure_name and contact_residues:
                    # Store as comma-separated string
                    residue_string = ", ".join(contact_residues)
                    all_results.append({
                        'structure': structure_name,
                        'contact_residues': residue_string
                    })
                    print(f"✓ {result_dir.name}: {len(contact_residues)} residues")
                    successful += 1
                else:
                    print(f"✗ {result_dir.name}: Failed")
                    failed += 1
        else:
            # Process directory directly
            structure_name, contact_residues = process_prediction(result_dir)
            if structure_name and contact_residues:
                residue_string = ", ".join(contact_residues)
                all_results.append({
                    'structure': structure_name,
                    'contact_residues': residue_string
                })
                print(f"✓ {result_dir.name}: {len(contact_residues)} residues")
                successful += 1
            else:
                print(f"✗ {result_dir.name}: Failed")
                failed += 1

    # Write single CSV file
    if all_results:
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=['structure', 'contact_residues'])
            writer.writeheader()
            writer.writerows(all_results)

        print("-" * 80)
        print(f"Complete! Processed {successful + failed} predictions")
        print(f"  Successful: {successful}")
        print(f"  Failed: {failed}")
        print(f"\nOutput saved to: {output_file}")
    else:
        print("-" * 80)
        print("No results to write!")


def main():
    parser = argparse.ArgumentParser(
        description="Extract residues within specified distance of ligands from AF3 predictions",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all predictions in af3_output
  python3 extract_ligand_contacts.py --input af3_output --output contacts.csv

  # Use custom distance threshold
  python3 extract_ligand_contacts.py --input af3_output_2 --output contacts_5A.csv --distance 5.0

Output format:
  CSV with two columns:
  - structure: Name of the prediction
  - contact_residues: Comma-separated list of residues in format "chain:number:type"
    Example: "A:37:ARG, A:38:GLY, B:53:ARG"
        """
    )

    parser.add_argument("--input", "-i", required=True,
                       help="Input directory containing AF3 prediction subdirectories")
    parser.add_argument("--output", "-o", required=True,
                       help="Output CSV file")
    parser.add_argument("--distance", "-d", type=float, default=5.0,
                       help="Distance threshold in Angstroms (default: 5.0)")

    args = parser.parse_args()

    # Validate distance
    if args.distance <= 0:
        print("Error: Distance must be positive")
        sys.exit(1)

    # Process predictions
    process_all_predictions(args.input, args.output, args.distance)


if __name__ == "__main__":
    main()
