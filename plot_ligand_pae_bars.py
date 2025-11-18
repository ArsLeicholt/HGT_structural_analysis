#!/usr/bin/env python3
"""
Create bar plots of ligand PAE values for structures in visualization_examples.
Bar plots are grey with ligand names on x-axis.
"""

import json
import os
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


def parse_cif_file(cif_file):
    """
    Parse CIF file to extract atom positions and B-factors (pLDDT).

    Returns:
        dict: Contains 'atoms' list with atom info and 'ligands'
    """
    atoms = []
    ligands = defaultdict(list)

    with open(cif_file, 'r') as f:
        in_atom_section = False
        for line in f:
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
                b_factor = float(parts[15])  # pLDDT stored as B-factor

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

                atoms.append(atom_info)

                # Store ligand atoms (non-protein, excluding MG and water)
                if record_type == 'HETATM' and comp_id not in ['MG', 'HOH', 'WAT', 'DOD']:
                    ligands[comp_id].append(atom_info)

            elif in_atom_section and line.startswith('#'):
                break

    return {
        'atoms': atoms,
        'ligands': ligands
    }


def identify_ligand_chains(cif_file):
    """Identify which chains contain substrate ligands."""
    data = parse_cif_file(cif_file)

    if not data['ligands']:
        return set()

    ligand_chains = set()
    for ligand_name, ligand_atoms in data['ligands'].items():
        for atom in ligand_atoms:
            ligand_chains.add(atom['chain_id'])

    return ligand_chains


def calculate_ligand_pae(confidences_file, cif_file, distance_threshold=4.0, min_contacts=3):
    """
    Calculate average PAE between binding pocket residues and substrate ligand.
    Lower PAE = better, so we invert it to higher = better.

    Returns:
        tuple: (mean_pae, sem_pae, n_values) or (None, None, None) if no data
        sem_pae = standard error of the mean = std / sqrt(n)
    """
    data = parse_cif_file(cif_file)

    if not data['ligands']:
        return None, None, None

    # Get all substrate ligand atoms
    all_ligand_atoms = []
    for ligand_name, ligand_atoms in data['ligands'].items():
        all_ligand_atoms.extend(ligand_atoms)

    if not all_ligand_atoms:
        return None, None, None

    # Identify ligand chains
    ligand_chains = identify_ligand_chains(cif_file)
    if not ligand_chains:
        return None, None, None

    # Load full PAE matrix
    try:
        with open(confidences_file, 'r') as f:
            conf_data = json.load(f)

        pae_matrix = np.array(conf_data['pae'])
        token_chain_ids = conf_data['token_chain_ids']
        token_res_ids = conf_data.get('token_res_ids', [None] * len(token_chain_ids))
    except (KeyError, FileNotFoundError):
        return None, None, None

    # Identify binding pocket residues
    residue_info = {}
    for atom in data['atoms']:
        if atom['record_type'] == 'ATOM' and atom['seq_id'] is not None:
            res_key = (atom['chain_id'], atom['seq_id'])

            is_contact = False
            for ligand_atom in all_ligand_atoms:
                dist = np.linalg.norm(atom['coords'] - ligand_atom['coords'])
                if dist <= distance_threshold:
                    is_contact = True
                    break

            if res_key not in residue_info:
                residue_info[res_key] = []
            residue_info[res_key].append({
                'atom_type': atom['atom_type'],
                'is_contact': is_contact
            })

    # Find binding pocket residues
    binding_pocket_residues = set()
    for res_key, atoms in residue_info.items():
        heavy_atom_contacts = sum(1 for a in atoms
                                   if a['is_contact'] and a['atom_type'] != 'H')
        if heavy_atom_contacts >= min_contacts:
            binding_pocket_residues.add(res_key)

    if not binding_pocket_residues:
        return None, None, None

    # Map binding pocket residues to token indices
    pocket_token_indices = []
    for i, (chain_id, res_id) in enumerate(zip(token_chain_ids, token_res_ids)):
        if res_id is not None:
            res_key = (chain_id, str(res_id))
            if res_key in binding_pocket_residues:
                pocket_token_indices.append(i)

    # Map ligand chains to token indices
    ligand_token_indices = []
    for i, chain_id in enumerate(token_chain_ids):
        if chain_id in ligand_chains:
            ligand_token_indices.append(i)

    if not pocket_token_indices or not ligand_token_indices:
        return None, None, None

    # Extract PAE values between binding pocket and ligand
    pae_values = []
    for pocket_idx in pocket_token_indices:
        for lig_idx in ligand_token_indices:
            pae_val = pae_matrix[pocket_idx, lig_idx]
            if not np.isnan(pae_val):
                pae_values.append(pae_val)

    if not pae_values:
        return None, None, None

    # Return mean PAE, standard error of the mean (SEM), and number of values
    mean_pae = np.mean(pae_values)
    std_pae = np.std(pae_values)
    n = len(pae_values)
    sem_pae = std_pae / np.sqrt(n)  # Standard error of the mean

    return mean_pae, sem_pae, n


def get_ligand_from_filename(filename):
    """Extract ligand name from filename like GAGA-0515.XGPRT_gene_CIT_dimer.cif"""
    # Pattern: extract the part after the last underscore (before .cif or _dimer.cif)
    # For monomer enzymes: GAGA-0515.XGPRT_gene_PurF.cif -> PurF
    # For dimer ligands: GAGA-0515.XGPRT_gene_CIT_dimer.cif -> CIT
    match = re.search(r'_([A-Za-z0-9]+?)(?:_dimer)?\.cif$', filename)
    if match:
        return match.group(1)
    return None


def main():
    """Main function to create bar plots of ligand PAE values."""

    print("Analyzing ligand PAE values from visualization_examples")
    print("=" * 80)

    viz_dir = Path("visualization_examples")

    if not viz_dir.exists():
        print(f"Error: {viz_dir} directory not found")
        return

    # Collect PAE values for each structure
    ligand_pae_data = []

    for cif_file in sorted(viz_dir.glob("*.cif")):
        # Get corresponding confidence file from original directories
        base_name = cif_file.stem  # e.g., GAGA-0515.XGPRT_gene_CIT_dimer

        # Try to find the source structure in af3_output directories
        confidences_file = None
        for output_dir in ["af3_output", "af3_output_2", "af3_output_3"]:
            # Check various possible locations - try both with and without _dimer suffix
            potential_names = [
                base_name,
                base_name + "_dimer",  # For monomers that came from dimers
                base_name.lower(),
                (base_name + "_dimer").lower(),
            ]

            potential_dirs = [Path(output_dir) / name for name in potential_names]

            for pot_dir in potential_dirs:
                if pot_dir.exists():
                    # Look for nested directory
                    nested_dirs = [d for d in pot_dir.iterdir()
                                   if d.is_dir() and not d.name.startswith('.')
                                   and not d.name.startswith('seed-')]

                    if nested_dirs:
                        nested_dir = nested_dirs[0]
                        conf_file = nested_dir / f"{nested_dir.name}_confidences.json"
                        if conf_file.exists():
                            confidences_file = conf_file
                            break

            if confidences_file:
                break

        if not confidences_file:
            print(f"Warning: Could not find confidences file for {cif_file.name}")
            continue

        # Calculate ligand PAE (returns mean, SEM, and n)
        ligand_pae_mean, ligand_pae_sem, n_values = calculate_ligand_pae(confidences_file, cif_file)

        if ligand_pae_mean is not None:
            ligand_name = get_ligand_from_filename(cif_file.name)
            if ligand_name:
                ligand_pae_data.append({
                    'ligand': ligand_name,
                    'pae': ligand_pae_mean,
                    'pae_sem': ligand_pae_sem,
                    'n': n_values,
                    'filename': cif_file.name
                })
                print(f"  {ligand_name}: PAE = {ligand_pae_mean:.2f} ± {ligand_pae_sem:.2f} Å (n={n_values})")
            else:
                print(f"  Warning: Could not extract ligand name from {cif_file.name}")
        else:
            print(f"  Warning: Could not calculate PAE for {cif_file.name}")

    if not ligand_pae_data:
        print("\nNo PAE data collected. Exiting.")
        return

    # Define functional groups and their order
    group_order = {
        'PRTases': ['APRT', 'HGPRT', 'HPRT', 'OPRT', 'UPRT', 'XGPRT', 'PRPP'],
        'TCA + α-oxo acids': ['AKG', 'CIT', 'CS', 'GXO', 'HPV', 'ME', 'MLT', 'OAA', 'SCS'],
        'Nucleotide/phosphate transfer': ['AK', 'NDK', 'PEPCK', 'PGK', 'PK', 'PurF']
    }

    # Create ordered list of ligands based on groups
    ordered_ligands = []
    for group_name, ligands_in_group in group_order.items():
        ordered_ligands.extend(ligands_in_group)

    # Sort ligand_pae_data according to the defined order
    ligand_pae_dict = {d['ligand']: d for d in ligand_pae_data}
    ligand_pae_data_ordered = []
    for lig in ordered_ligands:
        if lig in ligand_pae_dict:
            ligand_pae_data_ordered.append(ligand_pae_dict[lig])

    # Create bar plot
    print("\n" + "=" * 80)
    print("Creating bar plot...")

    ligands = [d['ligand'] for d in ligand_pae_data_ordered]
    pae_values = [d['pae'] for d in ligand_pae_data_ordered]
    pae_sem_values = [d['pae_sem'] for d in ligand_pae_data_ordered]

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))

    # Assign colors based on PAE thresholds
    # <5 Å: blue (high confidence)
    # 5-10 Å: orange (moderate confidence)
    # >10 Å: red (low confidence)
    colors = []
    for pae in pae_values:
        if pae < 5:
            colors.append('blue')
        elif pae <= 10:
            colors.append('orange')
        else:
            colors.append('red')

    # Create bars with 70% opacity and error bars (SEM)
    bars = ax.bar(range(len(ligands)), pae_values, color=colors,
                   edgecolor='black', linewidth=0.5, alpha=0.7,
                   yerr=pae_sem_values, error_kw={'ecolor': 'black', 'capsize': 3, 'elinewidth': 1.5, 'zorder': 2})

    # Customize plot
    ax.set_ylabel('Ligand-Binding Pocket PAE (Å)', fontsize=14, fontweight='bold')
    ax.set_title('Ligand PAE Values for GAGA-0515 XGPRT Structures', fontsize=16, fontweight='bold')
    ax.set_xticks(range(len(ligands)))
    ax.set_xticklabels(ligands, rotation=45, ha='right', fontsize=12)
    ax.tick_params(axis='y', labelsize=12)

    # Set y-axis maximum to 18
    ax.set_ylim(0, 18)

    # Add grid for better readability
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)

    # Add value labels on top of bars (with offset to avoid overlap with error bars)
    for i, (bar, value, sem) in enumerate(zip(bars, pae_values, pae_sem_values)):
        height = bar.get_height()
        # Position label above the error bar (height + SEM + small offset)
        label_y = height + sem + 0.3
        ax.text(bar.get_x() + bar.get_width() / 2., label_y,
                f'{value:.1f}',
                ha='center', va='bottom', fontsize=10, zorder=10)

    # Create legend box in top left with color coding
    from matplotlib.patches import Rectangle
    legend_elements = [
        Rectangle((0, 0), 1, 1, fc='blue', alpha=0.7, edgecolor='black', label='<5 Å (High)'),
        Rectangle((0, 0), 1, 1, fc='orange', alpha=0.7, edgecolor='black', label='5-10 Å (Moderate)'),
        Rectangle((0, 0), 1, 1, fc='red', alpha=0.7, edgecolor='black', label='>10 Å (Low)')
    ]

    # Add legend with matching box style
    legend = ax.legend(handles=legend_elements, loc='upper left',
                      fontsize=11, frameon=True, fancybox=False, shadow=False,
                      edgecolor='black', facecolor='white')
    legend.get_frame().set_alpha(0.9)
    legend.get_frame().set_boxstyle('round')

    # Add text box with statistics directly below legend (moved lower to avoid overlap)
    stats_text = f'Mean PAE: {np.mean(pae_values):.2f} Å\nMedian PAE: {np.median(pae_values):.2f} Å'
    ax.text(0.02, 0.70, stats_text, transform=ax.transAxes,
            fontsize=11, verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor='black'))

    # Add group labels below x-axis
    # Calculate group boundaries
    group_positions = []
    current_pos = 0
    for group_name, ligands_in_group in group_order.items():
        # Count how many ligands from this group are actually present
        group_ligands = [lig for lig in ligands_in_group if lig in ligands]
        group_size = len(group_ligands)
        if group_size > 0:
            # Calculate center position for group label
            center_pos = current_pos + (group_size - 1) / 2
            group_positions.append({
                'name': group_name,
                'start': current_pos,
                'end': current_pos + group_size - 1,
                'center': center_pos,
                'size': group_size
            })
            current_pos += group_size

    # Add group labels and separators
    y_offset = -0.18  # Position below x-axis labels
    for i, group_info in enumerate(group_positions):
        # Add group label
        ax.text(group_info['center'], y_offset, group_info['name'],
                transform=ax.get_xaxis_transform(),
                ha='center', va='top', fontsize=10, fontweight='bold')

        # Add horizontal line above group label
        line_y = y_offset + 0.03
        ax.plot([group_info['start'] - 0.5, group_info['end'] + 0.5],
                [line_y, line_y],
                transform=ax.get_xaxis_transform(),
                color='black', linewidth=1.5, clip_on=False)

        # Add vertical separators between groups (except after last group)
        if i < len(group_positions) - 1:
            separator_x = group_info['end'] + 0.5
            ax.axvline(x=separator_x, color='black', linestyle='--',
                      linewidth=1, alpha=0.5, zorder=0)

    # Adjust layout to make room for group labels
    plt.subplots_adjust(bottom=0.18)

    # Save figure
    output_file = "ligand_pae_barplot.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nBar plot saved to: {output_file}")

    # Summary statistics
    print("\n" + "=" * 80)
    print("Summary Statistics:")
    print(f"  Number of structures: {len(ligand_pae_data)}")
    print(f"  Mean PAE: {np.mean(pae_values):.2f} Å")
    print(f"  Median PAE: {np.median(pae_values):.2f} Å")
    print(f"  Min PAE: {np.min(pae_values):.2f} Å ({ligands[np.argmin(pae_values)]})")
    print(f"  Max PAE: {np.max(pae_values):.2f} Å ({ligands[np.argmax(pae_values)]})")
    print("=" * 80)


if __name__ == "__main__":
    main()
