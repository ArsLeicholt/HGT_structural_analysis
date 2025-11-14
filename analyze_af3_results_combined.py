#!/usr/bin/env python3
"""
Analyze AlphaFold3 multimer prediction results from multiple output directories.
Supports af3_output, af3_output_2, and af3_output_3 with different substrate sets.

The heatmap displays:
- Y-axis: Organism names (ordered by phylogenetic tree)
- X-axis: Substrate + assembly combinations
- Each cell: 4 tiles showing normalized metrics (0-100%): ipTM, PTM, PAE, Ligand pLDDT
- Overall pLDDT displayed below x-axis labels
"""

import json
import os
import re
import csv
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from Bio import Phylo
from io import StringIO
import pandas as pd
from collections import defaultdict


# Substrate definitions for different output directories
SUBSTRATES_OUTPUT_1_2 = {
    "CIT": "Citrate",
    "MLT": "Malate",
    "OAA": "Oxaloacetate",
    "AKG": "2-Oxoglutarate",
    "GXO": "Glyoxylate",
    "HPV": "Hydroxypyruvate"
}

SUBSTRATES_OUTPUT_3 = {
    "AK": "Adenylate kinase",
    "NDK": "Nucleoside diphosphate kinase",
    "PK": "Pyruvate kinase",
    "PGK": "Phosphoglycerate kinase",
    "SCS": "Succinyl-CoA synthetase",
    "CS": "Citrate synthase",
    "PEPCK": "Phosphoenolpyruvate carboxykinase",
    "ME": "Malic enzyme"
}

ASSEMBLIES = ["dimer", "tetramer"]


def parse_organism_mapping(mapping_file):
    """
    Parse the header to organism mapping file.

    Returns:
        dict: Mapping from header to organism name
    """
    mapping = {}
    with open(mapping_file, 'r') as f:
        lines = f.readlines()

    # Skip header lines (markdown table format)
    for line in lines[2:]:  # Skip header and separator
        line = line.strip()
        if not line or line.startswith('|'):
            if '|' in line:
                parts = [p.strip() for p in line.split('|')]
                if len(parts) >= 3 and parts[1] and parts[2]:
                    header = parts[1].strip()
                    organism = parts[2].strip()
                    mapping[header] = organism

    return mapping


def parse_newick_tree(tree_file):
    """
    Parse Newick tree file and extract organism order based on tree topology.

    Returns:
        list: Ordered list of sequence headers (as they appear in tree)
    """
    tree = Phylo.read(tree_file, "newick")

    # Get terminal nodes (leaves) in order
    terminals = tree.get_terminals()
    ordered_headers = [term.name for term in terminals]

    return ordered_headers


def normalize_header(header):
    """
    Normalize header to match between tree and file names.
    """
    # Remove special characters that might have been changed in filenames
    normalized = header.replace('[', '_').replace(']', '_').replace(' ', '_')
    # Remove trailing dots/underscores
    normalized = normalized.rstrip('._')
    return normalized


def calculate_overall_plddt(cif_file):
    """
    Calculate overall average pLDDT for all protein atoms in the structure.

    Args:
        cif_file: Path to model.cif file

    Returns:
        float: Average pLDDT score (0-100), or None if no data
    """
    plddt_values = []

    with open(cif_file, 'r') as f:
        in_atom_section = False
        for line in f:
            # Check if we're entering the atom section
            if '_atom_site' in line:
                in_atom_section = True
                continue

            if in_atom_section and line.startswith('ATOM'):
                parts = line.split()
                if len(parts) < 16:
                    continue

                b_factor = float(parts[15])  # pLDDT stored as B-factor
                plddt_values.append(b_factor)

            elif in_atom_section and line.startswith('#'):
                # End of atom section
                break

    if not plddt_values:
        return None

    return np.mean(plddt_values)


def parse_cif_file(cif_file):
    """
    Parse CIF file to extract atom positions and B-factors (pLDDT).

    Returns:
        dict: Contains 'atoms' list with atom info, 'mg_position', and 'ligands'
    """
    atoms = []
    mg_position = None
    ligands = defaultdict(list)

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

                # Store MG position
                if comp_id == 'MG':
                    mg_position = np.array([x, y, z])

                # Store ligand atoms (non-protein)
                if record_type == 'HETATM' and comp_id not in ['MG', 'HOH']:
                    ligands[comp_id].append(atom_info)

            elif in_atom_section and line.startswith('#'):
                # End of atom section
                break

    return {
        'atoms': atoms,
        'mg_position': mg_position,
        'ligands': ligands
    }


def identify_ligand_chains(cif_file):
    """
    Identify which chains contain substrate ligands (excluding MG).

    Args:
        cif_file: Path to model.cif file

    Returns:
        set: Chain IDs that contain substrate ligands
    """
    data = parse_cif_file(cif_file)

    if not data['ligands']:
        return set()

    ligand_chains = set()

    # Get all ligand chains (ligands are already filtered to exclude MG and HOH)
    for ligand_name, ligand_atoms in data['ligands'].items():
        for atom in ligand_atoms:
            ligand_chains.add(atom['chain_id'])

    return ligand_chains


def calculate_ligand_plddt(cif_file, distance_threshold=4.0, min_contacts=3):
    """
    Calculate average pLDDT of binding pocket residues.
    Binding pocket = residues with ≥3 heavy atoms within 4 Å of substrate.

    Args:
        cif_file: Path to model.cif file
        distance_threshold: Distance threshold in Angstrom (default 4.0)
        min_contacts: Minimum number of heavy atoms within threshold (default 3)

    Returns:
        float: Average pLDDT score (0-100), or None if no data
    """
    data = parse_cif_file(cif_file)

    if not data['ligands']:
        print(f"  Warning: No ligands found in {cif_file}")
        return None

    # Get all substrate ligand atoms (already excludes MG and HOH)
    all_ligand_atoms = []
    for ligand_name, ligand_atoms in data['ligands'].items():
        all_ligand_atoms.extend(ligand_atoms)

    if not all_ligand_atoms:
        print(f"  Warning: No substrate ligands found in {cif_file}")
        return None

    # For each protein residue, count heavy atoms within distance of ligand
    residue_contacts = {}  # res_key -> list of (atom, plddt, is_within_threshold)
    residue_info = {}  # res_key -> list of atom info

    for atom in data['atoms']:
        if atom['record_type'] == 'ATOM' and atom['seq_id'] is not None:
            res_key = (atom['chain_id'], atom['seq_id'])

            # Count contacts for this atom
            is_contact = False
            for ligand_atom in all_ligand_atoms:
                dist = np.linalg.norm(atom['coords'] - ligand_atom['coords'])
                if dist <= distance_threshold:
                    is_contact = True
                    break

            # Store atom info
            if res_key not in residue_info:
                residue_info[res_key] = []
            residue_info[res_key].append({
                'atom_type': atom['atom_type'],
                'plddt': atom['plddt'],
                'is_contact': is_contact
            })

    # Filter to binding pocket residues (≥3 heavy atoms in contact)
    binding_pocket_plddts = []

    for res_key, atoms in residue_info.items():
        # Count heavy atoms (not hydrogen) in contact
        heavy_atom_contacts = sum(1 for a in atoms
                                   if a['is_contact'] and a['atom_type'] != 'H')

        if heavy_atom_contacts >= min_contacts:
            # This residue is in the binding pocket
            # Average pLDDT across all atoms in this residue
            avg_plddt = np.mean([a['plddt'] for a in atoms])
            binding_pocket_plddts.append(avg_plddt)

    if not binding_pocket_plddts:
        print(f"  Warning: No binding pocket residues found in {cif_file}")
        return None

    # Average across binding pocket residues
    return np.mean(binding_pocket_plddts)


def calculate_contact_probability(confidences_file, cif_file):
    """
    Calculate average contact probability between protein and substrate ligands.
    Uses top contacts (max of 10 or 20%) to capture strong binding interactions.

    Args:
        confidences_file: Path to confidences.json
        cif_file: Path to model.cif

    Returns:
        float: Average contact probability of top contacts (0-100), or None if no data
    """
    # Identify substrate ligand chains from CIF structure
    ligand_chains = identify_ligand_chains(cif_file)

    if not ligand_chains:
        return None

    # Load contact probabilities
    with open(confidences_file, 'r') as f:
        data = json.load(f)

    contact_probs = np.array(data['contact_probs'])
    token_chain_ids = data['token_chain_ids']

    # Identify protein and ligand token indices
    protein_tokens = []
    ligand_tokens = []

    for i, chain_id in enumerate(token_chain_ids):
        if chain_id in ligand_chains:
            ligand_tokens.append(i)
        elif chain_id not in ['E', 'F', 'G', 'H', 'I', 'J', 'K', 'L']:  # Heuristic: early letters are typically protein
            protein_tokens.append(i)

    if not ligand_tokens or not protein_tokens:
        return None

    # Extract all contact probabilities between protein and substrate ligands
    contact_values = []
    for prot_idx in protein_tokens:
        for lig_idx in ligand_tokens:
            prob = contact_probs[prot_idx, lig_idx]
            contact_values.append(prob)

    if not contact_values:
        return None

    # Sort contacts and take top max(10, 20%)
    contact_values = sorted(contact_values, reverse=True)
    n_top = max(10, int(len(contact_values) * 0.2))
    top_contacts = contact_values[:n_top]

    # Return mean of top contacts as percentage
    return np.mean(top_contacts) * 100


def calculate_ligand_pae(confidences_file, cif_file, distance_threshold=4.0, min_contacts=3):
    """
    Calculate average PAE between binding pocket residues and substrate ligand.
    Uses residue-level PAE for fine-grained analysis.
    Lower PAE = better, so we invert it to higher = better.

    Args:
        confidences_file: Path to confidences.json (full PAE matrix)
        cif_file: Path to model.cif
        distance_threshold: Distance threshold in Angstrom (default 4.0)
        min_contacts: Minimum number of heavy atoms within threshold (default 3)

    Returns:
        float: Inverted PAE score (0-100), or None if no data
    """
    data = parse_cif_file(cif_file)

    if not data['ligands']:
        return None

    # Get all substrate ligand atoms
    all_ligand_atoms = []
    for ligand_name, ligand_atoms in data['ligands'].items():
        all_ligand_atoms.extend(ligand_atoms)

    if not all_ligand_atoms:
        return None

    # Identify ligand chains
    ligand_chains = identify_ligand_chains(cif_file)
    if not ligand_chains:
        return None

    # Load full PAE matrix
    try:
        with open(confidences_file, 'r') as f:
            conf_data = json.load(f)

        pae_matrix = np.array(conf_data['pae'])
        token_chain_ids = conf_data['token_chain_ids']
        token_res_ids = conf_data.get('token_res_ids', [None] * len(token_chain_ids))
    except (KeyError, FileNotFoundError):
        return None

    # Identify binding pocket residues (same logic as pLDDT)
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
        return None

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
        return None

    # Extract PAE values between binding pocket and ligand
    pae_values = []
    for pocket_idx in pocket_token_indices:
        for lig_idx in ligand_token_indices:
            pae_val = pae_matrix[pocket_idx, lig_idx]
            if not np.isnan(pae_val):
                pae_values.append(pae_val)

    if not pae_values:
        return None

    # Invert PAE: lower is better, convert to 0-100 scale
    # Typical good PAE < 10, poor PAE > 20
    mean_pae = np.mean(pae_values)
    inverted_pae = max(0, 100 - (mean_pae * 100 / 20))

    return inverted_pae


def extract_af3_metrics(result_dir):
    """
    Extract metrics from AlphaFold3 output directory.

    Looks for ranking_scores.csv to find best model, then extracts metrics.

    Returns:
        dict: Metrics including iptm, ligand_plddt, contact_prob, ligand_pae, overall_plddt
    """
    metrics = {}
    result_path = Path(result_dir)

    # Find best model from ranking_scores.csv
    ranking_file = result_path / "ranking_scores.csv"
    best_sample_dir = None

    if ranking_file.exists():
        # Parse ranking scores to find best model
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
        sample_dirs = sorted([d for d in result_path.iterdir() if d.is_dir() and d.name.startswith('seed-')])
        if not sample_dirs:
            return {
                'iptm': None,
                'ligand_plddt': None,
                'contact_prob': None,
                'ligand_pae': None,
                'overall_plddt': None,
            }
        best_sample_dir = sample_dirs[0]

    # Look for summary file, full confidences, and CIF file
    summary_file = best_sample_dir / "summary_confidences.json"
    confidences_file = best_sample_dir / "confidences.json"
    cif_file = best_sample_dir / "model.cif"

    # Extract ipTM from summary
    if summary_file.exists():
        with open(summary_file, 'r') as f:
            data = json.load(f)

        # Extract metrics (normalize to 0-100)
        metrics['iptm'] = data.get('iptm', 0) * 100 if 'iptm' in data else None
    else:
        metrics['iptm'] = None

    # Calculate overall pLDDT and ligand-specific metrics
    if cif_file.exists():
        # 0. Overall pLDDT (all protein atoms)
        overall_plddt = calculate_overall_plddt(cif_file)
        metrics['overall_plddt'] = overall_plddt

        # 1. Ligand pLDDT (protein residues around substrate ligands)
        ligand_plddt = calculate_ligand_plddt(cif_file)
        metrics['ligand_plddt'] = ligand_plddt

        # 2. Contact probability (protein-substrate ligand contacts)
        if confidences_file.exists():
            contact_prob = calculate_contact_probability(confidences_file, cif_file)
            metrics['contact_prob'] = contact_prob
        else:
            metrics['contact_prob'] = None

        # 3. Ligand PAE (binding pocket-substrate ligand PAE, residue-level)
        if confidences_file.exists():
            ligand_pae = calculate_ligand_pae(confidences_file, cif_file)
            metrics['ligand_pae'] = ligand_pae
        else:
            metrics['ligand_pae'] = None
    else:
        metrics['overall_plddt'] = None
        metrics['ligand_plddt'] = None
        metrics['contact_prob'] = None
        metrics['ligand_pae'] = None

    return metrics


def parse_all_results(output_base_dir, organism_mapping, substrate_dict):
    """
    Parse all AlphaFold3 results and organize by organism, substrate, and assembly.

    Args:
        output_base_dir: Directory containing results
        organism_mapping: Header to organism name mapping
        substrate_dict: Dictionary of substrates to search for

    Returns:
        dict: Nested dict {organism: {scenario: metrics}}
              where scenario is "SUBSTRATE_ASSEMBLY" (e.g., "CIT_dimer")
    """
    results = {}

    # Find all result directories
    output_path = Path(output_base_dir)
    if not output_path.exists():
        print(f"Warning: Output directory {output_base_dir} does not exist")
        return results

    result_dirs = [d for d in output_path.iterdir() if d.is_dir()]

    print(f"Found {len(result_dirs)} result directories in {output_base_dir}")

    # Create regex pattern from substrates
    substrate_pattern = '|'.join(substrate_dict.keys())

    for result_dir in result_dirs:
        dir_name = result_dir.name

        # Parse directory name to extract header, substrate, and assembly
        # Format: {header}_{substrate}_{assembly}
        match = re.match(rf'(.+)_({substrate_pattern})_(dimer|tetramer)$', dir_name)

        if not match:
            print(f"Warning: Could not parse directory name: {dir_name}")
            continue

        header_part = match.group(1)
        substrate = match.group(2)
        assembly = match.group(3)
        scenario = f"{substrate}_{assembly}"

        # Find matching organism
        organism = None
        for orig_header, org_name in organism_mapping.items():
            normalized_header = normalize_header(orig_header)
            # Try various matching strategies
            if header_part.replace('_', ' ').lower() in orig_header.lower():
                organism = org_name
                break
            # Also try matching start of header
            if normalized_header[:30].lower() in header_part[:30].lower():
                organism = org_name
                break

        if not organism:
            print(f"Warning: Could not find organism for header: {header_part}")
            continue

        # Find the nested directory (lowercase version)
        nested_dirs = [d for d in result_dir.iterdir() if d.is_dir()]
        if not nested_dirs:
            print(f"Warning: No nested directory found in {dir_name}")
            continue

        # Use first nested directory
        nested_dir = nested_dirs[0]

        # Extract metrics
        metrics = extract_af3_metrics(nested_dir)

        # Store results
        if organism not in results:
            results[organism] = {}

        results[organism][scenario] = metrics

        print(f"Processed: {organism} - {scenario}")

    return results


def create_split_tile_heatmap(results, organism_order, organism_mapping, output_file, tree_file, substrates_dict, assemblies):
    """
    Create split-tile heatmap with 4 metrics per cell, with phylogenetic tree.
    Includes overall pLDDT values below x-axis labels.

    Each cell is divided into 4 quadrants:
    - Top-left: ipTM
    - Top-right: Ligand pLDDT
    - Bottom-left: Contact Probability
    - Bottom-right: Ligand PAE (inverted)
    """
    # Prepare data - create all substrate_assembly combinations
    scenarios_ordered = []
    for substrate in substrates_dict.keys():
        for assembly in assemblies:
            scenarios_ordered.append(f"{substrate}_{assembly}")

    # Map tree headers to organism names
    organisms = []
    for header in organism_order:
        normalized = normalize_header(header)
        for orig_header, org_name in organism_mapping.items():
            if normalize_header(orig_header).startswith(normalized[:20]):
                organisms.append(org_name)
                break

    # Remove duplicates while preserving order
    seen = set()
    organisms = [x for x in organisms if not (x in seen or seen.add(x))]

    n_organisms = len(organisms)
    n_scenarios = len(scenarios_ordered)

    # Calculate overall pLDDT for each scenario (average across organisms)
    scenario_plddts = {}
    for scenario in scenarios_ordered:
        plddt_values = []
        for organism in organisms:
            if organism in results and scenario in results[organism]:
                metrics = results[organism][scenario]
                if metrics.get('overall_plddt') is not None:
                    plddt_values.append(metrics['overall_plddt'])

        if plddt_values:
            scenario_plddts[scenario] = np.mean(plddt_values)
        else:
            scenario_plddts[scenario] = None

    # Create figure with extra space at bottom for labels and pLDDT
    fig, ax = plt.subplots(figsize=(20, 12))

    # Color map with harsher discrimination at higher values
    colors = ['#d73027', '#f46d43', '#fdae61', '#fee090', '#e0f3f8', '#abd9e9', '#74add1', '#4575b4', '#313695']
    positions = [0.0, 0.15, 0.30, 0.45, 0.60, 0.75, 0.85, 0.92, 1.0]
    cmap = LinearSegmentedColormap.from_list('custom', list(zip(positions, colors)))

    # Plot each cell as 4 sub-tiles
    for i, organism in enumerate(organisms):
        for j, scenario in enumerate(scenarios_ordered):

            if organism not in results or scenario not in results[organism]:
                continue

            metrics = results[organism][scenario]

            # Define the 4 sub-tiles (top-left, top-right, bottom-left, bottom-right)
            sub_tiles = [
                ('iptm', 0, 0.5),              # top-left
                ('ligand_plddt', 0.5, 0.5),    # top-right
                ('contact_prob', 0, 0),        # bottom-left
                ('ligand_pae', 0.5, 0)         # bottom-right
            ]

            for metric_name, x_offset, y_offset in sub_tiles:
                value = metrics.get(metric_name)

                if value is not None:
                    # Normalize to 0-1 for colormap
                    normalized_value = value / 100.0
                    color = cmap(normalized_value)
                else:
                    # Gray for missing data
                    color = '#cccccc'

                # Draw rectangle
                rect = mpatches.Rectangle(
                    (j + x_offset, n_organisms - i - 1 + y_offset),
                    0.5, 0.5,
                    facecolor=color,
                    edgecolor='white',
                    linewidth=1
                )
                ax.add_patch(rect)

    # Set axis labels
    ax.set_xlim(0, n_scenarios)
    ax.set_ylim(0, n_organisms)

    # X-axis: substrate_assembly with overall pLDDT
    # Map assembly names to mer notation
    assembly_map = {'dimer': '2mer', 'tetramer': '4mer'}

    x_labels = []
    for scenario in scenarios_ordered:
        substrate, assembly = scenario.split('_')
        mer_label = assembly_map.get(assembly, assembly)
        plddt = scenario_plddts.get(scenario)

        if plddt is not None:
            label = f"{substrate}\n({mer_label})\npLDDT: {plddt:.1f}"
        else:
            label = f"{substrate}\n({mer_label})\npLDDT: N/A"

        x_labels.append(label)

    # Add grid FIRST (minor ticks for white grid)
    ax.set_xticks(range(n_scenarios + 1), minor=True)
    ax.set_yticks(range(n_organisms + 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=3)

    # X-axis ticks and labels (set AFTER grid to not be overwritten)
    ax.set_xticks([i + 0.5 for i in range(n_scenarios)])
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=10, fontweight='bold')

    # Y-axis: organisms ordered by tree
    ax.set_yticks([i + 0.5 for i in range(n_organisms)])
    ax.set_yticklabels(organisms[::-1], style='italic', fontsize=10)

    # Title
    ax.set_title('AlphaFold3 Predictions - Substrate and Assembly Comparison',
                 fontsize=16, fontweight='bold', pad=20)

    # Colorbar legend
    sm = plt.cm.ScalarMappable(cmap=cmap,
                               norm=plt.Normalize(vmin=0, vmax=100))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02, aspect=30, fraction=0.046)
    cbar.set_label('Score (%)', rotation=270, labelpad=20)

    # Add metric explanations below x-axis
    explanation_y = -2.8  # Below x-axis labels and pLDDT

    ax.text(n_scenarios / 2, explanation_y,
            'Each cell shows 4 metrics in 2×2 tiles:',
            ha='center', fontsize=9, fontweight='bold', transform=ax.transData)

    ax.text(n_scenarios / 2, explanation_y - 0.35,
            'Top-left: ipTM (interface pTM) | Top-right: pLDDT (binding pocket: ≥3 heavy atoms within 4 Å)',
            ha='center', fontsize=8, style='italic', transform=ax.transData)

    ax.text(n_scenarios / 2, explanation_y - 0.65,
            'Bottom-left: Contact probability (top max(10, 20%)) | Bottom-right: Ligand PAE (binding pocket-ligand, residue-level, inverted)',
            ha='center', fontsize=8, style='italic', transform=ax.transData)

    # Save with proper spacing for labels
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', pad_inches=0.5)
    print(f"\nHeatmap saved to: {output_file}")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Analyze AlphaFold3 results from multiple output directories",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("--results-dirs", nargs='+', required=True,
                       help="Directory(ies) containing AlphaFold3 output subdirectories (e.g., af3_output af3_output_2 af3_output_3)")
    parser.add_argument("--tree-file", required=True,
                       help="Newick tree file for organism ordering")
    parser.add_argument("--organism-map", required=True,
                       help="File mapping FASTA headers to organism names")
    parser.add_argument("--output-prefix", default="af3_combined",
                       help="Output figure filename prefix (default: af3_combined)")

    args = parser.parse_args()

    print("AlphaFold3 Combined Results Analysis")
    print("=" * 80)

    # Parse organism mapping
    print("\n1. Parsing organism mapping...")
    organism_mapping = parse_organism_mapping(args.organism_map)
    print(f"   Found {len(organism_mapping)} organisms")

    # Parse phylogenetic tree
    print("\n2. Parsing phylogenetic tree...")
    organism_order = parse_newick_tree(args.tree_file)
    print(f"   Tree contains {len(organism_order)} taxa")

    # Parse all AlphaFold3 results from each directory
    print("\n3. Parsing AlphaFold3 results...")

    for results_dir in args.results_dirs:
        # Determine which substrate set to use
        if 'output_3' in results_dir or 'af3_output_3' in results_dir:
            substrates = SUBSTRATES_OUTPUT_3
            assemblies = ["dimer"]  # output_3 only has dimers
            output_suffix = "output3_substrates"
        else:
            substrates = SUBSTRATES_OUTPUT_1_2
            assemblies = ASSEMBLIES
            output_suffix = results_dir.replace('/', '_').replace('af3_', '')

        print(f"\n   Processing {results_dir}...")
        results = parse_all_results(results_dir, organism_mapping, substrates)
        print(f"   Parsed results for {len(results)} organisms")

        # Create visualization
        print(f"\n4. Creating split-tile heatmap for {results_dir}...")
        output_file = f"{args.output_prefix}_{output_suffix}_heatmap.png"
        create_split_tile_heatmap(results, organism_order, organism_mapping,
                                 output_file, args.tree_file, substrates, assemblies)

    print("\n" + "=" * 80)
    print("Analysis complete!")


if __name__ == "__main__":
    main()
