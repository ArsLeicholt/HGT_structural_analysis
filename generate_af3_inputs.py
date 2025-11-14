#!/usr/bin/env python3
"""
Generate AlphaFold3 JSON input files for different substrate scenarios.
Supports multimer, TCA, and nucleotide/phosphate transfer predictions.
"""

import json
import os
import argparse
import re
from pathlib import Path


# Set 1: Original multimer scenarios
MULTIMER_SCENARIOS = [
    {"appendix": "APRT", "ligands": ["ADE", "PRP", "MG"], "copies": 2,
     "description": "Homodimer with ADE, PRP, MG"},
    {"appendix": "HPRT", "ligands": ["HXP", "PRP", "MG"], "copies": 4,
     "description": "Homotetramer with HXP, PRP, MG"},
    {"appendix": "HGPRT", "ligands": ["GUA", "PRP", "MG"], "copies": 4,
     "description": "Homotetramer with GUA, PRP, MG"},
    {"appendix": "UPRT", "ligands": ["URA", "PRP", "MG"], "copies": 2,
     "description": "Homodimer with URA, PRP, MG"},
    {"appendix": "OPRT", "ligands": ["ORT", "PRP", "MG"], "copies": 4,
     "description": "Homotetramer with ORT, PRP, MG"},
    {"appendix": "XGPRT", "ligands": ["XAN", "PRP", "MG"], "copies": 4,
     "description": "Homotetramer with XAN, PRP, MG"},
    {"appendix": "PRPP", "ligands": ["PRP", "MG"], "copies": 1,
     "description": "Monomer with PRP, MG"},
    {"appendix": "PPDK", "ligands": ["PRP", "MG", "PYR", "ATP"], "copies": 4,
     "description": "Homotetramer with PRP, MG, PYR, ATP"}
]

# Set 2: TCA cycle substrates
TCA_SUBSTRATES = {
    "CIT": "Citrate", "MLT": "Malate", "OAA": "Oxaloacetate",
    "AKG": "2-Oxoglutarate", "GXO": "Glyoxylate", "HPV": "Hydroxypyruvate"
}

TCA_ASSEMBLIES = ["dimer", "tetramer"]

# Set 3: Nucleotide/phosphate transfer reactions
NUCLEOTIDE_SUBSTRATES = {
    "AK": ("Adenylate_kinase", ["ATP", "AMP", "MG"]),
    "NDK": ("Nucleoside_diphosphate_kinase", ["ATP", "GDP", "MG"]),
    "PK": ("Pyruvate_kinase", ["PEP", "ADP", "MG"]),
    "PGK": ("Phosphoglycerate_kinase", ["3PG", "ATP", "MG"]),
    "SCS": ("Succinyl_CoA_synthetase", ["SIN", "ATP", "COA", "MG"]),
    "CS": ("Citrate_synthase", ["OAA", "ACO"]),
    "PEPCK": ("Phosphoenolpyruvate_carboxykinase", ["OAA", "GTP", "MG"]),
    "ME": ("Malic_enzyme", ["MAL", "NAP", "MG"])
}

# Hardcoded sequences (can also be read from FASTA)
SEQUENCES = {
    "GAGA-0515.XGPRT_gene": "MTELSSKKYIVTWEMLQTHTRTLAKRLLSSTERWKGIIAVSRGGLVPAAILARELDIRHVDTVCISSYDHDVQRDLSVIKRAEGDGEGFIVVDDLVDTGVTAAAIRDLYPKAHFITIFAKPAGRPLVNDYVIDVPQDTWIDLPWDTGVAFVPPMAEVNRVH",
    "HBC80610.1.MAG.TPA.xanthine.phosphoribosyltransfer": "MSEKYIVTWDMLQIHARKLAARLMPSEQWKGIIAVSRGGLVPAALLARELGIRHVDTVCISSYDHDNQRELKVLKRAEGDGEGFIVIDDLVDTGGTAVAIREMYPKAHFVTIFAKPAGRPLVNDYVIDIPQDTWIEQPWDMGVVFVPPISGR",
    "WP_103776845.1.Citrobacter_amalonaticus": "MSEKYVVTWDMLQIHARKLASRLMPSEQWKGIIAVSRGGLVPGALLARELGIRHVDTVCISSYDHDNQRELKVLKRAEGDGEGFIVIDDLVDTGGTAVAIREMYPKAHFVTIFAKPAGRPLVNDYVIDIPQDTWIEQPWDMGVAFVPPISGR",
    "WP_095280890.1.Lelliottia_aquatilis": "MSEKYVVTWDMLQIHARKLAARLMPSEQWKGIIAVSRGGLVPGALLARELGIRHVDTVCISSYDHDNQRELTVLKRAEGDGEGFIVIDDLVDTGGTAVAIREMYPKAHFVTIFAKPAGRPLVNDYVIDIPQDTWIEQPWDMGVAFIPPISSR",
    "MBS5775149.1.Enterobacter_cloacae": "MSEKYVVTWDMLQIHARKLAARLMPSEQWKGIIAVSRGGLVPGALLARELGIRHVDTVCISSYDHDNQRELKVLKRAEGDGEGFIVIDDLVDTGGTAVAIREMYPKAHFVTIFAKPAGRPLVNDYVIDIPQDTWIEQPWDMGVAFVPPISGR",
    "MDP9770468.1.Atlantibacter_hermannii": "MSEKYVVTWDMLQIHARKLAARLMPSEQWKGIIAVSRGGLVPGALLARELGIRHVDTVCISSYDHDNQRELKVLKRAEGDGEGFIVIDDLVDTGGTAVAIREMYPKAHFVTIFAKPAGRPLVNDYVIDIPQDTWIEQPWDMGVAFVPPISGR",
    "MDV5355389.1.Enterobacter_asburiae": "MSEKYVVTWDMLQIHARKLAARLMPSEQWKGIIAVSRGGLVPGALLARELGIRHVDTVCISSYDHDNQRELKVLKRAEGDGEGFIVIDDLVDTGGTAVAIREMYPKAHFVTIFAKPAGRPLVNDYVIDIPQDTWIEQPWDMGVAFIPPISGR",
    "EKU4266144.1_xanthine_phosphoribosyltransferase_Ps": "MSEKYIVTWDMLQIHARKLASRLMPSEQWKGIIAVSRGGLVPGALLARELGIRHVDTVCISSYDHDNQRELKVLKRAEGDGEGFIVIDDLVDTGGTAVAIREMYPKAHFVTIFAKPAGRPLVDDYVVDIPQNTWIEQPWDMGVVFVPPISGR"
}


def parse_fasta(fasta_file):
    """Parse FASTA file and return dict of {header: sequence}"""
    sequences = {}
    current_header = None
    current_sequence = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_sequence)
                current_header = line[1:].strip()
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_header:
            sequences[current_header] = ''.join(current_sequence)

    return sequences


def sanitize_name(name):
    """Clean up header names for use in filenames"""
    sanitized = re.sub(r'[^\w\-_\.]', '_', name)
    sanitized = re.sub(r'_+', '_', sanitized)
    return sanitized.strip('_')


def generate_multimer_json(seq_name, sequence, scenario):
    """Generate JSON for original multimer scenarios"""
    chain_ids = ["A", "B", "C", "D", "E", "F", "G", "H"]

    protein_chains = []
    for i in range(scenario["copies"]):
        protein_chains.append({
            "protein": {"id": [chain_ids[i]], "sequence": sequence}
        })

    ligand_entries = []
    for i, ligand_code in enumerate(scenario["ligands"]):
        ligand_chain_id = chain_ids[scenario["copies"] + i]
        ligand_entries.append({
            "ligand": {"id": [ligand_chain_id], "ccdCodes": [ligand_code]}
        })

    return {
        "name": f"{seq_name}_{scenario['appendix']}",
        "sequences": protein_chains + ligand_entries,
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }


def generate_tca_json(seq_name, sequence, substrate_code, num_chains):
    """Generate JSON for TCA substrate scenarios"""
    chain_ids = ["A", "B", "C", "D", "E", "F", "G", "H"]

    protein_chains = []
    for i in range(num_chains):
        protein_chains.append({
            "protein": {"id": [chain_ids[i]], "sequence": sequence}
        })

    ligand_chain_id = chain_ids[num_chains]
    mg_chain_id = chain_ids[num_chains + 1]

    ligand_entry = {"ligand": {"id": [ligand_chain_id], "ccdCodes": [substrate_code]}}
    mg_entry = {"ligand": {"id": [mg_chain_id], "ccdCodes": ["MG"]}}

    stoich_name = "dimer" if num_chains == 2 else "tetramer"

    return {
        "name": f"{seq_name}_{substrate_code}_{stoich_name}",
        "sequences": protein_chains + [ligand_entry, mg_entry],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }


def generate_nucleotide_json(seq_name, sequence, substrate_code, substrate_info):
    """Generate JSON for nucleotide/phosphate transfer scenarios (dimers)"""
    chain_ids = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"]
    substrate_name, ligand_codes = substrate_info

    protein_chains = [
        {"protein": {"id": ["A"], "sequence": sequence}},
        {"protein": {"id": ["B"], "sequence": sequence}}
    ]

    ligand_entries = []
    for i, ligand_code in enumerate(ligand_codes):
        ligand_chain_id = chain_ids[2 + i]
        ligand_entries.append({
            "ligand": {"id": [ligand_chain_id], "ccdCodes": [ligand_code]}
        })

    return {
        "name": f"{seq_name}_{substrate_code}_dimer",
        "sequences": protein_chains + ligand_entries,
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }


def generate_multimer_set(sequences, output_dir):
    """Generate Set 1: Original multimer scenarios"""
    os.makedirs(output_dir, exist_ok=True)

    file_count = 0
    print(f"\nGenerating Set 1: Original Multimer Scenarios -> {output_dir}")

    for seq_name, sequence in sequences.items():
        for scenario in MULTIMER_SCENARIOS:
            json_data = generate_multimer_json(seq_name, sequence, scenario)
            filename = f"{seq_name}_{scenario['appendix']}.json"
            filepath = os.path.join(output_dir, filename)

            with open(filepath, 'w') as f:
                json.dump(json_data, f, indent=2)

            file_count += 1

    print(f"Created {file_count} files")
    return file_count


def generate_tca_set(sequences, output_dir):
    """Generate Set 2: TCA substrate scenarios"""
    os.makedirs(output_dir, exist_ok=True)

    file_count = 0
    print(f"\nGenerating Set 2: TCA Substrates -> {output_dir}")

    for seq_name, sequence in sequences.items():
        for substrate_code in TCA_SUBSTRATES.keys():
            for assembly in TCA_ASSEMBLIES:
                num_chains = 2 if assembly == "dimer" else 4
                json_data = generate_tca_json(seq_name, sequence, substrate_code, num_chains)
                filename = f"{seq_name}_{substrate_code}_{assembly}.json"
                filepath = os.path.join(output_dir, filename)

                with open(filepath, 'w') as f:
                    json.dump(json_data, f, indent=2)

                file_count += 1

    print(f"Created {file_count} files")
    return file_count


def generate_nucleotide_set(sequences, output_dir):
    """Generate Set 3: Nucleotide/phosphate transfer scenarios"""
    os.makedirs(output_dir, exist_ok=True)

    file_count = 0
    print(f"\nGenerating Set 3: Nucleotide/Phosphate Transfer -> {output_dir}")

    for seq_name, sequence in sequences.items():
        for substrate_code, substrate_info in NUCLEOTIDE_SUBSTRATES.items():
            json_data = generate_nucleotide_json(seq_name, sequence, substrate_code, substrate_info)
            filename = f"{seq_name}_{substrate_code}_dimer.json"
            filepath = os.path.join(output_dir, filename)

            with open(filepath, 'w') as f:
                json.dump(json_data, f, indent=2)

            file_count += 1

    print(f"Created {file_count} files")
    return file_count


def main():
    parser = argparse.ArgumentParser(description="Generate AlphaFold3 JSON input files")

    parser.add_argument("--fasta", "-f", help="Input FASTA file (optional, uses hardcoded sequences if not provided)")
    parser.add_argument("--all", action="store_true", help="Generate all three sets")
    parser.add_argument("--multimer", action="store_true", help="Generate Set 1: Original multimer scenarios")
    parser.add_argument("--tca", action="store_true", help="Generate Set 2: TCA substrates")
    parser.add_argument("--nucleotide", action="store_true", help="Generate Set 3: Nucleotide/phosphate transfer")

    args = parser.parse_args()

    if not (args.all or args.multimer or args.tca or args.nucleotide):
        parser.error("Specify at least one set: --all, --multimer, --tca, or --nucleotide")

    # Get sequences
    if args.fasta:
        print(f"Reading sequences from: {args.fasta}")
        sequences = parse_fasta(args.fasta)
    else:
        print("Using hardcoded sequences")
        sequences = SEQUENCES

    print(f"Processing {len(sequences)} sequences")

    total_files = 0

    if args.all or args.multimer:
        total_files += generate_multimer_set(sequences, "json_outputs")

    if args.all or args.tca:
        total_files += generate_tca_set(sequences, "json_outputs_2")

    if args.all or args.nucleotide:
        total_files += generate_nucleotide_set(sequences, "json_outputs_4")

    print(f"\nTotal: {total_files} JSON files created")


if __name__ == "__main__":
    main()
