### 1. Alphafold3 JSON Input Generator
**`generate_af3_inputs.py`** - Complete AlphaFold3 JSON input generator


### 2. Analysis Scripts

**`analyze_af3_results_combined.py`** - Comprehensive analysis for all AF3 outputs
- Generates heatmaps with overall pLDDT values
- Extracts metrics: ipTM, pLDDT, PAE, contact probability

**`plot_ligand_pae_bars.py`** - Ligand:Contact residues PAE analysis for GAGA.
- Plotting and statistics for GAGA dimer predictions


**`extract_ligand_contacts.py`** - Ligand contact residue extraction
- Extracts residues within 5 Å (or custom distance) of all ligands
- Generates **single CSV file** with two columns: structure name and contact residues
- Residue format: `chain:number:type` (e.g., "A:37:ARG, B:53:ARG")
- Handles all ligand types: substrates, Mg²⁺, cofactors

### 3. Visualization
**`visualize_structures.py`** - Publication-ready PyMOL visualization
- Automatic ligand detection and coloring
- CMYK-safe colors
- Contact residue highlighting (5 Å threshold)
- Ray-tracing
