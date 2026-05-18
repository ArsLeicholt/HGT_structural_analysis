Code for structural analysis in "Janina L Rinke, Lukas Franke, Ding He, Maike L Fischer, Joel Vizueta, Lars A Eicholt, Rasmus S Larsen, Zijun Xiong, Phoebe Cunningham, Lee M Henry, Martin Kaltenpoth, Jürgen Gadau, Guojie Zhang, Jacobus J Boomsma, Lukas Schrader, Comparative analysis of 163 ant genomes reveals recurrent horizontal gene transfer from bacteria to ants, GigaScience, 2026;, giag043, https://doi.org/10.1093/gigascience/giag043"

### 1. Alphafold3 JSON Input Generator
**`generate_af3_inputs.py`** - Complete AlphaFold3 JSON input generator

**Input JSONs (when using 8 sequences)**
- Set 1 (PRTases): 64 files (8 sequences × 8 scenarios)
- Set 2 (TCA): 96 files (8 sequences × 6 substrates × 2 assemblies)
- Set 3 (nucleotide): 64 files (8 sequences × 8 substrates)
- **Total**: 224 JSON input files

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
- Mesh surfaces around ligands
- Contact residue highlighting (5 Å threshold)
- Ray-tracing
