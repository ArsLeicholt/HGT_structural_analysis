#!/usr/bin/env pymol
"""
PyMOL visualization script for protein structures with different substrates.
CMYK-safe colors for publication figures.
"""

from pymol import cmd, stored
import glob
import os

# Configuration
STRUCTURE_DIR = "."
FILE_PATTERNS = ["*.pdb", "*.cif", "*.ent"]

# CMYK-safe ligand colors (no grey or white)
LIGAND_COLORS = [
    [0.0, 1.0, 1.0],   # cyan
    [1.0, 0.0, 1.0],   # magenta
    [1.0, 1.0, 0.0],   # yellow
    [0.0, 0.0, 1.0],   # blue
    [1.0, 0.5, 0.0],   # orange
    [0.0, 1.0, 0.0],   # green
    [0.0, 0.8, 0.8],   # dark cyan
    [0.8, 0.0, 0.8],   # dark magenta
    [1.0, 0.8, 0.0],   # gold
    [0.5, 0.0, 0.5],   # purple
]

LIGAND_TRANSPARENCY = 0.4
CONTACT_DISTANCE = 5.0
CHAIN_B_TRANSPARENCY = 0.5

# Protein colors (no grey or white)
COLOR_HELIX = "marine"
COLOR_SHEET = "wheat"
COLOR_LOOP = "lightblue"
COLOR_CONTACT_RESIDUES = "brightorange"
COLOR_MOTIF = "tv_red"
MOTIF_RESIDUES = "39-45"  # AVSRGGL motif region


def find_structure_files(directory, patterns):
    """Find all structure files"""
    files = []
    for pattern in patterns:
        search_path = os.path.join(directory, pattern)
        files.extend(glob.glob(search_path))
    files.sort()
    return files


def get_base_name(filepath):
    """Extract base name without extension"""
    basename = os.path.basename(filepath)
    name_without_ext = os.path.splitext(basename)[0]
    name_without_ext = name_without_ext.split('.')[0]
    return name_without_ext


def create_color(name, rgb):
    """Create custom PyMOL color"""
    cmd.set_color(name, rgb)


def detect_ligands(object_name):
    """Detect non-polymer ligands in structure"""
    ligand_selection = f"{object_name} and (organic or inorganic) and not (resn HOH+WAT+DOD)"
    stored.ligands = []
    cmd.iterate(ligand_selection, "stored.ligands.append(resn)")
    unique_ligands = list(set(stored.ligands))
    return unique_ligands


def main():
    """Main visualization function"""

    print("PyMOL Visualization Script")
    print("=" * 60)

    cmd.delete("all")
    cmd.bg_color("white")

    structure_files = find_structure_files(STRUCTURE_DIR, FILE_PATTERNS)

    if not structure_files:
        print("No structure files found")
        return

    print(f"Found {len(structure_files)} structures")

    # Load structures
    loaded_objects = []
    for i, filepath in enumerate(structure_files):
        obj_name = get_base_name(filepath)
        if obj_name in loaded_objects:
            obj_name = f"{obj_name}_{i}"
        cmd.load(filepath, obj_name)
        loaded_objects.append(obj_name)

    reference_object = loaded_objects[0]
    print(f"Reference: {reference_object}")

    # Align structures
    for obj_name in loaded_objects[1:]:
        results = cmd.align(f"{obj_name} and name CA",
                           f"{reference_object} and name CA")
        rmsd = results[0]
        print(f"Aligned {obj_name}: RMSD = {rmsd:.3f} Å")

    # Detect ligands
    all_ligands = {}
    ligand_color_map = {}
    color_index = 0

    for obj_name in loaded_objects:
        ligands = detect_ligands(obj_name)
        all_ligands[obj_name] = ligands

        for lig in ligands:
            if lig not in ligand_color_map:
                color_name = f"ligand_color_{lig}"
                rgb = LIGAND_COLORS[color_index % len(LIGAND_COLORS)]
                create_color(color_name, rgb)
                ligand_color_map[lig] = color_name
                color_index += 1

    # Style protein backbone
    for obj_name in loaded_objects:
        cmd.hide("everything", obj_name)
        cmd.show("cartoon", f"{obj_name} and polymer")

        # Color chain A (main chain)
        cmd.color(COLOR_HELIX, f"{obj_name} and chain A and ss h")
        cmd.color(COLOR_SHEET, f"{obj_name} and chain A and ss s")
        cmd.color(COLOR_LOOP, f"{obj_name} and chain A and ss l+''")

        # Color and make transparent chains B, C, D
        for chain in ['B', 'C', 'D']:
            cmd.color(COLOR_HELIX, f"{obj_name} and chain {chain} and ss h")
            cmd.color(COLOR_SHEET, f"{obj_name} and chain {chain} and ss s")
            cmd.color(COLOR_LOOP, f"{obj_name} and chain {chain} and ss l+''")
            cmd.set("cartoon_transparency", CHAIN_B_TRANSPARENCY, f"{obj_name} and chain {chain}")

        cmd.set("cartoon_sampling", 14, obj_name)
        cmd.set("cartoon_fancy_helices", 1, obj_name)
        cmd.set("cartoon_smooth_loops", 1, obj_name)

    # Style ligands
    for obj_name in loaded_objects:
        for ligand in all_ligands.get(obj_name, []):
            ligand_sel = f"{obj_name} and resn {ligand}"
            cmd.show("sticks", ligand_sel)
            color_name = ligand_color_map[ligand]
            cmd.color(color_name, ligand_sel)
            cmd.set("stick_transparency", LIGAND_TRANSPARENCY, ligand_sel)

    # Create mesh surfaces
    for obj_name in loaded_objects:
        for ligand in all_ligands.get(obj_name, []):
            ligand_sel = f"{obj_name} and resn {ligand}"
            mesh_name = f"mesh_{obj_name}_{ligand}"

            cmd.set("surface_quality", 2)
            cmd.set("surface_type", 0)
            cmd.create(mesh_name, ligand_sel)
            cmd.show("mesh", mesh_name)

            color_name = ligand_color_map[ligand]
            cmd.color(color_name, mesh_name)
            cmd.set("mesh_quality", 4, mesh_name)
            cmd.set("mesh_width", 0.3, mesh_name)

    # Highlight contact residues
    for obj_name in loaded_objects:
        if all_ligands.get(obj_name):
            ligand_list = " or ".join([f"resn {lig}" for lig in all_ligands[obj_name]])
            contact_sel = f"contact_residues_{obj_name}"
            cmd.select(contact_sel,
                      f"{obj_name} and polymer and (byres (({obj_name} and ({ligand_list})) around {CONTACT_DISTANCE}))")
            cmd.show("sticks", contact_sel)
            cmd.color(COLOR_CONTACT_RESIDUES, contact_sel)

            stored.count = 0
            cmd.iterate(f"{contact_sel} and name CA", "stored.count += 1")
            print(f"{obj_name}: {stored.count} contact residues")

            cmd.delete(contact_sel)

    # Highlight AVSRGGL motif in red (both chains A and B)
    for obj_name in loaded_objects:
        for chain in ['A', 'B']:
            motif_sel = f"motif_{obj_name}_{chain}"
            cmd.select(motif_sel, f"{obj_name} and chain {chain} and resi {MOTIF_RESIDUES}")
            cmd.show("sticks", motif_sel)
            cmd.color(COLOR_MOTIF, motif_sel)
            cmd.delete(motif_sel)

    # Publication settings
    cmd.set("ray_shadows", 1)
    cmd.set("ray_trace_mode", 1)
    cmd.set("ambient", 0.4)
    cmd.set("specular", 0.5)
    cmd.set("shininess", 10)
    cmd.set("depth_cue", 0)
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("stick_radius", 0.15)
    cmd.set("stick_quality", 15)
    cmd.set("sphere_quality", 3)

    cmd.orient("all")
    cmd.zoom("all", buffer=5)

    print("\nVisualization complete")
    print(f"Loaded objects: {', '.join(loaded_objects)}")
    print("\nLigand colors:")
    for ligand, color in sorted(ligand_color_map.items()):
        print(f"  {ligand}: {color}")


if __name__ == "__main__" or __name__ == "pymol":
    main()
