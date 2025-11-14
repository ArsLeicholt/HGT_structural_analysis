#!/bin/bash
# Extract ligand contacts from all AF3 output directories

echo "Extracting ligand contacts (5 Å threshold) from all predictions"
echo "================================================================"
echo ""

# Process af3_output
if [ -d "af3_output" ]; then
    echo "Processing af3_output..."
    python3 extract_ligand_contacts.py \
        --input af3_output \
        --output contacts_af3_output.csv \
        --distance 5.0
    echo ""
fi

# Process af3_output_2
if [ -d "af3_output_2" ]; then
    echo "Processing af3_output_2..."
    python3 extract_ligand_contacts.py \
        --input af3_output_2 \
        --output contacts_af3_output_2.csv \
        --distance 5.0
    echo ""
fi

# Process af3_output_3
if [ -d "af3_output_3" ]; then
    echo "Processing af3_output_3..."
    python3 extract_ligand_contacts.py \
        --input af3_output_3 \
        --output contacts_af3_output_3.csv \
        --distance 5.0
    echo ""
fi

echo "================================================================"
echo "All contact extraction complete!"
echo ""
echo "Output files:"
echo "  - contacts_af3_output.csv"
echo "  - contacts_af3_output_2.csv"
echo "  - contacts_af3_output_3.csv"
