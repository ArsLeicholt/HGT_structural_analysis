#!/bin/bash
# Run combined analysis on all three output directories

python3 analyze_af3_results_combined.py \
    --results-dirs af3_output af3_output_2 af3_output_3 \
    --tree-file GAGA-0515_XGPRT.aln.treefile \
    --organism-map header_to_organism.txt \
    --output-prefix af3_combined

echo ""
echo "Analysis complete! Generated heatmaps:"
echo "  - af3_combined_output_heatmap.png (af3_output)"
echo "  - af3_combined_output_2_heatmap.png (af3_output_2)"
echo "  - af3_combined_output3_substrates_heatmap.png (af3_output_3)"
