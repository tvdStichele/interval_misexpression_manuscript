#!/bin/bash

star_fusion_rna_ids=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/6_misexp_dissect/star_fusion/star_fusion_rna_ids.txt
# Read the sample IDs from the text file line by line
while IFS= read -r sample_id; do
    echo "Processing sample ID: $sample_id"
    # launch STAR Fusion
    bash launch_star_fusion_max_sensitivity.sh ./star_fusion_max_sensitivity.sh $sample_id
done < $star_fusion_rna_ids