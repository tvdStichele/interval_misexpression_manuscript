#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/cadd_sv/bed_sort/sort_bed_%I.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/cadd_sv/bed_sort/sort_bed_%I.err
#BSUB -q "normal"
#BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
#BSUB -M10000
#BSUB -n 1
#BSUB -J sort_bed_files[1-25]

bed=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/5_misexp_vrnts/scores/cadd_sv/bed_paths_no_invs.txt
awk 'NR=='${LSB_JOBINDEX} $bed | while read -a line
do 
    bed_path=${line[0]}
    bed_root="$(cut -d'.' -f1 <<<"$bed_path")"
    sorted_bed_path="${bed_root}.sorted.bed"
    bedtools sort -i $bed_path > $sorted_bed_path
    rm $bed_path
    mv $sorted_bed_path $bed_path
done 