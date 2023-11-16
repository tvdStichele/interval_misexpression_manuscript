#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/1_data_download/download_ref_data.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/1_data_download/download_ref_data.out
#BSUB -q "normal"
#BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
#BSUB -M10000
#BSUB -n 1
#BSUB -J misexp_ref_data_download

script=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/INTERVAL-seq_rare_variants/4_misexpression_v3/1_data_download/download_data.sh
collapse_script=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/INTERVAL-seq_rare_variants/4_misexpression_v3/1_data_download/collapse_annotation.py
out_dir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/reference

bash $script $collapse_script $out_dir