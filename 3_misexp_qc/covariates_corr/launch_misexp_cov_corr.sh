#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/misexp_cov_corr/misexp_cov_corr.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/misexp_cov_corr/misexp_cov_corr.out
#BSUB -q "normal"
#BSUB -R "select[mem>100000] rusage[mem=100000] span[hosts=1]"
#BSUB -M100000
#BSUB -n 1
#BSUB -J misexp_corr

# conda activate tv5_base

wkdir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3
script=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/INTERVAL-seq_rare_variants/4_misexpression_v3/3_misexp_qc/covariates_corr/misexp_cov_corr.py

output=${wkdir}/2_misexp_qc/misexp_cov_corr
mkdir -p $output

expression=${wkdir}/1_rna_seq_qc/tpm_mtx_inactive/tpm_4568samples_8779genes_inactive.tsv
covariates=${wkdir}/2_misexp_qc/covariates/covariates_xcell.tsv

python $script --name "interval" --expression $expression --covariates $covariates --output $output
