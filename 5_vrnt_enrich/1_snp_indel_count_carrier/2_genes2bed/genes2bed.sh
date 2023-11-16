#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/genes2bed/genes2bed_%I.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/genes2bed/genes2bed_%I.out
#BSUB -q "normal"
#BSUB -R "select[mem>10000] rusage[mem=10000] span[hosts=1]"
#BSUB -M10000
#BSUB -n 1
#BSUB -J gene2bed[1-22]

# conda activate tv5_base
i=`expr ${LSB_JOBINDEX} - 1`
chroms=({1..22})
CHR=${chroms[$i]}

express=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/2_misexp_qc/misexp_gene_cov_corr/tpm_zscore_4568smpls_8610genes_flat_misexp_corr_qc.csv
gencode=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/reference/gencode/gencode.v31.annotation.sorted.gtf.gz
root=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/snp_indel_count_carriers

python genes2bed.py --express $express --chrom chr$CHR --gencode $gencode --root $root