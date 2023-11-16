#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/intersect_vrnts/chr%I_intersect_vrnts.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/logs/intersect_vrnts/chr%I_intersect_vrnts.out
#BSUB -q "basement"
#BSUB -R "select[mem>50000] rusage[mem=50000] span[hosts=1]"
#BSUB -M50000
#BSUB -n 1
#BSUB -J intersect_vrnts[1-22]

# conda activate tv5_base
i=`expr ${LSB_JOBINDEX} - 1`
chroms=({1..22})
CHR=${chroms[$i]}

gene_bed=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/snp_indel_count_carriers/genes_bed/chr${CHR}_genes.bed
vrnts_bed=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/snp_indel_count_carriers/vrnts_bed/chr${CHR}_vrnts.bed
root=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3/4_vrnt_enrich/snp_indel_count_carriers/

python intersect_genes_vrnts.py --chrom chr$CHR --gene_bed $gene_bed --vrnts_bed $vrnts_bed --window 1000000 --root $root