#!/bin/bash
#BSUB -o /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/logs/count_snp_indel_carriers_logs/chr%I_snp_indel_carriers_count.out
#BSUB -e /lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/logs/count_snp_indel_carriers_logs/chr%I_snp_indel_carriers_count.out
#BSUB -q "basement"
#BSUB -R "select[mem>50000] rusage[mem=50000] span[hosts=1]"
#BSUB -M50000
#BSUB -n 1
#BSUB -J intersect_vrnts[1-22]

# conda activate tv5_base
i=`expr ${LSB_JOBINDEX} - 1`
chroms=({1..22})
CHR=${chroms[$i]}

express=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/2_misexp_qc/misexp_gene_cov_corr/tpm_zscore_4568smpls_8610genes_flat_misexp_corr_qc.csv
vcf=/lustre/scratch126/humgen/projects/interval_wgs/final_release_freeze/gt_phased/interval_wgs.chr${CHR}.gt_phased.vcf.gz
gts_root=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/4_vrnt_enrich/snp_indel_count_carriers/vrnts_gts_intersect
paired_smpls=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/1_rna_seq_qc/wgs_rna_match/paired_wgs_rna_postqc_prioritise_wgs.tsv
genes_bed=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/4_vrnt_enrich/snp_indel_count_carriers/genes_bed/chr${CHR}_genes.bed
vrnts_bed=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/4_vrnt_enrich/snp_indel_count_carriers/vrnts_bed/chr${CHR}_vrnts.bed
out=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v2/4_vrnt_enrich/snp_indel_count_carriers/count_snp_indel_carriers

python snp_indel_carrier_count.py --express $express \
                                --chrom chr$CHR \
                                --vcf $vcf \
                                --gts_root $gts_root \
                                --paired_smpls $paired_smpls \
                                --genes_bed $genes_bed \
                                --vrnts_bed $vrnts_bed --out $out
