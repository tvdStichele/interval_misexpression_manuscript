#!/bin/bash
chrom=$1
gene_id_list=$2

wkdir=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/misexpression_v3
script=/lustre/scratch126/humgen/projects/interval_rna/interval_rna_seq/thomasVDS/lof_missense/INTERVAL-seq_rare_variants/4_misexpression_v3/5_vrnt_enrich/1_snp_indel_count_carrier/4_annotate_vrnts_by_gene/annotate_vrnts_carriers_only_by_gene.py
 
gene_id=`awk 'NR=='${LSB_JOBINDEX} $gene_id_list`
echo $gene_id

vcf=/lustre/scratch126/humgen/projects/interval_wgs/final_release_freeze/gt_phased/interval_wgs.chr${chrom}.gt_phased.vcf.gz
rna_ids=${wkdir}/1_rna_seq_qc/aberrant_smpl_qc/smpls_pass_qc_4568.csv
intersect_bed=${wkdir}/4_vrnt_enrich/snp_indel_count_carriers/intersect_beds/chr${chrom}_gene_vrnts_intersect.tsv
paired_smpls=${wkdir}/1_rna_seq_qc/wgs_rna_match/paired_wgs_rna_postqc_prioritise_wgs.tsv
root=${wkdir}/4_vrnt_enrich/snp_indel_count_carriers

python annotate_vrnts_carriers_only_by_gene.py --gene_id $gene_id \
                                               --chrom chr$chrom \
                                               --vcf $vcf \
                                               --rna_ids $rna_ids \
                                               --paired_smpls $paired_smpls \
                                               --intersect_bed $intersect_bed \
                                               --root $root